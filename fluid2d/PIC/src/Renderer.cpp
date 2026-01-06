/**
 * Renderer.cpp: 2D PIC 烟雾渲染器实现
 * 将粒子散布到密度场，使用着色器渲染烟雾效果
 */

#include "PIC/include/Renderer.h"
#include <vector>
#include <cmath>
#include <glad/glad.h>
#include "Configure.h"
#include "Logger.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        Renderer::Renderer()
            : smokeShader(nullptr), quadVAO(0), quadVBO(0), FBO(0),
              textureID(0), RBO(0), densityTexture(0)
        {
            // 密度场分辨率：可独立于模拟网格，保持 FBO + 全屏 quad 管线不变
            // 适度提高分辨率能显著减少“像素块/糊”的观感（代价：CPU散布更慢）
            const int renderScale = 2;
            gridResX = PIC2dPara::theDim2d[0] * renderScale;
            gridResY = PIC2dPara::theDim2d[1] * renderScale;
            densityData.resize(gridResX * gridResY, 0.0f);

            initGLResources();
        }

        Renderer::~Renderer()
        {
            if (quadVBO) glDeleteBuffers(1, &quadVBO);
            if (quadVAO) glDeleteVertexArrays(1, &quadVAO);
            if (RBO) glDeleteRenderbuffers(1, &RBO);
            if (textureID) glDeleteTextures(1, &textureID);
            if (densityTexture) glDeleteTextures(1, &densityTexture);
            if (FBO) glDeleteFramebuffers(1, &FBO);
            delete smokeShader;
        }

        /**
         * 将粒子散布到密度场
         */
        void Renderer::updateDensityTexture(const ParticleSystem &ps)
        {
            // 模拟域尺寸（世界坐标）
            const float simH = PIC2dPara::theCellSize2d;
            const float domainW = PIC2dPara::theDim2d[0] * simH;
            const float domainH = PIC2dPara::theDim2d[1] * simH;

            // 密度网格单元大小（世界坐标）
            const float densHx = domainW / (float)gridResX;
            const float densHy = domainH / (float)gridResY;
            const float invHx = 1.0f / densHx;
            const float invHy = 1.0f / densHy;

            // 清空密度场
            std::fill(densityData.begin(), densityData.end(), 0.0f);

            // 二次 B-spline 核（3x3）散布：权重非负且和为 1（避免负权重导致“钳制后质量变大”）
            // x: fractional in [0,1)
            // 覆盖 i0-1, i0, i0+1
            auto bspline2 = [](float x) {
                float w0 = 0.5f * (1.0f - x) * (1.0f - x);
                float w1 = 0.75f - (x - 0.5f) * (x - 0.5f);
                float w2 = 0.5f * x * x;
                // 数值保险（理论上不会为负）
                w0 = (w0 < 0.0f) ? 0.0f : w0;
                w1 = (w1 < 0.0f) ? 0.0f : w1;
                w2 = (w2 < 0.0f) ? 0.0f : w2;
                return glm::vec3(w0, w1, w2);
            };

            // 将每个粒子的贡献散布到周围网格
            // 这相当于“每粒子密度质量”；太大容易整屏发白，太小则看不见。
            const float densityPerParticle = 0.15f;
            for (const auto &p : ps.particles)
            {
                // 计算粒子在网格中的位置
                float fx = p.position.x * invHx - 0.5f;
                float fy = p.position.y * invHy - 0.5f;

                int i0 = (int)std::floor(fx);
                int j0 = (int)std::floor(fy);

                float tx = fx - i0;
                float ty = fy - j0;

                glm::vec3 wx = bspline2(tx);
                glm::vec3 wy = bspline2(ty);

                // 覆盖 i0-1..i0+1, j0-1..j0+1
                for (int dj = -1; dj <= 1; ++dj)
                {
                    int j = j0 + dj;
                    if (j < 0 || j >= gridResY) continue;
                    float wj = wy[dj + 1];
                    for (int di = -1; di <= 1; ++di)
                    {
                        int i = i0 + di;
                        if (i < 0 || i >= gridResX) continue;
                        float wi = wx[di + 1];
                        float w = wi * wj;

                        densityData[j * gridResX + i] += w * densityPerParticle;
                    }
                }
            }

            // 限制密度范围并应用平滑
            for (auto &d : densityData)
            {
                // 保留一定动态范围，把“映射/对比度”交给 shader 做（指数吸收更自然）
                if (d < 0.0f) d = 0.0f;
                if (d > 10.0f) d = 10.0f;
            }

            // 上传到纹理
            glBindTexture(GL_TEXTURE_2D, densityTexture);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, gridResX, gridResY,
                            GL_RED, GL_FLOAT, densityData.data());
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        void Renderer::draw(const ParticleSystem &ps)
        {
            // 更新密度场纹理
            updateDensityTexture(ps);

            // 绑定帧缓冲
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);
            glViewport(0, 0, imageWidth, imageHeight);
            glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            if (!smokeShader)
            {
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                return;
            }

            // 使用烟雾着色器
            smokeShader->use();
            smokeShader->setFloat("contrast", PIC2dPara::contrast);

            // 绑定密度纹理
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, densityTexture);
            smokeShader->setInt("densityTex", 0);

            // 绘制全屏四边形
            glBindVertexArray(quadVAO);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
            glBindVertexArray(0);

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        GLuint Renderer::getTextureID() const
        {
            return textureID;
        }

        void Renderer::initGLResources()
        {
            extern std::string shaderPath;

            // 加载烟雾渲染着色器
            std::string vertPath = shaderPath + "/DrawSmoke2d.vert";
            std::string fragPath = shaderPath + "/DrawSmoke2d.frag";
            smokeShader = new Glb::Shader();
            smokeShader->buildFromFile(vertPath, fragPath);

            // 创建全屏四边形 VAO/VBO
            float quadVertices[] = {
                // 位置        // 纹理坐标
                -1.0f, -1.0f,  0.0f, 0.0f,
                 1.0f, -1.0f,  1.0f, 0.0f,
                -1.0f,  1.0f,  0.0f, 1.0f,
                 1.0f,  1.0f,  1.0f, 1.0f,
            };

            glGenVertexArrays(1, &quadVAO);
            glGenBuffers(1, &quadVBO);

            glBindVertexArray(quadVAO);
            glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)(2 * sizeof(float)));
            glBindVertexArray(0);

            // 创建密度场纹理
            glGenTextures(1, &densityTexture);
            glBindTexture(GL_TEXTURE_2D, densityTexture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, gridResX, gridResY, 0, GL_RED, GL_FLOAT, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glBindTexture(GL_TEXTURE_2D, 0);

            // 创建帧缓冲和渲染目标纹理
            glGenFramebuffers(1, &FBO);
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);

            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

            // 创建渲染缓冲
            glGenRenderbuffers(1, &RBO);
            glBindRenderbuffer(GL_RENDERBUFFER, RBO);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, imageWidth, imageHeight);
            glBindRenderbuffer(GL_RENDERBUFFER, 0);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);

            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            {
                Glb::Logger::getInstance().addLog("PIC Smoke Renderer framebuffer incomplete!");
            }

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }
    }
}

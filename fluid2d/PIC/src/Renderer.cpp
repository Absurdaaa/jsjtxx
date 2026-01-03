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
            // 使用网格分辨率作为密度场分辨率
            gridResX = PIC2dPara::theDim2d[0];
            gridResY = PIC2dPara::theDim2d[1];
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
            const float h = PIC2dPara::theCellSize2d;
            const float invH = 1.0f / h;

            // 清空密度场
            std::fill(densityData.begin(), densityData.end(), 0.0f);

            // 将每个粒子的贡献散布到周围网格
            for (const auto &p : ps.particles)
            {
                // 计算粒子在网格中的位置
                float fx = p.position.x * invH - 0.5f;
                float fy = p.position.y * invH - 0.5f;

                int i0 = (int)std::floor(fx);
                int j0 = (int)std::floor(fy);

                float tx = fx - i0;
                float ty = fy - j0;

                // 双线性插值权重散布
                for (int di = 0; di <= 1; ++di)
                {
                    for (int dj = 0; dj <= 1; ++dj)
                    {
                        int i = i0 + di;
                        int j = j0 + dj;

                        if (i < 0 || i >= gridResX || j < 0 || j >= gridResY)
                            continue;

                        float wx = di ? tx : (1.0f - tx);
                        float wy = dj ? ty : (1.0f - ty);
                        float w = wx * wy;

                        densityData[j * gridResX + i] += w * 0.5f;  // 每个粒子贡献的密度
                    }
                }
            }

            // 限制密度范围并应用平滑
            for (auto &d : densityData)
            {
                // d = std::min(d, 1.0f);
                if (d > 1.0f) d = 1.0f;
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

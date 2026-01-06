/**
 * Renderer.cpp: 2D PIC 烟雾渲染器实现
 * 使用与欧拉2D相同的方式：在每个像素位置插值采样密度
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

        void Renderer::draw(const ParticleSystem &ps, PICGrid2d &grid)
        {
            // 像欧拉2D一样：在每个像素位置采样密度并生成颜色
            std::vector<float> imageData;
            imageData.reserve(imageWidth * imageHeight * 3);

            const float maxX = grid.dim[PICGrid2d::X] * grid.cellSize;
            const float maxY = grid.dim[PICGrid2d::Y] * grid.cellSize;

            for (int j = 1; j <= imageHeight; j++)
            {
                for (int i = 1; i <= imageWidth; i++)
                {
                    float pt_x = i * maxX / imageWidth;
                    float pt_y = j * maxY / imageHeight;
                    glm::vec2 pt(pt_x, pt_y);

                    if (grid.inSolid(pt))
                    {
                        // 固体渲染为绿色
                        imageData.push_back(0.0f);
                        imageData.push_back(1.0f);
                        imageData.push_back(0.0f);
                    }
                    else
                    {
                        // 根据密度获取渲染颜色
                        glm::vec4 color = grid.getRenderColor(pt);
                        imageData.push_back(color.x * PIC2dPara::contrast);
                        imageData.push_back(color.y * PIC2dPara::contrast);
                        imageData.push_back(color.z * PIC2dPara::contrast);
                    }
                }
            }

            // 更新纹理
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, imageWidth, imageHeight,
                            GL_RGB, GL_FLOAT, imageData.data());
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        GLuint Renderer::getTextureID() const
        {
            return textureID;
        }

        void Renderer::initGLResources()
        {
            // 创建输出纹理（RGB浮点格式）
            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, imageWidth, imageHeight, 0, GL_RGB, GL_FLOAT, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    }
}

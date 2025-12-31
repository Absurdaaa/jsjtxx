
#include "PIC/include/Renderer.h"
#include <vector>
#include <glad/glad.h>
#include "Configure.h"
#include "Logger.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        // 构造函数：初始化所有OpenGL资源
        // shader: 着色器对象
        // VAO: 顶点数组对象
        // VBO: 顶点缓冲对象
        // FBO: 帧缓冲对象
        // textureID: 渲染结果纹理
        // RBO: 渲染缓冲对象（深度/模板）
        Renderer::Renderer()
            : shader(nullptr), VAO(0), VBO(0), FBO(0), textureID(0), RBO(0)
        {
            initGLResources(); // 初始化OpenGL资源
        }

        // 析构函数：释放所有OpenGL资源
        Renderer::~Renderer()
        {
            if (VBO) glDeleteBuffers(1, &VBO); // 删除顶点缓冲
            if (VAO) glDeleteVertexArrays(1, &VAO); // 删除顶点数组
            if (RBO) glDeleteRenderbuffers(1, &RBO); // 删除渲染缓冲
            if (textureID) glDeleteTextures(1, &textureID); // 删除纹理
            if (FBO) glDeleteFramebuffers(1, &FBO); // 删除帧缓冲
            delete shader; shader = nullptr; // 删除着色器
        }

        // 绘制所有粒子到帧缓冲纹理
        // ps: 粒子系统，包含所有粒子位置
        void Renderer::draw(const ParticleSystem &ps)
        {
            // 绑定帧缓冲，设置视口和清空背景
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);
            glViewport(0, 0, imageWidth, imageHeight);
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // 没有粒子或着色器未初始化则直接返回
            if (!shader || ps.particles.empty())
            {
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                return;
            }

            // 收集所有粒子位置到数组
            std::vector<float> positions;
            positions.reserve(ps.particles.size() * 2);
            for (const auto &p : ps.particles)
            {
                positions.push_back(p.position.x); // x坐标
                positions.push_back(p.position.y); // y坐标
            }

            // 绑定VAO/VBO并上传数据
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_DYNAMIC_DRAW);

            // 激活着色器并设置uniform参数
            shader->use();
            float domainX = Eulerian2dPara::theDim2d[0] * Eulerian2dPara::theCellSize2d; // 域宽
            float domainY = Eulerian2dPara::theDim2d[1] * Eulerian2dPara::theCellSize2d; // 域高
            shader->setFloat("domainX", domainX);
            shader->setFloat("domainY", domainY);
            shader->setFloat("pointSize", 3.0f); // 点大小

            // 绘制所有粒子为GL_POINTS
            glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(ps.particles.size()));

            // 解绑VAO和帧缓冲
            glBindVertexArray(0);
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        // 获取渲染结果的纹理ID
        GLuint Renderer::getTextureID() const
        {
            return textureID;
        }

        // 初始化所有OpenGL资源（着色器、VAO/VBO、FBO/纹理、RBO）
        void Renderer::initGLResources()
        {
            extern std::string shaderPath;

            // 加载粒子渲染着色器
            std::string vertPath = shaderPath + "/DrawPICParticles2d.vert";
            std::string fragPath = shaderPath + "/DrawPICParticles2d.frag";
            shader = new Glb::Shader();
            shader->buildFromFile(vertPath, fragPath);

            // 创建VAO和VBO
            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);

            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            // 位置属性：2D float
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
            glBindVertexArray(0);

            // 创建帧缓冲和纹理
            glGenFramebuffers(1, &FBO);
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);

            // 绑定纹理到帧缓冲
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

            // 创建并绑定渲染缓冲（深度/模板）
            glGenRenderbuffers(1, &RBO);
            glBindRenderbuffer(GL_RENDERBUFFER, RBO);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, imageWidth, imageHeight);
            glBindRenderbuffer(GL_RENDERBUFFER, 0);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);

            // 检查帧缓冲完整性
            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            {
                Glb::Logger::getInstance().addLog("PIC Renderer framebuffer incomplete!");
            }

            glBindFramebuffer(GL_FRAMEBUFFER, 0);

            // 启用OpenGL点大小和禁用深度测试
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_PROGRAM_POINT_SIZE);
        }
    }
}

/**
 * Renderer3d.cpp: 3D PIC 渲染器实现
 * 实现粒子系统的 OpenGL 渲染功能
 */

#include "PIC/include/Renderer3d.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 初始化渲染器
         */
        void Renderer3d::init()
        {
            // 创建容器边框
            container = new Glb::Container();
            float h = Eulerian3dPara::theCellSize3d;
            int nx = Eulerian3dPara::theDim3d[0];
            int ny = Eulerian3dPara::theDim3d[1];
            int nz = Eulerian3dPara::theDim3d[2];
            container->resetSize(nx * h, ny * h, nz * h);
            container->init();

            // 创建着色器
            shader = new Glb::Shader();
            std::string vertPath = shaderPath + "/DrawParticles3d.vert";
            std::string fragPath = shaderPath + "/DrawParticles3d.frag";
            shader->buildFromFile(vertPath, fragPath);

            // 创建帧缓冲对象
            glGenFramebuffers(1, &FBO);
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            // 创建渲染目标纹理
            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);

            // 将纹理附加到帧缓冲
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

            // 创建渲染缓冲对象（深度+模板）
            glGenRenderbuffers(1, &RBO);
            glBindRenderbuffer(GL_RENDERBUFFER, RBO);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, imageWidth, imageHeight);
            glBindRenderbuffer(GL_RENDERBUFFER, 0);

            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
            glBindFramebuffer(GL_FRAMEBUFFER, 0);

            // 创建顶点缓冲
            glGenBuffers(1, &VBO);
            MakeVertexArrays();

            glEnable(GL_MULTISAMPLE);
            glViewport(0, 0, imageWidth, imageHeight);
        }

        /**
         * 创建顶点数组对象
         * 设置粒子数据的顶点属性布局
         */
        void Renderer3d::MakeVertexArrays()
        {
            glGenVertexArrays(1, &VAO);
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);

            // 属性 0：位置（vec3）
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void *)offsetof(Particle, position));
            glEnableVertexAttribArray(0);

            // 属性 1：速度（vec3）- 可用于着色
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void *)offsetof(Particle, velocity));
            glEnableVertexAttribArray(1);

            glBindVertexArray(0);
        }

        /**
         * 绘制粒子系统
         */
        void Renderer3d::draw(ParticleSystem3d &ps)
        {
            // 上传粒子数据到 GPU
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, VBO);
            glBufferData(GL_SHADER_STORAGE_BUFFER, ps.particles.size() * sizeof(Particle), ps.particles.data(), GL_DYNAMIC_COPY);
            particleNum = ps.particles.size();

            // 绑定帧缓冲并清空
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);
            glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glEnable(GL_PROGRAM_POINT_SIZE);

            // 设置着色器参数
            shader->use();
            shader->setMat4("view", Glb::Camera::getInstance().GetView());
            shader->setMat4("projection", Glb::Camera::getInstance().GetProjection());
            shader->setFloat("scale", 1.0f);

            // 绘制粒子（点精灵）
            glBindVertexArray(VAO);
            glDrawArrays(GL_POINTS, 0, particleNum);
            shader->unUse();

            // 绘制容器边框
            container->draw();

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        /**
         * 获取渲染结果纹理 ID
         */
        GLuint Renderer3d::getRenderedTexture()
        {
            return textureID;
        }
    }
}

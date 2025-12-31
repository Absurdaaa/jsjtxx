#pragma once
#ifndef __PIC_RENDERER_2D_H__
#define __PIC_RENDERER_2D_H__

#include <glad/glad.h>
#include "ParticleSystem.h"
#include "Shader.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        /**
         * 2D PIC 粒子渲染器
         * 负责将粒子系统渲染到帧缓冲纹理，供UI显示
         */
        class Renderer
        {
        public:
            /**
             * 构造函数：初始化OpenGL资源
             */
            Renderer();

            /**
             * 析构函数：释放所有OpenGL资源
             */
            ~Renderer();

            /**
             * 绘制所有粒子到帧缓冲纹理
             * @param ps 粒子系统，包含所有粒子位置
             */
            void draw(const ParticleSystem &ps);

            /**
             * 获取渲染结果的纹理ID
             * @return OpenGL 2D纹理ID
             */
            GLuint getTextureID() const;

        private:
            /**
             * 初始化所有OpenGL资源（着色器、VAO/VBO、FBO/纹理、RBO）
             */
            void initGLResources();

            Glb::Shader* shader;   ///< 着色器对象，负责粒子点渲染
            GLuint VAO;            ///< 顶点数组对象
            GLuint VBO;            ///< 顶点缓冲对象
            GLuint FBO;            ///< 帧缓冲对象
            GLuint textureID;      ///< 渲染结果纹理ID
            GLuint RBO;            ///< 渲染缓冲对象（深度/模板）
        };
    }
}

#endif

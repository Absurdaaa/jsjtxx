/**
 * Renderer3d.h: 3D PIC 渲染器头文件
 * 定义粒子系统的 OpenGL 渲染类
 */

#pragma once
#ifndef __PIC_3D_RENDERER_H__
#define __PIC_3D_RENDERER_H__

#include <glm/glm.hpp>
#include <glad/glad.h>
#include "glfw3.h"
#include "Shader.h"
#include "Container.h"
#include "Camera.h"
#include "Global.h"
#include "Configure.h"
#include "ParticleSystem3d.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 3D PIC 渲染器类
         * 负责将粒子系统渲染到帧缓冲
         */
        class Renderer3d
        {
        public:
            Renderer3d() {}

            /**
             * 初始化渲染器
             * 创建着色器、帧缓冲、顶点数组等 OpenGL 资源
             */
            void init();

            /**
             * 获取渲染结果纹理 ID
             */
            GLuint getRenderedTexture();

            /**
             * 绘制粒子系统
             * @param ps 粒子系统引用
             */
            void draw(ParticleSystem3d &ps);

        private:
            /**
             * 创建顶点数组对象
             */
            void MakeVertexArrays();

            Glb::Shader *shader = nullptr;      // 着色器：用于渲染粒子（顶点/片段着色器）
            Glb::Container *container = nullptr; // 容器边框：可视化仿真域边界

            GLuint FBO = 0;        // 帧缓冲对象：渲染目标的 FBO
            GLuint RBO = 0;        // 渲染缓冲对象：深度/模板缓冲
            GLuint VAO = 0;        // 顶点数组对象：描述顶点属性布局
            GLuint VBO = 0;        // 顶点缓冲对象：存放粒子数据的 GPU 缓冲
            GLuint textureID = 0;  // 渲染目标纹理：最终渲染输出的贴图 ID

            int32_t particleNum = 0;  // 粒子数量：当前上传到 GPU 的粒子数
        };
    }
}

#endif

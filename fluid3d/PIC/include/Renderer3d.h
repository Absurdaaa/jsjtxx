/**
 * Renderer3d.h: 3D PIC 渲染器头文件
 * 仅保留 Volume 体渲染：将密度上传为 3D texture，并在片元着色器中做 ray marching。
 */

#pragma once
#ifndef __PIC_3D_RENDERER_H__
#define __PIC_3D_RENDERER_H__

#include <vector>

#include <glm/glm.hpp>
#include <glad/glad.h>

#include "Configure.h"
#include "Shader.h"
#include "Container.h"
#include "Camera.h"

#include "ParticleSystem3d.h"
#include "PICGrid3d.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 3D PIC 渲染器类
         * 输出：一张 2D 纹理（FBO color attachment），SceneView 直接显示。
         */
        class Renderer3d
        {
        public:
            Renderer3d();
            ~Renderer3d();

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
            void draw(const ParticleSystem3d &ps, PICGrid3d &grid);

        private:
            void initGLResources();

            // --- volume rendering (true 3D) ---
            void ensureDensityTexture3D(PICGrid3d &grid);
            void uploadDensityTexture3D(PICGrid3d &grid);

            // 渲染输出（FBO -> textureID），SceneView 直接显示这张纹理
            GLuint FBO = 0;
            GLuint RBO = 0;
            GLuint textureID = 0;

            // 体渲染：3D 密度纹理 (R32F)
            GLuint densityTex3D = 0;
            int densityNx = 0;
            int densityNy = 0;
            int densityNz = 0;
            std::vector<float> densityData;

            // 体渲染：全屏三角形
            GLuint volumeVAO = 0;
            Glb::Shader *volumeShader = nullptr;
            Glb::Container *container = nullptr;

            bool inited = false;
        };
    }
}

#endif

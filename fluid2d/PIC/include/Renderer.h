/**
 * Renderer.h: 2D PIC 烟雾渲染器头文件
 * 使用密度场纹理渲染烟雾效果
 */

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
         * 2D PIC 烟雾渲染器
         * 将粒子散布到密度场，然后渲染为烟雾效果
         */
        class Renderer
        {
        public:
            Renderer();
            ~Renderer();

            /**
             * 绘制烟雾效果
             * @param ps 粒子系统
             */
            void draw(const ParticleSystem &ps);

            GLuint getTextureID() const;

        private:
            void initGLResources();
            void updateDensityTexture(const ParticleSystem &ps);

            Glb::Shader *smokeShader;  // 烟雾渲染着���器

            GLuint quadVAO, quadVBO;   // 全屏四边形
            GLuint FBO;                // 帧缓冲
            GLuint textureID;          // 渲染结果纹理
            GLuint RBO;                // 深度/模板缓冲
            GLuint densityTexture;     // 密度场纹理

            int gridResX, gridResY;    // 密度场分辨率
            std::vector<float> densityData;  // 密度场数据
        };
    }
}

#endif

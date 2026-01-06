/**
 * Renderer.h: 2D PIC 烟雾渲染器头文件
 * 使用与欧拉2D相同的方式渲染密度场
 */

#pragma once
#ifndef __PIC_RENDERER_2D_H__
#define __PIC_RENDERER_2D_H__

#include <glad/glad.h>
#include "ParticleSystem.h"
#include "PICGrid2d.h"
#include "Shader.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        class Renderer
        {
        public:
            Renderer();
            ~Renderer();

            void draw(const ParticleSystem &ps, PICGrid2d &grid);
            GLuint getTextureID() const;

        private:
            void initGLResources();

            // 不再需要的成员保留为0以便析构函数安全删除
            Glb::Shader *smokeShader;
            GLuint quadVAO, quadVBO;
            GLuint FBO;
            GLuint textureID;
            GLuint RBO;
            GLuint densityTexture;

            int gridResX, gridResY;
            std::vector<float> densityData;
        };
    }
}

#endif

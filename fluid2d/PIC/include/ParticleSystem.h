#pragma once
#ifndef __PIC_PARTICLE_SYSTEM_2D_H__
#define __PIC_PARTICLE_SYSTEM_2D_H__

#include <vector>
#include <glm/glm.hpp>
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        /**
         * 单个粒子结构体
         * position: 粒子世界坐标
         * velocity: 粒子速度
         */
        struct Particle
        {
            glm::vec2 position{0.0f}; ///< 粒子位置
            glm::vec2 velocity{0.0f}; ///< 粒子速度
        };

        /**
         * 2D PIC 粒子系统
         * 管理所有粒子的存储与发射
         */
        class ParticleSystem
        {
        public:
            std::vector<Particle> particles; ///< 所有粒子

            /**
             * 根据配置的源添加粒子（每步调用，类似持续发射）
             * 粒子初始位置和速度由源参数决定
             */
            void emitFromSources();

            /**
             * 将所有粒子限制在容器边界内
             * @param lower 边界左下角
             * @param upper 边界右上角
             */
            void clampToDomain(const glm::vec2 &lower, const glm::vec2 &upper);
        };
    }
}

#endif

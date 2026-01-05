/**
 * ParticleSystem3d.h: 3D PIC 粒子系统头文件
 * 定义粒子结构体和粒子系统类
 */

#pragma once
#ifndef __PIC_PARTICLE_SYSTEM_3D_H__
#define __PIC_PARTICLE_SYSTEM_3D_H__

#include <vector>
#include <glm/glm.hpp>
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 粒子结构体
         * 存储粒子的位置和速度
         * 成员变量：
         * - position: 粒子在世界空间中的位置 (x,y,z)
         * - velocity: 粒子的速度向量 (x,y,z)
         */
        struct Particle
        {
            glm::vec3 position{0.0f};  // 粒子位置
            glm::vec3 velocity{0.0f};  // 粒子速度
        };

        /**
         * 3D 粒子系统类
         * 管理所有粒子的发射和边界约束
         */
        class ParticleSystem3d
        {
        public:
            std::vector<Particle> particles;  // 粒子容器：保存当前活跃的粒子列表

            /**
             * 从烟雾源发射新粒子
             */
            void emitFromSources();

            /**
             * 将粒子位置限制在指定域内
             * @param lower 域的下边界
             * @param upper 域的上边界
             */
            void clampToDomain(const glm::vec3 &lower, const glm::vec3 &upper);
        };
    }
}

#endif

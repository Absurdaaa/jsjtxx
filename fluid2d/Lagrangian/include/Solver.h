#pragma once
#ifndef __LAGRANGIAN_2D_SOLVER_H__
#define __LAGRANGIAN_2D_SOLVER_H__

#include "ParticleSystem2d.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace Lagrangian2d
    {
        // 拉格朗日法流体求解器类
        // 负责求解基于粒子的流体动力学方程
        class Solver
        {
        public:
            float velocityAttenuation_; // 碰撞后速度衰减系数
            float viscosity_;          // 粘度系数
            int substeps_;            // 子步数
            float dt_;                // 时间步长
            float density_;          // 流体密度
            float stiffness_;       // 刚度
            float exponent_;         // 压力计算公式中的指数
            float gravityX_;        // x轴重力
            float gravityY_;        // y轴重力
            float maxVelocity_;     // 最大允许速度
            Solver(ParticleSystem2d &ps);

            void solve();                         // 求解流体方程

        private:
            ParticleSystem2d &mPs;               // 粒子系统引用
            
            float cul_density(ParticleInfo2d& pi);    // 计算粒子密度
            float cul_pressure(ParticleInfo2d& pi);   // 计算粒子压力

            glm::vec2 cul_acceleration(ParticleInfo2d& pi); // 计算粒子zong加速度
            
            void update_velocity_position(ParticleInfo2d& pi, float dt); // 更新粒子速度和位置
            void boundary_check(ParticleInfo2d& pi); // 边界检查
            void update_blockId(ParticleInfo2d& pi); // 更新粒子块ID
            
            float PI = 3.14159265358979323846f;
            // Poly6核函数（2D，密度计算专用）
            float poly6Kernel(float dist2, float supportRadius, float supportRadius2);
            float SpikyKernel(float dist, float supportRadius);

           
        };
    }
}

#endif

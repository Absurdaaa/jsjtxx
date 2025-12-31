#pragma once
#ifndef __PIC_SOLVER_2D_H__
#define __PIC_SOLVER_2D_H__

#include "ParticleSystem.h"
#include "PICGrid2d.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        /**
         * 2D PIC 求解器：粒子-网格混合
         * 负责粒子与网格的相互转换、力学计算与推进
         */
        class Solver
        {
        public:
            /**
             * 构造函数
             * @param ps 粒子系统
             * @param grid 网格对象
             */
            Solver(ParticleSystem &ps, PICGrid2d &grid);

            /**
             * 执行一次完整的仿真推进（发射、P2G、加力、投影、G2P、推进）
             */
            void solve();

        private:
            ParticleSystem &mPs; ///< 引用的粒子系统
            PICGrid2d &mGrid;    ///< 引用的网格对象

            Glb::GridData2dX mU_prev;
            Glb::GridData2dY mV_prev;

            /**
             * 粒子速度投射到网格（P2G）
             */
            void particleToGrid();

            /**
             * 给网格加外力（如重力）
             * @param dt 时间步长
             */
            void addForces(double dt);

            /**
             * 压力投影，保证不可压缩
             * @param dt 时间步长
             */
            void pressureProjection(double dt);

            /**
             * 网格速度插值回粒子（G2P）
             * @param dt 时间步长
             */
            void gridToParticle(double dt);

            /**
             * 推进所有粒子位置
             * @param dt 时间步长
             */
            void advectParticles(double dt);
        };
    }
}

#endif

/**
 * Solver3d.h: 3D PIC 求解器头文件
 * 定义 PIC/FLIP 混合方法的求解器类
 */

#pragma once
#ifndef __PIC_SOLVER_3D_H__
#define __PIC_SOLVER_3D_H__

#include "ParticleSystem3d.h"
#include "PICGrid3d.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 3D PIC 求解器类
         * 实现粒子-网格混合方法的完整求解流程
         */
        class Solver3d
        {
        public:
            /**
             * 构造函数
             * @param ps 粒子系统引用
             * @param grid 网格引用
             */
            Solver3d(ParticleSystem3d &ps, PICGrid3d &grid);

            /**
             * 执行一次完整的仿真步进
             * 流程：发射粒子 -> 更新源 -> 标量平流 -> P2G -> 加力 -> 压力投影 -> G2P -> 粒子移动
             */
            void solve();

        private:
            ParticleSystem3d &mPs;   // 粒子系统引用：操作与查询粒子数据
            PICGrid3d &mGrid;        // 网格引用：访问/修改网格速度、密度、温度等字段

            // 保存上一步的速度场，用于 FLIP 计算（计算网格速度变化量）
            Glb::GridData3dX mU_prev;  // 上一步 X 方向速度（用于 FLIP 差值）
            Glb::GridData3dY mV_prev;  // 上一步 Y 方向速度
            Glb::GridData3dZ mW_prev;  // 上一步 Z 方向速度

            void particleToGrid();           // P2G：粒子速度转移到网格
            void addForces(double dt);       // 添加外力（浮力）
            void advectScalars(double dt);   // 平流温度和密度场
            void pressureProjection(double dt);  // 压力投影（保证不可压缩）
            void gridToParticle(double dt);  // G2P：网格速度转移回粒子
            void advectParticles(double dt); // 移动粒子
        };
    }
}

#endif

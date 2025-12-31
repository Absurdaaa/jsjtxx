#pragma once
#ifndef __LAGRANGIAN_3D_SOLVER_H__
#define __LAGRANGIAN_3D_SOLVER_H__

#include "ParticleSystem3d.h"
#include "Global.h"
#include "Configure.h"
#include <iostream>

namespace FluidSimulation
{
    namespace Lagrangian3d
    {
        // 拉格朗日法流体求解器类
        // 负责求解基于粒子的三维流体动力学方程
        class Solver
        {
        public:
            Solver(ParticleSystem3d &ps);

            void solve();                         // 求解流体方程

        private:
            float computeDensity(particle3d &pi);
            float computePressure(particle3d &pi);
            glm::vec3 computeAcceleration(particle3d &pi);

            void integrateParticle(particle3d &pi, float dt);
            void boundaryCheck(particle3d &pi);
            void updateBlockId(particle3d &pi);

            float poly6Kernel(float dist2) const;
            glm::vec3 spikyGradient(const glm::vec3 &rVec, float dist) const;
            float viscosityLaplacian(float dist) const;
            glm::vec3 xsphVelocityCorrection(const particle3d &pi) const;

            static constexpr float PI = 3.14159265358979323846f;

            ParticleSystem3d &mPs; // 粒子系统引用
        };
    }
}

#endif

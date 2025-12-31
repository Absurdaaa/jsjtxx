#pragma once
#ifndef __LAGRANGIAN_2D_FOUNTAIN_COMPONENT_H__
#define __LAGRANGIAN_2D_FOUNTAIN_COMPONENT_H__

#include "Renderer.h"
#include "Solver.h"
#include "ParticleSystem2d.h"

#include "Component.h"
#include "Configure.h"
#include "Logger.h"
#include "Global.h"

namespace FluidSimulation
{
    namespace Lagrangian2d
    {
        // 2D喷泉专用组件
        // 提供持续粒子发射与渲染逻辑
        class Lagrangian2dFountainComponent : public Glb::Component
        {
        public:
            Renderer *renderer;
            Solver *solver;
            ParticleSystem2d *ps;

            Lagrangian2dFountainComponent(char *description, int id)
            {
                this->description = description;
                this->id = id;
                renderer = NULL;
                solver = NULL;
                ps = NULL;
            }

            virtual void shutDown() override;
            virtual void init() override;
            virtual void simulate() override;
            virtual GLuint getRenderedTexture() override;

        private:
            Glb::RandomGenerator mRand;
            void emitParticles();
            ParticleInfo2d buildParticle();
            void removeSettledParticles();
        };
    }
}

#endif

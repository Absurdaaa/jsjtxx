#pragma once
#ifndef __PIC_COMPONENT_2D_H__
#define __PIC_COMPONENT_2D_H__

#include "Component.h"
#include "Renderer.h"
#include "Solver.h"
#include "PICGrid2d.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        class PICComponent : public Glb::Component
        {
        public:
            PICComponent(const char *desc, int id)
            {
                this->description = const_cast<char *>(desc);
                this->id = id;
            }

            virtual void shutDown() override;
            virtual void init() override;
            virtual void simulate() override;
            virtual GLuint getRenderedTexture() override;

        private:
            ParticleSystem *ps = nullptr;
            PICGrid2d *grid = nullptr;
            Solver *solver = nullptr;
            Renderer *renderer = nullptr;
        };
    }
}

#endif

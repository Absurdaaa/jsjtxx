#include "PIC/include/PICComponent.h"
#include "Logger.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
        void PICComponent::shutDown()
        {
            delete renderer; renderer = nullptr;
            delete solver; solver = nullptr;
            delete ps; ps = nullptr;
            delete grid; grid = nullptr;
        }

        void PICComponent::init()
        {
            shutDown();
            grid = new PICGrid2d();
            ps = new ParticleSystem();
            solver = new Solver(*ps, *grid);
            renderer = new Renderer();
            Glb::Logger::getInstance().addLog("PIC 2d initialized");
        }

        void PICComponent::simulate()
        {
            if (solver)
                solver->solve();
        }

        GLuint PICComponent::getRenderedTexture()
        {
            if (renderer && ps && grid)
                renderer->draw(*ps, *grid);
            return renderer ? renderer->getTextureID() : 0;
        }
    }
}

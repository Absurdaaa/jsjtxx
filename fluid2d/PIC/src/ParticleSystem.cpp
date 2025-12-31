#include "PIC/include/ParticleSystem.h"
#include <random>

namespace FluidSimulation
{
    namespace PIC2d
    {
        void ParticleSystem::emitFromSources()
        {
            static std::mt19937 rng(42);
            std::uniform_real_distribution<float> uni(-0.5f, 0.5f);

            const float h = Eulerian2dPara::theCellSize2d;
            const int emitN = PIC2dPara::particlesPerStep;
            const float jitter = PIC2dPara::emissionJitter * h;

            for (const auto &src : Eulerian2dPara::source)
            {
                // 源位置：网格坐标转世界坐标（单元中心）
                float baseX = (src.position.x + 0.5f) * h;
                float baseY = (src.position.y + 0.5f) * h;

                for (int n = 0; n < emitN; ++n)
                {
                    Particle p;
                    p.position.x = baseX + uni(rng) * jitter;
                    p.position.y = baseY + uni(rng) * jitter;
                    p.velocity = src.velocity;
                    particles.push_back(p);
                }
            }
        }

        void ParticleSystem::clampToDomain(const glm::vec2 &lower, const glm::vec2 &upper)
        {
            for (auto &p : particles)
            {
                p.position = glm::clamp(p.position, lower, upper);
            }
        }
    }
}

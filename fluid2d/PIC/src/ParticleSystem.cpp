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

            const float h = PIC2dPara::theCellSize2d;
            const int emitN = PIC2dPara::particlesPerStep;
            const float jitter = PIC2dPara::emissionJitter * h;
            const int radius = PIC2dPara::emitterRadius; // emitter radius in cells
            std::uniform_int_distribution<int> idist(radius > 0 ? -radius : 0, radius > 0 ? radius : 0);

            std::uniform_int_distribution<int> xdist(-1, 1);

            for (const auto &src : PIC2dPara::source)
            {
                // 发射位置：可以在源周围 radius 格子内随机分布，形成更粗的喷射
                for (int n = 0; n < emitN; ++n)
                {
                    int dx = xdist(rng);
                    int dy = idist(rng);
                    int cellX = src.position.x + dx;
                    int cellY = src.position.y + dy;
                    // 跳过越界的粒子，而不是 clamp（避免边界堆积）
                    if (cellX < 0 || cellX >= PIC2dPara::theDim2d[0] ||
                        cellY < 0 || cellY >= PIC2dPara::theDim2d[1])
                        continue;

                    float baseX = (cellX + 0.5f) * h;
                    float baseY = (cellY + 0.5f) * h;

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

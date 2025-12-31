#include "PIC/include/ParticleSystem.h"
#include <random>

namespace FluidSimulation
{
    namespace PIC2d
    {
        void ParticleSystem::emitFromSources()
        {
            std::mt19937 rng(std::random_device{}());
            // 均匀分布 [0,1)
            std::uniform_real_distribution<float> uni01(0.0f, 1.0f);

            // 每个源每步发射的粒子数
            const int emitN = PIC2dPara::particlesPerStep;
            for (const auto &src : Eulerian2dPara::source)
            {
                for (int n = 0; n < emitN; ++n)
                {
                    Particle p;
                    // 采样公式：x = (i + 0.5 + ξ) * h，y = (j + 0.5 + η) * h
                    // 其中 ξ,η ∈ [-0.5 * jitter, 0.5 * jitter]，保证粒子均匀分布在单元中心附近
                    glm::vec2 base((src.position.x + 0.5f) * Eulerian2dPara::theCellSize2d,
                                   (src.position.y + 0.5f) * Eulerian2dPara::theCellSize2d);
                    float jx = (uni01(rng) - 0.5f) * PIC2dPara::emissionJitter * Eulerian2dPara::theCellSize2d;
                    float jy = (uni01(rng) - 0.5f) * PIC2dPara::emissionJitter * Eulerian2dPara::theCellSize2d;
                    p.position = base + glm::vec2(jx, jy);
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

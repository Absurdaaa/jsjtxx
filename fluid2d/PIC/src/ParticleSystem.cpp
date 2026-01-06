#include "PIC/include/ParticleSystem.h"
#include <random>
#include <algorithm>
#include <cmath>

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

            const int nx = PIC2dPara::theDim2d[0];
            const int ny = PIC2dPara::theDim2d[1];

            // emitterRadius: 以“格子数”为半径的喷口（用于在源周围发射一团粒子）
            const int radius = 0 > PIC2dPara::emitterRadius ? 0 : PIC2dPara::emitterRadius;

            for (const auto &src : PIC2dPara::source)
            {
                // 关键：不要用 clamp 把越界采样“夹回”边界格子（会在边界堆叠，从而看起来偏在某一半）
                // 正确做法：根据源位置算出合法的 dx/dy 取值范围。
                const int sx = src.position.x;
                const int sy = src.position.y;

                const int minDx = (radius > 0) ? -(sx < radius ? sx : radius) : 0;
                const int maxDx = (radius > 0) ? ((nx - 1) - sx < radius ? (nx - 1) - sx : radius) : 0;
                const int minDy = (radius > 0) ? -(sy < radius ? sy : radius) : 0;
                const int maxDy = (radius > 0) ?  ((ny - 1) - sy < radius ? (ny - 1) - sy : radius) : 0;

                // 你描述的喷口是“左边界上的一段开口”。
                // 所以当源在 x=0 时：
                // - 只在 y 方向扩展（形成一段竖直开口）
                // - x 方向只给很小的厚度（避免生成在墙外/或被边界堆叠）
                const bool boundaryJet = (sx == 0);

                std::uniform_int_distribution<int> dyDist(minDy, maxDy);
                std::uniform_int_distribution<int> dxDist(
                    boundaryJet ? 0 : minDx,
                    boundaryJet ? (radius > 0 ? (maxDx < 2 ? maxDx : 2) : 0) : maxDx);

                // 发射位置：可以在源周围 radius 格子内随机分布，形成更粗的喷射
                for (int n = 0; n < emitN; ++n)
                {
                    int dx = dxDist(rng);
                    int dy = dyDist(rng);

                    // 如果不是边界喷口（例如你将来把源放到域内部），才使用“圆形块”采样
                    if (!boundaryJet && radius > 0)
                    {
                        const int r2 = radius * radius;
                        int tries = 0;
                        while (tries < 8 && (dx * dx + dy * dy > r2))
                        {
                            dx = dxDist(rng);
                            dy = dyDist(rng);
                            ++tries;
                        }
                    }

                    int cellX = sx + dx;
                    int cellY = sy + dy;

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

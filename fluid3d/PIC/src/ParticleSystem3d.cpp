/**
 * ParticleSystem3d.cpp: 3D PIC 粒子系统实现
 * 实现粒子发射和边界约束功能
 */

#include "PIC/include/ParticleSystem3d.h"
#include <random>

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 从烟雾源发射新粒子
         * 在每个源位置周围随机抖动生成粒子
         * 相关成员：
         * - Particle::position / Particle::velocity: 新粒子的属性
         * - ParticleSystem3d::particles: 将新粒子 push_back 到该容器
         */
        void ParticleSystem3d::emitFromSources()
        {
            // 使用固定种子的随机数生成器，保证可重复性
            static std::mt19937 rng(42);
            std::uniform_real_distribution<float> uni{-0.5f, 0.5f};

            const float h = PIC3dPara::theCellSize3d;        // 网格单元大小
            const int emitN = PIC3dPara::particlesPerStep;        // 每步发射粒子数
            const float jitter = PIC3dPara::emissionJitter * h;   // 发射抖动范围
            const int radius = PIC3dPara::emitterRadius;         // 喷口半径（格子数）：圆盘半径
            std::uniform_int_distribution<int> idist{radius > 0 ? -radius : 0, radius > 0 ? radius : 0};

            // 遍历所有烟雾源
            for (const auto &src : PIC3dPara::source)
            {
                // 发射指定数量的粒子：面上的圆形喷口（在 y-z 平面内为圆盘）
                for (int n = 0; n < emitN; ++n)
                {
                    // 固定在一个面附近：x 方向只给一个很薄的厚度（默认 0）
                    const int cellX = src.position.x;

                    int cellY = src.position.y;
                    int cellZ = src.position.z;

                    // 在 y-z 平面做圆盘采样；越界就重采样，避免 clamp 带来的偏置
                    if (radius > 0)
                    {
                        bool ok = false;
                        for (int tries = 0; tries < 16 && !ok; ++tries)
                        {
                            int dy = idist(rng);
                            int dz = idist(rng);
                            if (dy * dy + dz * dz > radius * radius)
                                continue;

                            int y = src.position.y + dy;
                            int z = src.position.z + dz;
                            if (y < 0 || y >= PIC3dPara::theDim3d[1] || z < 0 || z >= PIC3dPara::theDim3d[2])
                                continue;

                            cellY = y;
                            cellZ = z;
                            ok = true;
                        }
                    }

                    // 兜底：确保合法
                    cellY = (cellY < 0) ? 0 : (cellY >= PIC3dPara::theDim3d[1] ? PIC3dPara::theDim3d[1] - 1 : cellY);
                    cellZ = (cellZ < 0) ? 0 : (cellZ >= PIC3dPara::theDim3d[2] ? PIC3dPara::theDim3d[2] - 1 : cellZ);

                    float baseX = (cellX + 0.5f) * h;
                    float baseY = (cellY + 0.5f) * h;
                    float baseZ = (cellZ + 0.5f) * h;

                    Particle p;
                    // 在单元内 jitter
                    p.position.x = baseX + uni(rng) * jitter;
                    p.position.y = baseY + uni(rng) * jitter;
                    p.position.z = baseZ + uni(rng) * jitter;
                    p.velocity = src.velocity;  // 继承源的初始速度
                    particles.push_back(p);
                }
            }
        }

        /**
         * 将所有粒子位置限制在指定域内
         */
        void ParticleSystem3d::clampToDomain(const glm::vec3 &lower, const glm::vec3 &upper)
        {
            for (auto &p : particles)
            {
                p.position = glm::clamp(p.position, lower, upper);
            }
        }
    }
}

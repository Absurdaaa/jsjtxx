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
            std::uniform_real_distribution<float> uni(-0.5f, 0.5f);

            const float h = PIC3dPara::theCellSize3d;        // 网格单元大小
            const int emitN = PIC3dPara::particlesPerStep;        // 每步发射粒子数
            const float jitter = PIC3dPara::emissionJitter * h;   // 发射抖动范围
            const int radius = PIC3dPara::emitterRadius;         // 发射源半径（格子数）
            std::uniform_int_distribution<int> idist(radius > 0 ? -radius : 0, radius > 0 ? radius : 0);

            // 遍历所有烟雾源
            for (const auto &src : PIC3dPara::source)
            {
                // 发射指定数量的粒子，发射源可为球体（以格子为单位的半径）
                for (int n = 0; n < emitN; ++n)
                {
                    int dx = idist(rng);
                    int dy = idist(rng);
                    int dz = idist(rng);

                    // 若使用球形发射区，需要满足 dx^2+dy^2+dz^2 <= radius^2
                    if (radius > 0 && (dx * dx + dy * dy + dz * dz > radius * radius))
                    {
                        // 若不满足，重试一次（简单策略），否则退回到中心
                        dx = 0; dy = 0; dz = 0;
                    }

                    int cellX = src.position.x + dx;
                    int cellY = src.position.y + dy;
                    int cellZ = src.position.z + dz;

                    // clamp 到网格范围
                    cellX = cellX < PIC3dPara::theDim3d[0] ? cellX : PIC3dPara::theDim3d[0] - 1;
                    cellX = cellX > 0 ? cellX : 0;
                    cellY = cellY < PIC3dPara::theDim3d[1] ? cellY : PIC3dPara::theDim3d[1] - 1;
                    cellY = cellY > 0 ? cellY : 0;
                    cellZ = cellZ < PIC3dPara::theDim3d[2] ? cellZ : PIC3dPara::theDim3d[2] - 1;
                    cellZ = cellZ > 0 ? cellZ : 0;

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

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
         */
        void ParticleSystem3d::emitFromSources()
        {
            // 使用固定种子的随机数生成器，保证可重复性
            static std::mt19937 rng(42);
            std::uniform_real_distribution<float> uni(-0.5f, 0.5f);

            const float h = Eulerian3dPara::theCellSize3d;        // 网格单元大小
            const int emitN = PIC3dPara::particlesPerStep;        // 每步发射粒子数
            const float jitter = PIC3dPara::emissionJitter * h;   // 发射抖动范围

            // 遍历所有烟雾源
            for (const auto &src : Eulerian3dPara::source)
            {
                // 计算源的世界坐标（单元中心）
                float baseX = (src.position.x + 0.5f) * h;
                float baseY = (src.position.y + 0.5f) * h;
                float baseZ = (src.position.z + 0.5f) * h;

                // 发射指定数量的粒子
                for (int n = 0; n < emitN; ++n)
                {
                    Particle p;
                    // 在源位置周围随机抖动
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

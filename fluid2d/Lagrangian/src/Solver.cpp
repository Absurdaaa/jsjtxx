/**
 * Solver.cpp: 2D拉格朗日流体求解器实现文件
 * 实现基于粒子的2D流体仿真算法
 */

#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        /**
         * 构造函数，保存对粒子系统的引用
         * @param ps 粒子系统引用
         */
        Solver::Solver(ParticleSystem2d &ps) : mPs(ps)
        {
        }

        /**
         * 求解流体方程
         * 实现一步粒子流体的仿真计算
         */
        void Solver::solve()
        {
            // TODO: 实现求解器
            // 求解流体模拟的主要步骤:
            // 1. 计算密度 - 使用SPH核函数计算每个粒子的密度
            // 2. 计算压力 - 根据密度计算压力
            // 3. 计算加速度 - 计算压力梯度和粘性力
            // 4. 更新速度和位置 - 使用欧拉积分更新粒子状态
            // 5. 边界检查 - 处理粒子与边界的碰撞
            // 6. 更新块ID - 更新粒子的空间索引
            // ...

            // 使用引用直接操作粒子系统中的粒子数据
            // 子步数
            int substeps = substeps_ > 0 ? substeps_ :1;
            for (int s = 0; s < substeps; s++)
            {
                // 更新块信息以便进行邻域查询
                mPs.updateBlockInfo();

                // ========== 第二步：创建临时缓存（存储本轮计算结果） ==========
                std::vector<glm::vec2> temp_acceleration(mPs.particles.size());

                // ========== 第三步：计算所有粒子的密度 ==========
                for (size_t i = 0; i < mPs.particles.size(); i++)
                {
                  auto &p = mPs.particles[i];
                  p.density = cul_density(p);
                }

                // ========== 第四步：计算所有粒子的压力（存入缓存） ==========
                for (size_t i = 0; i < mPs.particles.size(); i++)
                {
                  auto &p = mPs.particles[i];
                  p.pressure = cul_pressure(p);
                }

                // ========== 第五步：计算所有粒子的加速度（存入缓存） ==========
                for (size_t i = 0; i < mPs.particles.size(); i++)
                {
                  auto &p = mPs.particles[i];
                  temp_acceleration[i] = cul_acceleration(p); // 基于原位置/密度/压力计算
                }

                float dt = dt_ / float(substeps);
                for (size_t i = 0; i < mPs.particles.size(); i++)
                {
                  auto &p = mPs.particles[i];
                  p.accleration = temp_acceleration[i];
                }

                // 4) 更新速度、位置并做边界检测、更新块ID
                for (auto &p : mPs.particles)
                {
                    update_velocity_position(p, dt);
                    boundary_check(p);
                    update_blockId(p);
                }
            }

          // 最后更新一次块信息，保证索引正确
          mPs.updateBlockInfo();
        }

        float Solver::poly6Kernel(float dist2, float supportRadius, float supportRadius2)
        {
          // 传入距离平方，避免重复计算
          if (dist2 >= supportRadius2)
            return 0.0f;

          float h = supportRadius;
          float h8 = pow(h, 8); // 这里可以踢出去减少计算
          float h2_minus_r2 = supportRadius2 - dist2;

          // Poly6公式：4/(πh^8) * (h? - r?)?
          float numerator = 4.0f * pow(h2_minus_r2, 3);
          float denominator = PI * h8;

          return numerator / denominator;
        }

        // 计算粒子密度（使用邻域搜索加速）
        float Solver::cul_density(ParticleInfo2d &pi)
        {
            float density = 0.0f;

            int blocks = (int)mPs.blockExtens.size();
            int curBlock = -1;
            // 检查 blockId 是否为无效标记 UINT32_MAX（getBlockIdByPosition 在越界时可能返回 -1 -> uint32_t 最大值）
            if (pi.blockId != UINT32_MAX)
                curBlock = (int)pi.blockId;

            if (curBlock < 0 || curBlock >= blocks)
            {
                // 粒子不在任何块中，返回恢复密度（使用参考密度，避免后续计算异常）
                return density_;
            }

            const float h = mPs.supportRadius;
            const float h2 = mPs.supportRadius2;

            // 遍历相邻 9 个块
            for (int k = 0; k < (int)mPs.blockIdOffs.size(); k++)
            {
                int nb = curBlock + mPs.blockIdOffs[k];
                if (nb < 0 || nb >= blocks)
                    continue;

                glm::uvec2 range = mPs.blockExtens[nb];
                int left = (int)range.x;
                int right = (int)range.y;

                for (int idx = left; idx < right; idx++)
                {
                    const auto &pj = mPs.particles[idx];
                    glm::vec2 dist = pj.position - pi.position;
                    float dist2 = glm::dot(dist, dist);
                    if (dist2 >= h2)
                        continue;
                    // 使用每粒子的质量来累加密度：m_j * W(r)
                    // 这里的粒子质量取为 ρ0 * particleVolume（在其它地方也使用该定义）
                    float mass_j = density_ * mPs.particleVolume;
                    density += mass_j * poly6Kernel(dist2, h, h2);
                }
            }

            // 最小限保护，避免密度为 0
            if (density <= 1e-8f)
                density = density_ * 0.5f;

            return density;
        }


        // Spiky核函数（2D，压力梯度计算专用）
        float Solver::SpikyKernel(float dist, float supportRadius)
        {
            if (dist >= supportRadius)
                return 0.0f;

            float h = supportRadius;
            float h2 = pow(h, 2);
            float h_minus_r = h - dist;

            // Spiky公式：10/(πh^2) * (h - r)^2
            float numerator = 10.0f * h_minus_r * h_minus_r;
            float denominator = PI * h2;

            return numerator / denominator;
          
        }
        // 根据密度计算压力（使用 Tait/幂律形式：p = stiffness * ((rho/rho0)^exponent - 1)）
        float Solver::cul_pressure(ParticleInfo2d &pi)
        {
            const float rho0 = density_;
            const float k = stiffness_;
            const float gamma = exponent_;

            // 防止 rho0 为 0
            if (rho0 <= 0.0f)
            {
                pi.pressure = 0.0f;
            }
            else
            {
                float ratio = pi.density / rho0;
                // gongshi  p = k * (ratio^gamma - 1)
                pi.pressure = k * (pow(ratio, gamma) - 1.0f);
                
            }

            // 存储 pressure / density^2 用于力计算的优化
            float denom = pi.density * pi.density + 1e-8f;
            pi.pressDivDens2 = pi.pressure / denom;
            return pi.pressure;
        }
    
        // 计算粒子受力并返回加速度
        glm::vec2 Solver::cul_acceleration(ParticleInfo2d &pi)
        {
            // 初始化加速度为零
            glm::vec2 acc(0.0f);

            int blocks = (int)mPs.blockExtens.size();
            int curBlock = -1;
            if (pi.blockId != UINT32_MAX)
                curBlock = (int)pi.blockId;
            if (curBlock < 0 || curBlock >= blocks)
                curBlock = -1; // 会导致跳过邻域

            const float h = mPs.supportRadius;
            const float h2 = mPs.supportRadius2;
            const float mass = density_ * mPs.particleVolume; // 统一使用与密度计算一致的粒子质量

            // 预计算 poly6 常数相关项用于 gradient
            const float h8 = pow(h, 8);
            const float h6 = pow(h, 6);
            
            const float nu = viscosity_; // 粘度系数

            // std::cout << "blockIdOffs size: " << mPs.blockIdOffs.size() << std::endl; // 应输出9

            for (int k = 0; k < (int)mPs.blockIdOffs.size(); k++)
            {
                if (curBlock < 0)
                    {
                      // std::cout << "Warning: Particle out of bounds in cul_acceleration()" << std::endl;
                      break;
                    }
                int nb = curBlock + mPs.blockIdOffs[k];
                if (nb < 0 || nb >= blocks)
                    continue;

                glm::uvec2 range = mPs.blockExtens[nb];
                int left = (int)range.x;
                int right = (int)range.y;
                
                float rho_i = pi.density;

                // std::cout << "Block " << nb << " range: " << left << " ~ " << right << std::endl;

                for (int idx = left; idx < right; idx++)
                {
                    const auto &pj = mPs.particles[idx];
                    if (&pj == &pi)continue;

                    glm::vec2 r_vec = pi.position - pj.position;
                    float dist2 = glm::dot(r_vec, r_vec);
                    if (dist2 >= h2)
                        {
                          // std::cout<< "Skip particle due to dist2 >= h2: " << dist2 << " >= " << h2 << std::endl;
                          continue;
                        }

                    // 距离
                    float dist = sqrt(dist2 + 1e-12f);

                    // --- 压力梯度 ---
                    // 2D Spiky 核梯度 计算
                    // ?W(r,h) = -30/(πh^5) * (h - r)^2 * (r_vec / r)
                    float spiky_coeff = -30.0f / (PI * pow(h, 5));
                    float factor = spiky_coeff * (h - dist) * (h - dist) / (dist + 1e-12f); // 避免除0
                    glm::vec2 gradW = factor * r_vec;

                    float pij = pi.pressDivDens2 + pj.pressDivDens2; // p_i/rho_i^2 + p_j/rho_j^2
                    glm::vec2 press_acc = -mass * pij * gradW;
                    
                    
                    // std::cout << "press_acc: " << press_acc.x << "," << press_acc.y << std::endl;
                    acc += press_acc;

                    // --- 粘性项 ---
                    // 1. 计算粒子i-j的矢量和距离
                    float r = dist;
                    if (r > h || r < 1e-8)
                    { // 超出核支撑域/粒子重合，无贡献
                      // std::cout << "Skip particle due to r > h or r < 1e-8: " << r << std::endl;
                      continue;
                    }
                    // 2. 2D 粘性核的 Laplacian
                    // ??W(r,h) = 40/(πh^5) * (h - r)
                    float laplacian_coeff = 40.0f / (PI * pow(h, 5));
                    float laplacianW = laplacian_coeff * (h - r);

                    // 3. 速度差（二维矢量）
                    glm::vec2 vel_diff = pj.velocity - pi.velocity;

                    // 4. 对称形式的粘性加速度贡献（核心计算）
                    // float rho_i = pi.density;
                    float rho_j = pj.density;
                    // 粘性加速度贡献（矢量）
                    glm::vec2 visc_acc = nu * mass * (pj.velocity - pi.velocity) / (pj.density + 1e-8f) * laplacianW;

                    acc += visc_acc;
                    
                }
            }
            // 加上重力
            // std::cout << "acc - g: " << acc.x << "," << acc.y << std::endl;

            acc -= glm::vec2(gravityX_, gravityY_);
            
            return acc;
        }

        // 更新速度和位置（欧拉积分）
        void Solver::update_velocity_position(ParticleInfo2d &pi, float dt)
        {
            // 显式欧拉,更新速度
            pi.velocity += pi.accleration * dt;

            // 限制最大速度
            float maxV = maxVelocity_;
            float vlen2 = glm::dot(pi.velocity, pi.velocity);
            if (maxV > 0.0f && vlen2 > maxV * maxV)
            {
                pi.velocity = glm::normalize(pi.velocity) * maxV;
            }

            pi.position += pi.velocity * dt;
        }

        // 简单的边界碰撞处理（反弹并衰减速度）
        void Solver::boundary_check(ParticleInfo2d &pi)
        {
            const glm::vec2 &lb = mPs.lowerBound;
            const glm::vec2 &ub = mPs.upperBound;
            float r = mPs.particleRadius;

            // X 方向
            if (pi.position.x < lb.x + r)
            {
                pi.position.x = lb.x + r;
                if (pi.velocity.x < 0.0f)
                    pi.velocity.x = -pi.velocity.x * velocityAttenuation_;
            }
            else if (pi.position.x > ub.x - r)
            {
                pi.position.x = ub.x - r;
                if (pi.velocity.x > 0.0f)
                    pi.velocity.x = -pi.velocity.x * velocityAttenuation_;
            }

            // Y 方向
            if (pi.position.y < lb.y + r)
            {
                pi.position.y = lb.y + r;
                if (pi.velocity.y < 0.0f)
                    pi.velocity.y = -pi.velocity.y * velocityAttenuation_;
            }
            else if (pi.position.y > ub.y - r)
            {
                pi.position.y = ub.y - r;
                if (pi.velocity.y > 0.0f)
                    pi.velocity.y = -pi.velocity.y * velocityAttenuation_;
            }
        }

        // 根据位置更新粒子的块ID
        void Solver::update_blockId(ParticleInfo2d &pi)
        {
            uint32_t id = mPs.getBlockIdByPosition(pi.position);
            pi.blockId = id;
        }
    }
}

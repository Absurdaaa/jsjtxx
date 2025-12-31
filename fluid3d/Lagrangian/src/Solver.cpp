/**
 * Solver.cpp: 3D拉格朗日流体求解器实现文件
 * 实现基于粒子的3D流体仿真算法
 */

#include "fluid3d/Lagrangian/include/Solver.h"

#include <algorithm>
#include <cmath>

namespace FluidSimulation
{

	namespace Lagrangian3d
	{
		/**
		 * 构造函数，保存粒子系统引用
		 * @param ps 粒子系统引用
		 */
		Solver::Solver(ParticleSystem3d &ps) : mPs(ps)
		{
			// 构造函数,保存粒子系统引用
		}

		/**
		 * 求解流体方程
		 * 实现一步3D粒子流体的仿真计算
		 */
		void Solver::solve()
		{
			if (mPs.particles.empty())
			{
				return;
			}

			// === 统一时间步设置 ===
			// substeps 用于细化时间步，提升稳定性；dt 是每个子步的积分长度
			const int substeps = (Lagrangian3dPara::substep > 0) ? Lagrangian3dPara::substep : 1;
			const float dt = Lagrangian3dPara::dt / static_cast<float>(substeps);
			const bool useXsph = Lagrangian3dPara::xsph_c > 0.0f;

			for (int step = 0; step < substeps; ++step)
			{
				// 每个子步开始时刷新块信息，保证邻域查询使用最新的排序结果
				mPs.updateBlockInfo();

				const size_t count = mPs.particles.size();
				// 加速度缓存，避免在密度/压力迭代期间覆盖原值
				std::vector<glm::vec3> tempAcceleration(count, glm::vec3(0.0f));
				std::vector<glm::vec3> xsphDelta;
				if (useXsph)
				{
					// 只有开启 XSPH 才分配修正缓存，节约内存
					xsphDelta.resize(count, glm::vec3(0.0f));
				}

				// 1. 基于核函数累加周围粒子质量，得到当前粒子密度
				for (size_t i = 0; i < count; ++i)
				{
					auto &p = mPs.particles[i];
					p.density = computeDensity(p);
				}

				// 2. 使用 Tait 方程从密度推导压力
				for (size_t i = 0; i < count; ++i)
				{
					auto &p = mPs.particles[i];
					p.pressure = computePressure(p);
				}

				// 3. 计算压力梯度、粘性项以及外力，得到加速度
				for (size_t i = 0; i < count; ++i)
				{
					auto &p = mPs.particles[i];
					tempAcceleration[i] = computeAcceleration(p);
				}

				if (useXsph)
				{
					// 4. 可选的 XSPH 粘性，通过速度平滑减少粒子噪声
					for (size_t i = 0; i < count; ++i)
					{
						xsphDelta[i] = xsphVelocityCorrection(mPs.particles[i]);
					}
				}

				// 5. 写回加速度与速度，积分更新位置，并处理边界/索引
				for (size_t i = 0; i < count; ++i)
				{
					auto &p = mPs.particles[i];
					p.accleration = tempAcceleration[i];
					if (useXsph)
					{
						p.velocity += xsphDelta[i];
					}
					integrateParticle(p, dt);
					boundaryCheck(p);
					updateBlockId(p);
				}
			}

			mPs.updateBlockInfo();
		}

		float Solver::computeDensity(particle3d &pi)
		{
			float density = 0.0f;
			const float h2 = mPs.supportRadius2;
			const int blocks = static_cast<int>(mPs.blockExtens.size());
			int curBlock = (pi.blockId != UINT32_MAX) ? static_cast<int>(pi.blockId) : -1;
			if (curBlock < 0 || curBlock >= blocks)
			{
				return Lagrangian3dPara::density;
			}

			// 遍历当前块及 26 个邻居块，累加核权重，实现近邻搜索剪枝
			for (size_t k = 0; k < mPs.blockIdOffs.size(); ++k)
			{
				int nb = curBlock + mPs.blockIdOffs[k];
				if (nb < 0 || nb >= blocks)
				{
					continue;
				}

				glm::uvec2 range = mPs.blockExtens[nb];
				int left = static_cast<int>(range.x);
				int right = static_cast<int>(range.y);

				for (int idx = left; idx < right; ++idx)
				{
					const auto &pj = mPs.particles[idx];
					glm::vec3 delta = pj.position - pi.position;
					float dist2 = glm::dot(delta, delta);
					if (dist2 >= h2)
					{
						continue;
					}
					density += mPs.particleMass * poly6Kernel(dist2);
				}
			}

			const float minDensity = 0.5f * Lagrangian3dPara::density;
      return minDensity < density ? density : minDensity;
    }

		float Solver::computePressure(particle3d &pi)
		{
			const float rho0 = Lagrangian3dPara::density > 1e-8f ? Lagrangian3dPara::density : 1000.0f;
      const float stiffness = Lagrangian3dPara::stiffness;
			const float exponent = Lagrangian3dPara::exponent;

			float ratio = pi.density / rho0;
			pi.pressure = stiffness * (std::pow(ratio, exponent) - 1.0f);
			float denom = pi.density * pi.density + 1e-8f;
			pi.pressDivDens2 = pi.pressure / denom;
			return pi.pressure;
		}

		glm::vec3 Solver::computeAcceleration(particle3d &pi)
		{
			glm::vec3 acc(0.0f);
			const float h2 = mPs.supportRadius2;
			const int blocks = static_cast<int>(mPs.blockExtens.size());
			int curBlock = (pi.blockId != UINT32_MAX) ? static_cast<int>(pi.blockId) : -1;
			if (curBlock < 0 || curBlock >= blocks)
			{
				curBlock = -1;
			}

			for (size_t k = 0; k < mPs.blockIdOffs.size(); ++k)
			{
				if (curBlock < 0)
				{
					break;
				}
				int nb = curBlock + mPs.blockIdOffs[k];
				if (nb < 0 || nb >= blocks)
				{
					continue;
				}

				glm::uvec2 range = mPs.blockExtens[nb];
				int left = static_cast<int>(range.x);
				int right = static_cast<int>(range.y);

				for (int idx = left; idx < right; ++idx)
				{
					const auto &pj = mPs.particles[idx];
					if (&pj == &pi)
					{
						continue;
					}

					glm::vec3 rVec = pi.position - pj.position;
					float dist2 = glm::dot(rVec, rVec);
					if (dist2 >= h2)
					{
						continue;
					}

					float dist = std::sqrt(dist2 + 1e-12f);
					glm::vec3 gradW = spikyGradient(rVec, dist);
					// 压力项：对称形式 (p_i/ρ_i^2 + p_j/ρ_j^2) * ?W，确保动量守恒
					float pressureTerm = pi.pressDivDens2 + pj.pressDivDens2;
					acc += -mPs.particleMass * pressureTerm * gradW;

					float lapW = viscosityLaplacian(dist);
					// 粘性项：使用拉普拉斯核平衡速度差，使流动更光滑
					glm::vec3 visc = Lagrangian3dPara::viscosity * mPs.particleMass * (pj.velocity - pi.velocity) / (pj.density + 1e-8f) * lapW;
					acc += visc;
				}
			}

			// 外力：简单重力模型，方向遵循配置文件
			acc -= glm::vec3(Lagrangian3dPara::gravityX, Lagrangian3dPara::gravityY, Lagrangian3dPara::gravityZ);
			return acc;
		}

		void Solver::integrateParticle(particle3d &pi, float dt)
		{
			pi.velocity += pi.accleration * dt;
			float maxV = Lagrangian3dPara::maxVelocity;
			if (maxV > 0.0f)
			{
				float speed2 = glm::dot(pi.velocity, pi.velocity);
				float maxV2 = maxV * maxV;
				if (speed2 > maxV2)
				{
					pi.velocity = glm::normalize(pi.velocity) * maxV;
				}
			}

			pi.position += pi.velocity * dt;
		}

		void Solver::boundaryCheck(particle3d &pi)
		{
			const glm::vec3 &lb = mPs.lowerBound;
			const glm::vec3 &ub = mPs.upperBound;
			const float r = mPs.particleRadius;
			const float eps = Lagrangian3dPara::eps;
			const float attenuation = Lagrangian3dPara::velocityAttenuation;

			auto reflectAxis = [&](float &pos, float &vel, float lower, float upper)
			{
				// 超出边界时，将粒子拉回并反向速度，同时乘以衰减系数实现“弹性但有损”的碰撞
				if (pos < lower + r)
				{
					pos = lower + r + eps;
					if (vel < 0.0f)
					{
						vel = -vel * attenuation;
					}
				}
				else if (pos > upper - r)
				{
					pos = upper - r - eps;
					if (vel > 0.0f)
					{
						vel = -vel * attenuation;
					}
				}
			};

			reflectAxis(pi.position.x, pi.velocity.x, lb.x, ub.x);
			reflectAxis(pi.position.y, pi.velocity.y, lb.y, ub.y);
			reflectAxis(pi.position.z, pi.velocity.z, lb.z, ub.z);
		}

		void Solver::updateBlockId(particle3d &pi)
		{
			pi.blockId = mPs.getBlockIdByPosition(pi.position);
		}

		float Solver::poly6Kernel(float dist2) const
		{
			if (dist2 >= mPs.supportRadius2)
			{
				return 0.0f;
			}

			const float h = mPs.supportRadius;
			const float h2 = mPs.supportRadius2;
			// Poly6 核：W(r,h) = 315/(64πh^9) * (h^2 - r^2)^3，适合密度计算
			const float coeff = 315.0f / (64.0f * PI * std::pow(h, 9));
			float t = h2 - dist2;
			return coeff * t * t * t;
		}

		glm::vec3 Solver::spikyGradient(const glm::vec3 &rVec, float dist) const
		{
			const float h = mPs.supportRadius;
			if (dist < 1e-6f || dist >= h)
			{
				return glm::vec3(0.0f);
			}

			// Spiky 核梯度：?W(r,h) = -45/(πh^6) * (h - r)^2 * r?，主要用于压力
			const float coeff = -45.0f / (PI * std::pow(h, 6));
			float factor = coeff * (h - dist) * (h - dist) / dist;
			return factor * rVec;
		}

		float Solver::viscosityLaplacian(float dist) const
		{
			const float h = mPs.supportRadius;
			if (dist >= h)
			{
				return 0.0f;
			}

			// Viscosity 核拉普拉斯：??W(r,h) = 45/(πh^6) * (h - r)，与速度差配合形成粘性力
			const float coeff = 45.0f / (PI * std::pow(h, 6));
			return coeff * (h - dist);
		}

		glm::vec3 Solver::xsphVelocityCorrection(const particle3d &pi) const
		{
			glm::vec3 correction(0.0f);
			const float h2 = mPs.supportRadius2;
			const int blocks = static_cast<int>(mPs.blockExtens.size());
			int curBlock = (pi.blockId != UINT32_MAX) ? static_cast<int>(pi.blockId) : -1;
			if (curBlock < 0 || curBlock >= blocks)
			{
				return correction;
			}

			// XSPH：对速度做局部平均，能明显抑制粒子“杂点”现象
			for (size_t k = 0; k < mPs.blockIdOffs.size(); ++k)
			{
				int nb = curBlock + mPs.blockIdOffs[k];
				if (nb < 0 || nb >= blocks)
				{
					continue;
				}

				glm::uvec2 range = mPs.blockExtens[nb];
				int left = static_cast<int>(range.x);
				int right = static_cast<int>(range.y);

				for (int idx = left; idx < right; ++idx)
				{
					const auto &pj = mPs.particles[idx];
					if (&pj == &pi)
					{
						continue;
					}

					glm::vec3 delta = pi.position - pj.position;
					float dist2 = glm::dot(delta, delta);
					if (dist2 >= h2)
					{
						continue;
					}

					float kernel = poly6Kernel(dist2);
					float rhoAvg = 0.5f * (pi.density + pj.density) + 1e-8f;
					correction += (mPs.particleMass / rhoAvg) * (pj.velocity - pi.velocity) * kernel;
				}
			}

			return correction * Lagrangian3dPara::xsph_c;
		}
	}
}
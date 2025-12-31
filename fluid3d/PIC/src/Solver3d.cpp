/**
 * Solver3d.cpp: 3D PIC 求解器实现
 * 实现 PIC/FLIP 混合方法的完整求解流程
 */

#include "PIC/include/Solver3d.h"
#include <cmath>

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 构造函数：初始化求解器
         */
        Solver3d::Solver3d(ParticleSystem3d &ps, PICGrid3d &grid)
            : mPs(ps), mGrid(grid)
        {
            // 初始化上一步速度场为零
            mU_prev.initialize(0.0);
            mV_prev.initialize(0.0);
            mW_prev.initialize(0.0);
        }

        /**
         * 执行一次完整的仿真步进
         */
        void Solver3d::solve()
        {
            const double dt = Eulerian3dPara::dt;

            // 1. 发射新粒子
            mPs.emitFromSources();

            // 2. 更新网格源（温度/密度）
            mGrid.updateSources();

            // 3. 平流温度和密度场（半拉格朗日）
            advectScalars(dt);

            // 4. P2G：粒子速度转移到网格
            particleToGrid();

            // 5. 保存旧速度（用于 FLIP）
            mU_prev = mGrid.mU;
            mV_prev = mGrid.mV;
            mW_prev = mGrid.mW;

            // 6. 添加外力（浮力）
            addForces(dt);

            // 7. 压力投影（保证不可压缩）
            pressureProjection(dt);

            // 8. G2P：网格速度转移回粒子
            gridToParticle(dt);

            // 9. 移动粒子
            advectParticles(dt);
        }

        /**
         * P2G：粒子速度转移到网格
         * 使用三线性插值将粒子速度散布到周围的网格节点
         */
        void Solver3d::particleToGrid()
        {
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float h = mGrid.cellSize;

            // 清零速度场
            mGrid.mU.initialize(0.0);
            mGrid.mV.initialize(0.0);
            mGrid.mW.initialize(0.0);

            // 权重累加器
            Glb::GridData3dX wU; wU.initialize(0.0);
            Glb::GridData3dY wV; wV.initialize(0.0);
            Glb::GridData3dZ wW; wW.initialize(0.0);

            for (const auto &p : mPs.particles)
            {
                // === U 分量（位于 YZ 面中心）===
                {
                    float fx = p.position.x / h;
                    float fy = p.position.y / h - 0.5f;
                    float fz = p.position.z / h - 0.5f;
                    int i0 = (int)std::floor(fx), j0 = (int)std::floor(fy), k0 = (int)std::floor(fz);
                    float tx = fx - i0, ty = fy - j0, tz = fz - k0;

                    // 遍历周围 8 个网格节点
                    for (int di = 0; di <= 1; ++di)
                        for (int dj = 0; dj <= 1; ++dj)
                            for (int dk = 0; dk <= 1; ++dk)
                            {
                                int i = i0 + di, j = j0 + dj, k = k0 + dk;
                                if (i < 0 || i > nx || j < 0 || j >= ny || k < 0 || k >= nz) continue;
                                // 三线性插值权重
                                float w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) * (dk ? tz : 1 - tz);
                                mGrid.mU(i, j, k) += w * p.velocity.x;
                                wU(i, j, k) += w;
                            }
                }

                // === V 分量（位于 XZ 面中心）===
                {
                    float fx = p.position.x / h - 0.5f;
                    float fy = p.position.y / h;
                    float fz = p.position.z / h - 0.5f;
                    int i0 = (int)std::floor(fx), j0 = (int)std::floor(fy), k0 = (int)std::floor(fz);
                    float tx = fx - i0, ty = fy - j0, tz = fz - k0;

                    for (int di = 0; di <= 1; ++di)
                        for (int dj = 0; dj <= 1; ++dj)
                            for (int dk = 0; dk <= 1; ++dk)
                            {
                                int i = i0 + di, j = j0 + dj, k = k0 + dk;
                                if (i < 0 || i >= nx || j < 0 || j > ny || k < 0 || k >= nz) continue;
                                float w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) * (dk ? tz : 1 - tz);
                                mGrid.mV(i, j, k) += w * p.velocity.y;
                                wV(i, j, k) += w;
                            }
                }

                // === W 分量（位于 XY 面中心）===
                {
                    float fx = p.position.x / h - 0.5f;
                    float fy = p.position.y / h - 0.5f;
                    float fz = p.position.z / h;
                    int i0 = (int)std::floor(fx), j0 = (int)std::floor(fy), k0 = (int)std::floor(fz);
                    float tx = fx - i0, ty = fy - j0, tz = fz - k0;

                    for (int di = 0; di <= 1; ++di)
                        for (int dj = 0; dj <= 1; ++dj)
                            for (int dk = 0; dk <= 1; ++dk)
                            {
                                int i = i0 + di, j = j0 + dj, k = k0 + dk;
                                if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k > nz) continue;
                                float w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) * (dk ? tz : 1 - tz);
                                mGrid.mW(i, j, k) += w * p.velocity.z;
                                wW(i, j, k) += w;
                            }
                }
            }

            // 归一化：除以权重总和
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i <= nx; ++i)
                        if (wU(i, j, k) > 1e-6f) mGrid.mU(i, j, k) /= wU(i, j, k);

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j <= ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (wV(i, j, k) > 1e-6f) mGrid.mV(i, j, k) /= wV(i, j, k);

            for (int k = 0; k <= nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (wW(i, j, k) > 1e-6f) mGrid.mW(i, j, k) /= wW(i, j, k);
        }

        /**
         * 平流温度和密度场
         * 使用半拉格朗日方法回溯并插值
         */
        void Solver3d::advectScalars(double dt)
        {
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float h = mGrid.cellSize;

            Glb::CubicGridData3d newT, newD;
            newT.initialize(Eulerian3dPara::ambientTemp);
            newD.initialize(0.0);

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k)) continue;
                        // 单元中心位置
                        glm::vec3 pos((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
                        // 半拉格朗日回溯
                        glm::vec3 backPos = mGrid.semiLagrangian(pos, dt);
                        // 插值获取旧值
                        newT(i, j, k) = mGrid.getTemperature(backPos);
                        newD(i, j, k) = mGrid.getDensity(backPos);
                    }

            // 复制回网格
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        mGrid.mT(i, j, k) = newT(i, j, k);
                        mGrid.mD(i, j, k) = newD(i, j, k);
                    }
        }

        /**
         * 添加外力（浮力）
         * 使用 Boussinesq 近似计算浮力
         */
        void Solver3d::addForces(double dt)
        {
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float h = mGrid.cellSize;

            // 在 Z 方向添加浮力
            for (int k = 1; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        glm::vec3 pos((i + 0.5f) * h, (j + 0.5f) * h, k * h);
                        double buoy = mGrid.getBoussinesqForce(pos);
                        mGrid.mW(i, j, k) += dt * buoy;
                    }

            // 边界速度设为 0
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                {
                    mGrid.mU(0, j, k) = 0.0;
                    mGrid.mU(nx, j, k) = 0.0;
                }
            for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                {
                    mGrid.mV(i, 0, k) = 0.0;
                    mGrid.mV(i, ny, k) = 0.0;
                }
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    mGrid.mW(i, j, 0) = 0.0;
                    mGrid.mW(i, j, nz) = 0.0;
                }
        }

        /**
         * 压力投影
         * 求解压力泊松方程，保证速度场无散度（不可压缩）
         */
        void Solver3d::pressureProjection(double dt)
        {
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float h = mGrid.cellSize;
            const float scale = dt / (Eulerian3dPara::airDensity * h * h);

            // 将固体边界的法向速度设为 0
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (mGrid.isSolidCell(i, j, k))
                        {
                            mGrid.mU(i, j, k) = mGrid.mU(i + 1, j, k) = 0;
                            mGrid.mV(i, j, k) = mGrid.mV(i, j + 1, k) = 0;
                            mGrid.mW(i, j, k) = mGrid.mW(i, j, k + 1) = 0;
                        }

            // 计算散度
            Glb::GridData3d div;
            div.initialize(0.0);

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k)) continue;
                        // 处理固体边界
                        float uR = mGrid.isSolidCell(i + 1, j, k) ? 0 : mGrid.mU(i + 1, j, k);
                        float uL = mGrid.isSolidCell(i - 1, j, k) ? 0 : mGrid.mU(i, j, k);
                        float vT = mGrid.isSolidCell(i, j + 1, k) ? 0 : mGrid.mV(i, j + 1, k);
                        float vB = mGrid.isSolidCell(i, j - 1, k) ? 0 : mGrid.mV(i, j, k);
                        float wF = mGrid.isSolidCell(i, j, k + 1) ? 0 : mGrid.mW(i, j, k + 1);
                        float wK = mGrid.isSolidCell(i, j, k - 1) ? 0 : mGrid.mW(i, j, k);
                        div(i, j, k) = (uR - uL + vT - vB + wF - wK) / h;
                    }

            // Gauss-Seidel 迭代求解压力
            Glb::GridData3d p;
            p.initialize(0.0);

            for (int iter = 0; iter < 50; ++iter)
                for (int k = 0; k < nz; ++k)
                    for (int j = 0; j < ny; ++j)
                        for (int i = 0; i < nx; ++i)
                        {
                            if (mGrid.isSolidCell(i, j, k)) continue;
                            // 统计流体邻居数量和压力和
                            int n = 0;
                            float pSum = 0.0f;
                            if (i > 0 && !mGrid.isSolidCell(i - 1, j, k)) { pSum += p(i - 1, j, k); n++; }
                            if (i < nx - 1 && !mGrid.isSolidCell(i + 1, j, k)) { pSum += p(i + 1, j, k); n++; }
                            if (j > 0 && !mGrid.isSolidCell(i, j - 1, k)) { pSum += p(i, j - 1, k); n++; }
                            if (j < ny - 1 && !mGrid.isSolidCell(i, j + 1, k)) { pSum += p(i, j + 1, k); n++; }
                            if (k > 0 && !mGrid.isSolidCell(i, j, k - 1)) { pSum += p(i, j, k - 1); n++; }
                            if (k < nz - 1 && !mGrid.isSolidCell(i, j, k + 1)) { pSum += p(i, j, k + 1); n++; }
                            if (n > 0) p(i, j, k) = (pSum - div(i, j, k) / scale) / n;
                        }

            // 应用压力梯度（只在两侧都是流体时）
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 1; i < nx; ++i)
                        if (!mGrid.isSolidCell(i - 1, j, k) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mU(i, j, k) -= scale * (p(i, j, k) - p(i - 1, j, k)) * h;

            for (int k = 0; k < nz; ++k)
                for (int j = 1; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (!mGrid.isSolidCell(i, j - 1, k) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mV(i, j, k) -= scale * (p(i, j, k) - p(i, j - 1, k)) * h;

            for (int k = 1; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (!mGrid.isSolidCell(i, j, k - 1) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mW(i, j, k) -= scale * (p(i, j, k) - p(i, j, k - 1)) * h;
        }

        /**
         * G2P：网格速度转移回粒子
         * 使用 FLIP/PIC 混合方法更新粒子速度
         */
        void Solver3d::gridToParticle(double dt)
        {
            const float flipRatio = 0.95f;  // FLIP 比例（0=纯PIC，1=纯FLIP）

            for (auto &p : mPs.particles)
            {
                // PIC：直接从新网格插值
                glm::vec3 vPIC = mGrid.getVelocity(p.position);

                // FLIP：粒子速度 + 网格速度变化量
                glm::vec3 vOld = mGrid.getVelocityFromGrid(p.position, mU_prev, mV_prev, mW_prev);
                glm::vec3 vFLIP = p.velocity + (vPIC - vOld);

                // 混合 FLIP 和 PIC
                p.velocity = flipRatio * vFLIP + (1.0f - flipRatio) * vPIC;
            }
        }

        /**
         * 移动粒子
         * 使用子步法防止高速粒子穿透边界
         */
        void Solver3d::advectParticles(double dt)
        {
            const float h = mGrid.cellSize;
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float eps = h * 0.01f;  // 边界偏移量
            const glm::vec3 minP(eps), maxP(nx * h - eps, ny * h - eps, nz * h - eps);

            const float restitution = PIC3dPara::wallRestitution;  // 反弹系数
            const float friction = 0.9f;  // 切向摩擦保留

            for (auto &p : mPs.particles)
            {
                // 计算需要的子步数：确保每步移动不超过半个格子
                float dist = glm::length(p.velocity) * (float)dt;
                int numSubsteps = dist < h * 0.5f ? 1 : (int)std::ceil(dist / (h * 0.5f));
                float subDt = (float)dt / numSubsteps;

                for (int s = 0; s < numSubsteps; ++s)
                {
                    glm::vec3 oldPos = p.position;
                    p.position += subDt * p.velocity;

                    // 域边界碰撞（带反弹）
                    if (p.position.x < minP.x) { p.position.x = minP.x; p.velocity.x *= -restitution; }
                    if (p.position.x > maxP.x) { p.position.x = maxP.x; p.velocity.x *= -restitution; }
                    if (p.position.y < minP.y) { p.position.y = minP.y; p.velocity.y *= -restitution; }
                    if (p.position.y > maxP.y) { p.position.y = maxP.y; p.velocity.y *= -restitution; }
                    if (p.position.z < minP.z) { p.position.z = minP.z; p.velocity.z *= -restitution; }
                    if (p.position.z > maxP.z) { p.position.z = maxP.z; p.velocity.z *= -restitution; }

                    // 固体碰撞检测
                    int ci, cj, ck;
                    if (mGrid.inSolid(p.position, ci, cj, ck))
                    {
                        // 计算碰撞法线（从旧位置指向固体中心的反方向）
                        glm::vec3 solidCenter = mGrid.getCenter(ci, cj, ck);
                        glm::vec3 normal = glm::normalize(oldPos - solidCenter);

                        // 分解速度为法向和切向
                        float vn = glm::dot(p.velocity, normal);
                        glm::vec3 vNormal = vn * normal;
                        glm::vec3 vTangent = p.velocity - vNormal;

                        // 反弹：法向反转并衰减，切向保留（带摩擦）
                        if (vn < 0)  // 只有朝向固体时才反弹
                            p.velocity = -restitution * vNormal + friction * vTangent;

                        // 将粒子推回到固体外
                        p.position = oldPos;
                    }
                }
            }
        }
    }
}

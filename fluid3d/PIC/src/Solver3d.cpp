/**
 * Solver3d.cpp: 3D PIC 求解器实现
 * 实现 PIC/FLIP 混合方法的完整求解流程
 */

#include "PIC/include/Solver3d.h"
#include "../../../common/include/GridData3d.h"
#include <cmath>
#include <vector>

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 构造函数：初始化求解器
         * 成员变量（在 Solver3d 类中）：
         * - mPs: 对 ParticleSystem3d 的引用，用于访问粒子
         * - mGrid: 对 PICGrid3d 的引用，用于访问和修改网格数据
         * - mU_prev/mV_prev/mW_prev: 保存上一帧网格速度以支持 FLIP 更新
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
            const double dt = PIC3dPara::dt;

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

            mGrid.updateSources();

            // 6. 添加外力（浮力）
            addForces(dt);

            // 6.5 涡量增强（可选，vorticityConst=0 时不生效）
            if (PIC3dPara::vorticityConst > 0.0f)
                applyVorticityConfinement(dt);

            // 7. 压力投影（保证不可压缩）
            pressureProjection(dt);

            // 8. G2P：网格速度转移回粒子
            gridToParticle(dt);

            // 9. 移动粒子
            advectParticles(dt);
        }

        // ==================== 涡量增强（Vorticity Confinement, 3D） ====================
        void Solver3d::applyVorticityConfinement(double dt)
        {
            const int nx = mGrid.dim[PICGrid3d::X];
            const int ny = mGrid.dim[PICGrid3d::Y];
            const int nz = mGrid.dim[PICGrid3d::Z];
            const float h = mGrid.cellSize;
            const float eps = PIC3dPara::vorticityConst;
            if (eps <= 0.0f || nx <= 1 || ny <= 1 || nz <= 1)
                return;

            auto idx = [nx, ny](int i, int j, int k) -> int { return i + j * nx + k * nx * ny; };

            // cell-centered velocities
            std::vector<float> uC((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> vC((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> wC((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);

            auto uCenter = [&](int i, int j, int k) -> float {
                float uL = mGrid.isSolidFace(i, j, k, PICGrid3d::X) ? 0.0f : (float)mGrid.mU(i, j, k);
                float uR = mGrid.isSolidFace(i + 1, j, k, PICGrid3d::X) ? 0.0f : (float)mGrid.mU(i + 1, j, k);
                return 0.5f * (uL + uR);
            };
            auto vCenter = [&](int i, int j, int k) -> float {
                float vB = mGrid.isSolidFace(i, j, k, PICGrid3d::Y) ? 0.0f : (float)mGrid.mV(i, j, k);
                float vT = mGrid.isSolidFace(i, j + 1, k, PICGrid3d::Y) ? 0.0f : (float)mGrid.mV(i, j + 1, k);
                return 0.5f * (vB + vT);
            };
            auto wCenter = [&](int i, int j, int k) -> float {
                float wK = mGrid.isSolidFace(i, j, k, PICGrid3d::Z) ? 0.0f : (float)mGrid.mW(i, j, k);
                float wF = mGrid.isSolidFace(i, j, k + 1, PICGrid3d::Z) ? 0.0f : (float)mGrid.mW(i, j, k + 1);
                return 0.5f * (wK + wF);
            };

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k))
                            continue;
                        uC[idx(i, j, k)] = uCenter(i, j, k);
                        vC[idx(i, j, k)] = vCenter(i, j, k);
                        wC[idx(i, j, k)] = wCenter(i, j, k);
                    }

            // curl ω and |ω|
            std::vector<float> wx((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> wy((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> wz((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> mag((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);

            auto sample = [&](const std::vector<float> &a, int i, int j, int k) -> float {
                if (i < 0) i = 0; if (i > nx - 1) i = nx - 1;
                if (j < 0) j = 0; if (j > ny - 1) j = ny - 1;
                if (k < 0) k = 0; if (k > nz - 1) k = nz - 1;
                return a[idx(i, j, k)];
            };

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k))
                            continue;

                        const float dwdy = (j == 0) ? (sample(wC, i, 1, k) - sample(wC, i, 0, k)) / h
                                                   : (j == ny - 1) ? (sample(wC, i, ny - 1, k) - sample(wC, i, ny - 2, k)) / h
                                                                  : (sample(wC, i, j + 1, k) - sample(wC, i, j - 1, k)) / (2.0f * h);
                        const float dvdz = (k == 0) ? (sample(vC, i, j, 1) - sample(vC, i, j, 0)) / h
                                                   : (k == nz - 1) ? (sample(vC, i, j, nz - 1) - sample(vC, i, j, nz - 2)) / h
                                                                  : (sample(vC, i, j, k + 1) - sample(vC, i, j, k - 1)) / (2.0f * h);

                        const float dudz = (k == 0) ? (sample(uC, i, j, 1) - sample(uC, i, j, 0)) / h
                                                   : (k == nz - 1) ? (sample(uC, i, j, nz - 1) - sample(uC, i, j, nz - 2)) / h
                                                                  : (sample(uC, i, j, k + 1) - sample(uC, i, j, k - 1)) / (2.0f * h);
                        const float dwdx = (i == 0) ? (sample(wC, 1, j, k) - sample(wC, 0, j, k)) / h
                                                   : (i == nx - 1) ? (sample(wC, nx - 1, j, k) - sample(wC, nx - 2, j, k)) / h
                                                                  : (sample(wC, i + 1, j, k) - sample(wC, i - 1, j, k)) / (2.0f * h);

                        const float dvdx = (i == 0) ? (sample(vC, 1, j, k) - sample(vC, 0, j, k)) / h
                                                   : (i == nx - 1) ? (sample(vC, nx - 1, j, k) - sample(vC, nx - 2, j, k)) / h
                                                                  : (sample(vC, i + 1, j, k) - sample(vC, i - 1, j, k)) / (2.0f * h);
                        const float dudy = (j == 0) ? (sample(uC, i, 1, k) - sample(uC, i, 0, k)) / h
                                                   : (j == ny - 1) ? (sample(uC, i, ny - 1, k) - sample(uC, i, ny - 2, k)) / h
                                                                  : (sample(uC, i, j + 1, k) - sample(uC, i, j - 1, k)) / (2.0f * h);

                        const float ox = dwdy - dvdz;
                        const float oy = dudz - dwdx;
                        const float oz = dvdx - dudy;

                        const int id = idx(i, j, k);
                        wx[id] = ox;
                        wy[id] = oy;
                        wz[id] = oz;
                        mag[id] = std::sqrt(ox * ox + oy * oy + oz * oz);
                    }

            // force at cell centers
            std::vector<float> fx((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> fy((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
            std::vector<float> fz((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);

            auto magAt = [&](int i, int j, int k) -> float {
                if (i < 0) i = 0; if (i > nx - 1) i = nx - 1;
                if (j < 0) j = 0; if (j > ny - 1) j = ny - 1;
                if (k < 0) k = 0; if (k > nz - 1) k = nz - 1;
                return mag[idx(i, j, k)];
            };

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k))
                            continue;

                        const float dmx = (i == 0) ? (magAt(1, j, k) - magAt(0, j, k)) / h
                                                 : (i == nx - 1) ? (magAt(nx - 1, j, k) - magAt(nx - 2, j, k)) / h
                                                                : (magAt(i + 1, j, k) - magAt(i - 1, j, k)) / (2.0f * h);
                        const float dmy = (j == 0) ? (magAt(i, 1, k) - magAt(i, 0, k)) / h
                                                 : (j == ny - 1) ? (magAt(i, ny - 1, k) - magAt(i, ny - 2, k)) / h
                                                                : (magAt(i, j + 1, k) - magAt(i, j - 1, k)) / (2.0f * h);
                        const float dmz = (k == 0) ? (magAt(i, j, 1) - magAt(i, j, 0)) / h
                                                 : (k == nz - 1) ? (magAt(i, j, nz - 1) - magAt(i, j, nz - 2)) / h
                                                                : (magAt(i, j, k + 1) - magAt(i, j, k - 1)) / (2.0f * h);

                        const float len = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
                        if (len < 1e-8f)
                            continue;

                        const float Nx = dmx / len;
                        const float Ny = dmy / len;
                        const float Nz = dmz / len;

                        const int id = idx(i, j, k);
                        const float ox = wx[id], oy = wy[id], oz = wz[id];

                        // f = eps * h * (N × ω)
                        const float cx = Ny * oz - Nz * oy;
                        const float cy = Nz * ox - Nx * oz;
                        const float cz = Nx * oy - Ny * ox;

                        fx[id] = eps * h * cx;
                        fy[id] = eps * h * cy;
                        fz[id] = eps * h * cz;
                    }

            // add to staggered velocities by averaging adjacent cell-centered forces
            auto cellForce = [&](const std::vector<float> &f, int i, int j, int k) -> float {
                if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz)
                    return 0.0f;
                if (mGrid.isSolidCell(i, j, k))
                    return 0.0f;
                return f[idx(i, j, k)];
            };

            // U faces
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i <= nx; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::X))
                            continue;
                        float fL = cellForce(fx, i - 1, j, k);
                        float fR = cellForce(fx, i, j, k);
                        float fFace = (i == 0) ? fR : (i == nx) ? fL : 0.5f * (fL + fR);
                        mGrid.mU(i, j, k) += (double)dt * (double)fFace;
                    }

            // V faces
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j <= ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::Y))
                            continue;
                        float fB = cellForce(fy, i, j - 1, k);
                        float fT = cellForce(fy, i, j, k);
                        float fFace = (j == 0) ? fT : (j == ny) ? fB : 0.5f * (fB + fT);
                        mGrid.mV(i, j, k) += (double)dt * (double)fFace;
                    }

            // W faces
            for (int k = 0; k <= nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::Z))
                            continue;
                        float fK = cellForce(fz, i, j, k - 1);
                        float fF = cellForce(fz, i, j, k);
                        float fFace = (k == 0) ? fF : (k == nz) ? fK : 0.5f * (fK + fF);
                        mGrid.mW(i, j, k) += (double)dt * (double)fFace;
                    }
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
            newT.initialize(PIC3dPara::ambientTemp);
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

            // // 在 Z 方向添加浮力
            // for (int k = 1; k < nz; ++k)
            //     for (int j = 0; j < ny; ++j)
            //         for (int i = 0; i < nx; ++i)
            //         {
            //             glm::vec3 pos((i + 0.5f) * h, (j + 0.5f) * h, k * h);
            //             double buoy = mGrid.getBoussinesqForce(pos);
            //             mGrid.mW(i, j, k) += dt * buoy;
            //         }

            // wind (X direction)
            // windX > 0: 向 +X；windX < 0: 向 -X
            const double windX = (double)PIC3dPara::windX;
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 1; i < nx; ++i) // U: i=0 与 i=nx 通常是边界面
                        mGrid.mU(i, j, k) += dt * windX;

            // 边界速度设为 0
            // 左侧 X 边界是墙，但“源区域”除外（避免把喷口速度直接抹掉）
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                {
                    bool isSource = false;
                    for (const auto &src : PIC3dPara::source)
                    {
                        const int cy = src.position.y;
                        const int cz = src.position.z;
                        const int r = PIC3dPara::emitterRadius;
                        if (j >= cy - r && j <= cy + r && k >= cz - r && k <= cz + r)
                        {
                            isSource = true;
                            break;
                        }
                    }
                    if (!isSource)
                        mGrid.mU(0, j, k) = 0.0;
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
            const float scale = dt / (PIC3dPara::airDensity * h * h);

            // 1) 固体面法向速度清零（含容器边界面）
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i <= nx; ++i)
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::X))
                            mGrid.mU(i, j, k) = 0.0;

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j <= ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::Y))
                            mGrid.mV(i, j, k) = 0.0;

            for (int k = 0; k <= nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (mGrid.isSolidFace(i, j, k, PICGrid3d::Z))
                            mGrid.mW(i, j, k) = 0.0;

            // 2) 计算散度（基于 solid face 判定，避免 “solid cell 邻居” 导致的不一致）
            Glb::GridData3d div;
            div.initialize(0.0);

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        if (mGrid.isSolidCell(i, j, k))
                            continue;

                        const float uR = mGrid.isSolidFace(i + 1, j, k, PICGrid3d::X) ? 0.0f : (float)mGrid.mU(i + 1, j, k);
                        const float uL = mGrid.isSolidFace(i, j, k, PICGrid3d::X) ? 0.0f : (float)mGrid.mU(i, j, k);
                        const float vT = mGrid.isSolidFace(i, j + 1, k, PICGrid3d::Y) ? 0.0f : (float)mGrid.mV(i, j + 1, k);
                        const float vB = mGrid.isSolidFace(i, j, k, PICGrid3d::Y) ? 0.0f : (float)mGrid.mV(i, j, k);
                        const float wF = mGrid.isSolidFace(i, j, k + 1, PICGrid3d::Z) ? 0.0f : (float)mGrid.mW(i, j, k + 1);
                        const float wK = mGrid.isSolidFace(i, j, k, PICGrid3d::Z) ? 0.0f : (float)mGrid.mW(i, j, k);

                        div(i, j, k) = (uR - uL + vT - vB + wF - wK) / h;
                    }

            // 3) Gauss-Seidel 迭代求解压力
            Glb::GridData3d p;
            p.initialize(0.0);

            const int iters = (PIC3dPara::pressureIters > 0) ? PIC3dPara::pressureIters : 0;
            for (int iter = 0; iter < iters; ++iter)
                for (int k = 0; k < nz; ++k)
                    for (int j = 0; j < ny; ++j)
                        for (int i = 0; i < nx; ++i)
                        {
                            if (mGrid.isSolidCell(i, j, k))
                                continue;

                            int n = 0;
                            float pSum = 0.0f;
                            if (i > 0 && !mGrid.isSolidCell(i - 1, j, k)) { pSum += (float)p(i - 1, j, k); ++n; }
                            if (i < nx - 1 && !mGrid.isSolidCell(i + 1, j, k)) { pSum += (float)p(i + 1, j, k); ++n; }
                            if (j > 0 && !mGrid.isSolidCell(i, j - 1, k)) { pSum += (float)p(i, j - 1, k); ++n; }
                            if (j < ny - 1 && !mGrid.isSolidCell(i, j + 1, k)) { pSum += (float)p(i, j + 1, k); ++n; }
                            if (k > 0 && !mGrid.isSolidCell(i, j, k - 1)) { pSum += (float)p(i, j, k - 1); ++n; }
                            if (k < nz - 1 && !mGrid.isSolidCell(i, j, k + 1)) { pSum += (float)p(i, j, k + 1); ++n; }

                            if (n > 0)
                                p(i, j, k) = (pSum - (float)div(i, j, k) / scale) / (float)n;
                        }

            // 4) 应用压力梯度到面速度（仅当面两侧都是流体 cell 时）
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 1; i < nx; ++i)
                        if (!mGrid.isSolidCell(i - 1, j, k) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mU(i, j, k) -= scale * ((float)p(i, j, k) - (float)p(i - 1, j, k)) * h;

            for (int k = 0; k < nz; ++k)
                for (int j = 1; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (!mGrid.isSolidCell(i, j - 1, k) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mV(i, j, k) -= scale * ((float)p(i, j, k) - (float)p(i, j - 1, k)) * h;

            for (int k = 1; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        if (!mGrid.isSolidCell(i, j, k - 1) && !mGrid.isSolidCell(i, j, k))
                            mGrid.mW(i, j, k) -= scale * ((float)p(i, j, k) - (float)p(i, j, k - 1)) * h;

            // 5) 出流边界（x = nx*h）：对 u 施加零法向梯度（简单拷贝）
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    mGrid.mU(nx, j, k) = mGrid.mU(nx - 1, j, k);
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

          std::vector<Particle> survivors;
          survivors.reserve(mPs.particles.size());

          const float h = mGrid.cellSize;
          const int nx = mGrid.dim[PICGrid3d::X];
          const int ny = mGrid.dim[PICGrid3d::Y];
          const int nz = mGrid.dim[PICGrid3d::Z];
          const float eps = h * 0.01f; // 边界偏移量
          const glm::vec3 minP(eps), maxP(nx * h - eps, ny * h - eps, nz * h - eps);
          const float outflowX = nx * h; // 右侧出流平面

          const float restitution = PIC3dPara::wallRestitution; // 反弹系数
          const float friction = 0.9f;                          // 切向摩擦保留

                    for (auto &p : mPs.particles)
          {
                        bool alive = true;

            // 计算需要的子步数：确保每步移动不超过半个格子
            float dist = glm::length(p.velocity) * (float)dt;
            int numSubsteps = dist < h * 0.5f ? 1 : (int)std::ceil(dist / (h * 0.5f));
            float subDt = (float)dt / numSubsteps;

                        for (int s = 0; s < numSubsteps; ++s)
            {
              glm::vec3 oldPos = p.position;
              p.position += subDt * p.velocity;

                            // 域边界碰撞（带反弹）
              if (p.position.x < minP.x)
              {
                p.position.x = minP.x;
                p.velocity.x *= -restitution;
              }
                            // 右侧为出流：超过出流面则直接移除粒子
                            if (p.position.x >= outflowX)
                            {
                                alive = false;
                                break;
                            }
              if (p.position.y < minP.y)
              {
                p.position.y = minP.y;
                p.velocity.y *= -restitution;
              }
              if (p.position.y > maxP.y)
              {
                p.position.y = maxP.y;
                p.velocity.y *= -restitution;
              }
              if (p.position.z < minP.z)
              {
                p.position.z = minP.z;
                p.velocity.z *= -restitution;
              }
              if (p.position.z > maxP.z)
              {
                p.position.z = maxP.z;
                p.velocity.z *= -restitution;
              }

                            // 中心球体碰撞（SDF 法线，更圆）
                            if (Eulerian3dPara::addSolid && mGrid.inSphere(p.position))
                            {
                                glm::vec3 n = mGrid.sphereNormal(p.position);
                                float vn = glm::dot(p.velocity, n);
                                glm::vec3 vN = vn * n;
                                glm::vec3 vT = p.velocity - vN;

                                if (vn < 0.0f)
                                    p.velocity = -restitution * vN + friction * vT;

                                p.position = mGrid.projectOutOfSphere(p.position, eps);
                            }
            }

                        if (!alive)
                            continue;

            survivors.push_back(p);
          }

            // 用幸存者替换粒子列表
            mPs.particles.swap(survivors);
        }
    }
}

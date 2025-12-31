/**
 * Solver.cpp: 3D欧拉流体求解器实现文件
 * 实现流体仿真的核心算法
 */

#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"
#include <algorithm>

namespace FluidSimulation
{
    namespace Eulerian3d
    {
        /**
         * 构造函数，初始化求解器并重置网格
         * @param grid MAC网格引用
         */
        Solver::Solver(MACGrid3d &grid) : mGrid(grid)
        {
            // 初始化时重置网格
            mGrid.reset();
        }

        /**
         * 求解流体方程
         * 实现一步3D流体仿真的主要步骤
         */
        void Solver::solve()
        {
            const double dt = Eulerian3dPara::dt;
            const double invRho = 1.0 / Eulerian3dPara::airDensity;

            const float h = mGrid.cellSize;
            const float h2 = h * h;

            // ------------------------------------------------------------------
            // 1) 平流：半拉格朗日对速度、密度、温度进行平流
            // ------------------------------------------------------------------
            Glb::GridData3dX advU;
            Glb::GridData3dY advV;
            Glb::GridData3dZ advW;
            Glb::CubicGridData3d advD;
            Glb::CubicGridData3d advT;

            advU.initialize(0.0);
            advV.initialize(0.0);
            advW.initialize(0.0);
            advD.initialize(0.0);
            advT.initialize(Eulerian3dPara::ambientTemp);

            // u 分量：位置 (i*h, (j+0.5)h, (k+0.5)h)
            for (int k = 0; k < mGrid.dim[MACGrid3d::Z]; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y]; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X] + 1; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::X))
                        {
                            advU(i, j, k) = 0.0;
                            continue;
                        }

                        glm::vec3 facePos((float)i * h, ((float)j + 0.5f) * h, ((float)k + 0.5f) * h);
                        glm::vec3 backPos = mGrid.semiLagrangian(facePos, dt);
                        advU(i, j, k) = mGrid.getVelocityX(backPos);
                    }
                }
            }

            // v 分量：位置 ((i+0.5)h, j*h, (k+0.5)h)
            for (int k = 0; k < mGrid.dim[MACGrid3d::Z]; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y] + 1; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X]; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::Y))
                        {
                            advV(i, j, k) = 0.0;
                            continue;
                        }

                        glm::vec3 facePos(((float)i + 0.5f) * h, (float)j * h, ((float)k + 0.5f) * h);
                        glm::vec3 backPos = mGrid.semiLagrangian(facePos, dt);
                        advV(i, j, k) = mGrid.getVelocityY(backPos);
                    }
                }
            }

            // w 分量：位置 ((i+0.5)h, (j+0.5)h, k*h)
            for (int k = 0; k < mGrid.dim[MACGrid3d::Z] + 1; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y]; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X]; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z))
                        {
                            advW(i, j, k) = 0.0;
                            continue;
                        }

                        glm::vec3 facePos(((float)i + 0.5f) * h, ((float)j + 0.5f) * h, (float)k * h);
                        glm::vec3 backPos = mGrid.semiLagrangian(facePos, dt);
                        advW(i, j, k) = mGrid.getVelocityZ(backPos);
                    }
                }
            }

            // 标量场
            FOR_EACH_CELL
            {
                glm::vec3 center = mGrid.getCenter(i, j, k);
                glm::vec3 backPos = mGrid.semiLagrangian(center, dt);
                advD(i, j, k) = mGrid.getDensity(backPos);
                advT(i, j, k) = mGrid.getTemperature(backPos);
            }

            mGrid.mU = advU;
            mGrid.mV = advV;
            mGrid.mW = advW;
            mGrid.mD = advD;
            mGrid.mT = advT;

            // ------------------------------------------------------------------
            // 2) 外力：Boussinesq 浮力（沿 Z 方向，作用在 w 分量）
            // ------------------------------------------------------------------
            for (int k = 0; k < mGrid.dim[MACGrid3d::Z] + 1; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y]; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X]; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z))
                        {
                            mGrid.mW(i, j, k) = 0.0;
                            continue;
                        }

                        glm::vec3 pos(((float)i + 0.5f) * h, ((float)j + 0.5f) * h, (float)k * h);
                        double buoy = mGrid.getBoussinesqForce(pos);
                        mGrid.mW(i, j, k) += dt * buoy;
                    }
                }
            }

            // ------------------------------------------------------------------
            // 3) 压力求解（Jacobi）：?? p = (ρ/dt) ?·u
            // ------------------------------------------------------------------
            Glb::GridData3d divergence;
            divergence.initialize(0.0);

            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j, k))
                {
                    divergence(i, j, k) = 0.0;
                    continue;
                }
                divergence(i, j, k) = mGrid.getDivergence(i, j, k);
            }

            Glb::GridData3d pressure;
            Glb::GridData3d pressureNew;
            pressure.initialize(0.0);
            pressureNew.initialize(0.0);

            const int maxIter = 60;
            const double rhsScale = Eulerian3dPara::airDensity / dt;

            for (int iter = 0; iter < maxIter; ++iter)
            {
                FOR_EACH_CELL
                {
                    if (mGrid.isSolidCell(i, j, k))
                    {
                        pressureNew(i, j, k) = 0.0;
                        continue;
                    }

                    double rhs = rhsScale * divergence(i, j, k);
                    double sum = 0.0;
                    int count = 0;

                    if (!mGrid.isSolidCell(i - 1, j, k))
                    {
                        sum += pressure(i - 1, j, k);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i + 1, j, k))
                    {
                        sum += pressure(i + 1, j, k);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j - 1, k))
                    {
                        sum += pressure(i, j - 1, k);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j + 1, k))
                    {
                        sum += pressure(i, j + 1, k);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j, k - 1))
                    {
                        sum += pressure(i, j, k - 1);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j, k + 1))
                    {
                        sum += pressure(i, j, k + 1);
                        ++count;
                    }

                    if (count == 0)
                    {
                        pressureNew(i, j, k) = 0.0;
                    }
                    else
                    {
                        pressureNew(i, j, k) = (sum - rhs * h2) / (double)count;
                    }
                }

                pressure = pressureNew;
            }

            // ------------------------------------------------------------------
            // 4) 投影：从速度中减去压力梯度
            // ------------------------------------------------------------------
            for (int k = 0; k < mGrid.dim[MACGrid3d::Z]; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y]; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X] + 1; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::X))
                        {
                            mGrid.mU(i, j, k) = 0.0;
                            continue;
                        }

                        double pR = mGrid.isSolidCell(i, j, k) ? 0.0 : pressure(i, j, k);
                        double pL = mGrid.isSolidCell(i - 1, j, k) ? 0.0 : pressure(i - 1, j, k);
                        mGrid.mU(i, j, k) -= dt * invRho * (pR - pL) / h;
                    }
                }
            }

            for (int k = 0; k < mGrid.dim[MACGrid3d::Z]; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y] + 1; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X]; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::Y))
                        {
                            mGrid.mV(i, j, k) = 0.0;
                            continue;
                        }

                        double pT = mGrid.isSolidCell(i, j, k) ? 0.0 : pressure(i, j, k);
                        double pB = mGrid.isSolidCell(i, j - 1, k) ? 0.0 : pressure(i, j - 1, k);
                        mGrid.mV(i, j, k) -= dt * invRho * (pT - pB) / h;
                    }
                }
            }

            for (int k = 0; k < mGrid.dim[MACGrid3d::Z] + 1; ++k)
            {
                for (int j = 0; j < mGrid.dim[MACGrid3d::Y]; ++j)
                {
                    for (int i = 0; i < mGrid.dim[MACGrid3d::X]; ++i)
                    {
                        if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z))
                        {
                            mGrid.mW(i, j, k) = 0.0;
                            continue;
                        }

                        double pF = mGrid.isSolidCell(i, j, k) ? 0.0 : pressure(i, j, k);
                        double pB = mGrid.isSolidCell(i, j, k - 1) ? 0.0 : pressure(i, j, k - 1);
                        mGrid.mW(i, j, k) -= dt * invRho * (pF - pB) / h;
                    }
                }
            }
        }
    }
}

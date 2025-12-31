/**
 * Solver.cpp: 2D欧拉流体求解器实现文件
 * 实现流体仿真的核心算法
 */

#include "Eulerian/include/Solver.h"
#include "Configure.h"
#include <algorithm>

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        /**
         * 构造函数，初始化求解器并重置网格
         * @param grid MAC网格引用
         */
        Solver::Solver(MACGrid2d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        /**
         * 求解流体方程
         * 实现一步流体仿真的主要步骤
         */
        void Solver::solve()
        {
            const double dt = Eulerian2dPara::dt;
            const double invRho = 1.0 / Eulerian2dPara::airDensity;

            const float h = mGrid.cellSize;
            const float h2 = h * h;

            // ------------------------------------------------------------------
            // 1) 平流：使用半拉格朗日法对速度、温度、密度进行平流
            // ------------------------------------------------------------------
            Glb::GridData2dX advU;
            Glb::GridData2dY advV;
            Glb::CubicGridData2d advD;
            Glb::CubicGridData2d advT;

            advU.initialize(0.0);
            advV.initialize(0.0);
            advD.initialize(0.0);
            advT.initialize(Eulerian2dPara::ambientTemp);

            // 平流速度场（u 存在垂直网格线中点，v 存在水平网格线中点）
            for (int j = 0; j < mGrid.dim[MACGrid2d::Y]; ++j)
            {
                for (int i = 0; i < mGrid.dim[MACGrid2d::X] + 1; ++i)
                {
                    if (mGrid.isSolidFace(i, j, MACGrid2d::X))
                    {
                        advU(i, j) = 0.0;
                        continue;
                    }

                    glm::vec2 facePos((float)i * h, ((float)j + 0.5f) * h);
                    glm::vec2 backPos = mGrid.semiLagrangian(facePos, dt);
                    advU(i, j) = mGrid.getVelocityX(backPos);
                }
            }

            for (int j = 0; j < mGrid.dim[MACGrid2d::Y] + 1; ++j)
            {
                for (int i = 0; i < mGrid.dim[MACGrid2d::X]; ++i)
                {
                    if (mGrid.isSolidFace(i, j, MACGrid2d::Y))
                    {
                        advV(i, j) = 0.0;
                        continue;
                    }

                    glm::vec2 facePos(((float)i + 0.5f) * h, (float)j * h);
                    glm::vec2 backPos = mGrid.semiLagrangian(facePos, dt);
                    advV(i, j) = mGrid.getVelocityY(backPos);
                }
            }

            // 平流标量场（密度、温度）
            FOR_EACH_CELL
            {
                glm::vec2 center = mGrid.getCenter(i, j);
                glm::vec2 backPos = mGrid.semiLagrangian(center, dt);
                advD(i, j) = mGrid.getDensity(backPos);
                advT(i, j) = mGrid.getTemperature(backPos);
            }

            mGrid.mU = advU;
            mGrid.mV = advV;
            mGrid.mD = advD;
            mGrid.mT = advT;

            // ------------------------------------------------------------------
            // 2) 外力：Boussinesq 浮力（作用在 v 分量，即 Y 方向）
            // ------------------------------------------------------------------
            for (int j = 0; j < mGrid.dim[MACGrid2d::Y] + 1; ++j)
            {
                for (int i = 0; i < mGrid.dim[MACGrid2d::X]; ++i)
                {
                    if (mGrid.isSolidFace(i, j, MACGrid2d::Y))
                    {
                        mGrid.mV(i, j) = 0.0;
                        continue;
                    }

                    glm::vec2 pos(((float)i + 0.5f) * h, (float)j * h);
                    double buoy = mGrid.getBoussinesqForce(pos);
                    mGrid.mV(i, j) += dt * buoy;
                }
            }

            // ------------------------------------------------------------------
            // 3) 求解压力泊松方程（Jacobi 迭代）
            //    ?? p = (ρ/dt) ?·u
            // ------------------------------------------------------------------
            Glb::GridData2d divergence;
            divergence.initialize(0.0);

            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j))
                {
                    divergence(i, j) = 0.0;
                    continue;
                }
                divergence(i, j) = mGrid.getDivergence(i, j);
            }

            Glb::GridData2d pressure;
            Glb::GridData2d pressureNew;
            pressure.initialize(0.0);
            pressureNew.initialize(0.0);

            const int maxIter = 60;
            const double rhsScale = Eulerian2dPara::airDensity / dt;

            for (int iter = 0; iter < maxIter; ++iter)
            {
                FOR_EACH_CELL
                {
                    if (mGrid.isSolidCell(i, j))
                    {
                        pressureNew(i, j) = 0.0;
                        continue;
                    }

                    double rhs = rhsScale * divergence(i, j);
                    double sum = 0.0;
                    int count = 0;

                    if (!mGrid.isSolidCell(i - 1, j))
                    {
                        sum += pressure(i - 1, j);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i + 1, j))
                    {
                        sum += pressure(i + 1, j);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j - 1))
                    {
                        sum += pressure(i, j - 1);
                        ++count;
                    }
                    if (!mGrid.isSolidCell(i, j + 1))
                    {
                        sum += pressure(i, j + 1);
                        ++count;
                    }

                    if (count == 0)
                    {
                        pressureNew(i, j) = 0.0;
                    }
                    else
                    {
                        pressureNew(i, j) = (sum - rhs * h2) / (double)count;
                    }
                }

                pressure = pressureNew;
            }

            // ------------------------------------------------------------------
            // 4) 投影：速度减去压力梯度，使速度场无散
            // ------------------------------------------------------------------
            for (int j = 0; j < mGrid.dim[MACGrid2d::Y]; ++j)
            {
                for (int i = 0; i < mGrid.dim[MACGrid2d::X] + 1; ++i)
                {
                    if (mGrid.isSolidFace(i, j, MACGrid2d::X))
                    {
                        mGrid.mU(i, j) = 0.0;
                        continue;
                    }

                    double pR = mGrid.isSolidCell(i, j) ? 0.0 : pressure(i, j);
                    double pL = mGrid.isSolidCell(i - 1, j) ? 0.0 : pressure(i - 1, j);
                    mGrid.mU(i, j) -= dt * invRho * (pR - pL) / h;
                }
            }

            for (int j = 0; j < mGrid.dim[MACGrid2d::Y] + 1; ++j)
            {
                for (int i = 0; i < mGrid.dim[MACGrid2d::X]; ++i)
                {
                    if (mGrid.isSolidFace(i, j, MACGrid2d::Y))
                    {
                        mGrid.mV(i, j) = 0.0;
                        continue;
                    }

                    double pT = mGrid.isSolidCell(i, j) ? 0.0 : pressure(i, j);
                    double pB = mGrid.isSolidCell(i, j - 1) ? 0.0 : pressure(i, j - 1);
                    mGrid.mV(i, j) -= dt * invRho * (pT - pB) / h;
                }
            }
        }
    }
}

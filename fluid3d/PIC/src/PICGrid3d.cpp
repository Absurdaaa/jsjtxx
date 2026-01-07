/**
 * PICGrid3d.cpp: 3D PIC 网格实现
 * 实现：
 * - 从指定速度网格插值获取速度（FLIP 支持）
 * - 3D 场景：中心球体障碍 + 右侧出流边界（X 正方向）
 */

#include "PIC/include/PICGrid3d.h"
#include "GridData3d.h"
#include <cmath>
namespace FluidSimulation
{
    namespace PIC3d
    {
        PICGrid3d::PICGrid3d() : Eulerian3d::MACGrid3d()
        {
            // 重新走一遍 PIC3d 的固体初始化（基类构造会先调用基类 createSolids）
            initialize();
        }

        void PICGrid3d::initialize()
        {
            reset();

            const int nx = dim[Eulerian3d::MACGrid3d::X];
            const int ny = dim[Eulerian3d::MACGrid3d::Y];
            const int nz = dim[Eulerian3d::MACGrid3d::Z];
            const float h = cellSize;

            // 场景：中心球体
            mSphereCenter = glm::vec3(nx * 0.5f * h, ny * 0.5f * h, nz * 0.5f * h);
            const int minDim = (nx < ny) ? ((nx < nz) ? nx : nz) : ((ny < nz) ? ny : nz);

            float rCells = PIC3dPara::sphereRadiusCells;
            if (rCells < 0.0f)
                rCells = 0.0f;
            const float maxRC = 0.49f * (float)minDim;
            if (rCells > maxRC)
                rCells = maxRC;
            // 若 rCells == 0，则等价于“关闭障碍物”
            mSphereRadius = rCells * h;

            createSolids();
        }

        void PICGrid3d::createSolids()
        {
            mSolid.initialize(0.0);

            if (!Eulerian3dPara::addSolid)
                return;

            const int nx = dim[Eulerian3d::MACGrid3d::X];
            const int ny = dim[Eulerian3d::MACGrid3d::Y];
            const int nz = dim[Eulerian3d::MACGrid3d::Z];

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        const glm::vec3 c = getCenter(i, j, k);
                        if (sphereSDF(c) < 0.0f)
                            mSolid(i, j, k) = 1;
                    }
        }

        int PICGrid3d::isSolidFace(int i, int j, int k, Direction d)
        {
            const int nx = dim[Eulerian3d::MACGrid3d::X];
            const int ny = dim[Eulerian3d::MACGrid3d::Y];
            const int nz = dim[Eulerian3d::MACGrid3d::Z];
            const float h = cellSize;

            // 容器边界：右侧 X 边界(i==nx)作为出流，不视为固体面
            if (d == X && (i == 0))
                return 1;
            if (d == X && i == nx)
                return 0;
            if (d == Y && (j == 0 || j == ny))
                return 1;
            if (d == Z && (k == 0 || k == nz))
                return 1;

            // 中心球体：用 face-center 的 SDF 判定更圆
            glm::vec3 faceCenter(0.0f);
            if (d == X)
                faceCenter = glm::vec3(i * h, (j + 0.5f) * h, (k + 0.5f) * h);
            else if (d == Y)
                faceCenter = glm::vec3((i + 0.5f) * h, j * h, (k + 0.5f) * h);
            else
                faceCenter = glm::vec3((i + 0.5f) * h, (j + 0.5f) * h, k * h);

            if (Eulerian3dPara::addSolid && sphereSDF(faceCenter) < 0.0f)
                return 1;

            return 0;
        }

        void PICGrid3d::updateSources()
        {
            const int nx = dim[Eulerian3d::MACGrid3d::X];
            const int ny = dim[Eulerian3d::MACGrid3d::Y];
            const int nz = dim[Eulerian3d::MACGrid3d::Z];

            const int radius = PIC3dPara::emitterRadius;

            for (const auto &src : PIC3dPara::source)
            {
                const int cx = src.position.x;
                const int cy = src.position.y;
                const int cz = src.position.z;

                // 做一个薄的 x 厚度（类似 2D 的 dx=-1..1），在 y-z 平面做圆盘覆盖
                for (int dx = -1; dx <= 1; ++dx)
                {
                    for (int dy = -radius; dy <= radius; ++dy)
                    {
                        for (int dz = -radius; dz <= radius; ++dz)
                        {
                            if (radius > 0 && (dy * dy + dz * dz > radius * radius))
                                continue;

                            const int x = cx + dx;
                            const int y = cy + dy;
                            const int z = cz + dz;
                            if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
                                continue;
                            if (isSolidCell(x, y, z))
                                continue;

                            mT(x, y, z) = src.temp;
                            mD(x, y, z) = src.density;
                            mU(x, y, z) = src.velocity.x;
                            mV(x, y, z) = src.velocity.y;
                            mW(x, y, z) = src.velocity.z;
                        }
                    }
                }
            }
        }

        float PICGrid3d::sphereSDF(const glm::vec3 &pt) const
        {
            return glm::length(pt - mSphereCenter) - mSphereRadius;
        }

        bool PICGrid3d::inSphere(const glm::vec3 &pt) const
        {
            return sphereSDF(pt) < 0.0f;
        }

        glm::vec3 PICGrid3d::sphereNormal(const glm::vec3 &pt) const
        {
            glm::vec3 n = pt - mSphereCenter;
            float len = glm::length(n);
            if (len < 1e-8f)
                return glm::vec3(1.0f, 0.0f, 0.0f);
            return n / len;
        }

        glm::vec3 PICGrid3d::projectOutOfSphere(const glm::vec3 &pt, float extra) const
        {
            const float d = sphereSDF(pt);
            if (d >= 0.0f)
                return pt;
            const glm::vec3 n = sphereNormal(pt);
            return pt + n * (-d + extra);
        }

        /**
         * 从指定的速度网格插值获取速度
         * 使用三线性插值从 MAC 网格的交错位置获取速度分量
         * 该类主要使用父类 MACGrid3d 中的成员（如 mU/mV/mW, cellSize 等）
         */
        glm::vec3 PICGrid3d::getVelocityFromGrid(const glm::vec3 &pt,
                                                  Glb::GridData3dX &uGrid,
                                                  Glb::GridData3dY &vGrid,
                                                  Glb::GridData3dZ &wGrid)
        {
            float u = uGrid.interpolate(pt);  // X方向速度（在YZ面上）
            float v = vGrid.interpolate(pt);  // Y方向速度（在XZ面上）
            float w = wGrid.interpolate(pt);  // Z方向速度（在XY面上）
            return glm::vec3(u, v, w);
        }
    }
}

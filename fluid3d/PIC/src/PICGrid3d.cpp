/**
 * PICGrid3d.cpp: 3D PIC 网格实现
 * 实现从指定速度网格插值获取速度的功能
 */

#include "PIC/include/PICGrid3d.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
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

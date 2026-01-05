/**
 * PICGrid3d.h: 3D PIC 网格头文件
 * 继承自 MACGrid3d，提供 PIC 方法特有的网格操作
 */

#pragma once
#ifndef __PIC_GRID_3D_H__
#define __PIC_GRID_3D_H__

#include "MACGrid3d.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 3D PIC 网格类
         * 继承自 MACGrid3d，复用其速度场、温度场、密度场等数据结构
         * 该类主要复用父类的成员变量（在 MACGrid3d 中定义），例如：
         * - mU/mV/mW: 网格在三个叉格面上的速度分量
         * - mT/mD: 温度和密度的标量场
         * - dim: 网格尺寸 (x,y,z)
         * - cellSize: 网格单元大小
         */
        class PICGrid3d : public Eulerian3d::MACGrid3d
        {
        public:
            PICGrid3d() : Eulerian3d::MACGrid3d() {}

            /**
             * 从指定的速度网格插值获取速度
             * 用于 FLIP 方法中获取旧速度场的值
             * @param pt 查询点的世界坐标
             * @param uGrid X方向速度网格
             * @param vGrid Y方向速度网格
             * @param wGrid Z方向速度网格
             * @return 插值得到的速度向量
             */
            glm::vec3 getVelocityFromGrid(const glm::vec3 &pt,
                                          Glb::GridData3dX &uGrid,
                                          Glb::GridData3dY &vGrid,
                                          Glb::GridData3dZ &wGrid);
        };

/** 遍历所有网格单元的宏 */
#define PIC3D_FOR_EACH_CELL                                              \
    for (int k = 0; k < Eulerian3dPara::theDim3d[MACGrid3d::Z]; k++)     \
        for (int j = 0; j < Eulerian3dPara::theDim3d[MACGrid3d::Y]; j++) \
            for (int i = 0; i < Eulerian3dPara::theDim3d[MACGrid3d::X]; i++)
    }
}

#endif

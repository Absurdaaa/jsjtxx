#pragma once
#ifndef __PIC_GRID_2D_H__
#define __PIC_GRID_2D_H__

#include <windows.h>
#include <glm/glm.hpp>
#include "GridData2d.h"
#include <Logger.h>
#include "Configure.h"

namespace FluidSimulation
{
    namespace PIC2d
    {
      
        
            /**
             * 2D PIC 网格类
             * 存储速度、密度、温度、固体信息，提供采样、几何、边界等操作
             */
            class PICGrid2d
            {
            public:
                /**
                 * 网格方向枚举（X/Y）
                 */
                enum Direction
                {
                    X, ///< X方向
                    Y  ///< Y方向
                };

                PICGrid2d(); ///< 构造函数，初始化网格
                ~PICGrid2d(); ///< 析构函数
                PICGrid2d(const PICGrid2d &orig); ///< 拷贝构造
                PICGrid2d &operator=(const PICGrid2d &orig); ///< 赋值

                /**
                 * 重置所有网格数据为初始值
                 */
                void reset();

                // ----------- 初始化与场景设置 -----------

                /**
                 * 初始化网格（重置、创建固体、断言散度）
                 */
                void initialize();
                /**
                 * 创建固体障碍物
                 */
                void createSolids();
                /**
                 * 根据源参数更新温度/密度/速度
                 */
                void updateSources();

                // ----------- 物理采样与辅助 -----------

                /**
                 * 半拉格朗日追踪辅助
                 * @param pt 采样点
                 * @param dt 时间步长
                 * @return 回溯后的位置
                 */
                glm::vec2 semiLagrangian(const glm::vec2 &pt, double dt);

                // ----------- 物理量采样 -----------
                //
                // 从指定的交错网格（uGrid/vGrid）中在位置 pt 插值得到速度（用于 FLIP 的旧网格查询）
                // 参数 uGrid/vGrid 采用非常规引用以允许访问 operator()(i,j)
                glm::vec2 getVelocityFromGrid(const glm::vec2 &pt, Glb::GridData2dX &uGrid, Glb::GridData2dY &vGrid);
                /**
                 * 采样点速度
                 */
                glm::vec2 getVelocity(const glm::vec2 &pt);
                /**
                 * 采样点X方向速度
                 */
                double getVelocityX(const glm::vec2 &pt);
                /**
                 * 采样点Y方向速度
                 */
                double getVelocityY(const glm::vec2 &pt);
                /**
                 * 采样点温度
                 */
                double getTemperature(const glm::vec2 &pt);
                /**
                 * 采样点密度
                 */
                double getDensity(const glm::vec2 &pt);

                // ----------- 网格几何与索引 -----------

                /**
                 * 获取(i,j)单元中心坐标
                 */
                glm::vec2 getCenter(int i, int j);
                glm::vec2 getLeft(int i, int j);
                glm::vec2 getRight(int i, int j);
                glm::vec2 getTop(int i, int j);
                glm::vec2 getBottom(int i, int j);

                /**
                 * 由一维索引获取(i,j)
                 */
                void getCell(int index, int &i, int &j);
                /**
                 * 由(i,j)获取一维索引
                 */
                int getIndex(int i, int j);
                /**
                 * 判断(i0,j0)与(i1,j1)是否为邻居
                 */
                bool isNeighbor(int i0, int j0, int i1, int j1);
                /**
                 * 判断(i,j)在指定方向上是否合法
                 */
                bool isValid(int i, int j, Direction d);

                // ----------- 边界与固体 -----------

                /**
                 * 判断(i,j)是否为固体单元
                 */
                int isSolidCell(int i, int j);
                /**
                 * 判断(i,j)在d方向是否为固体面
                 */
                int isSolidFace(int i, int j, Direction d);
                /**
                 * 判断点pt是否在固体内
                 */
                bool inSolid(const glm::vec2 &pt);
                /**
                 * 判断点pt是否在固体内，并返回(i,j)
                 */
                bool inSolid(const glm::vec2 &pt, int &i, int &j);
                /**
                 * 判断射线是否与(i,j)单元相交
                 */
                bool intersects(const glm::vec2 &pt, const glm::vec2 &dir, int i, int j, double &time);
                /**
                 * 统计固体单元数
                 */
                int numSolidCells();

                // ----------- 压力/散度 -----------

                /**
                 * 获取(i0,j0)-(i1,j1)之间的压力系数
                 */
                double getPressureCoeffBetweenCells(int i0, int j0, int i1, int j1);
                /**
                 * 计算(i,j)处的散度
                 */
                double getDivergence(int i, int j);
                /**
                 * 检查(i,j)处的散度
                 */
                double checkDivergence(int i, int j);
                /**
                 * 检查全局散度
                 */
                bool checkDivergence();

                // ----------- 物理力 -----------

                /**
                 * 采样点的Boussinesq浮力
                 */
                double getBoussinesqForce(const glm::vec2 &pt);

                // ----------- 渲染辅助 -----------

                /**
                 * 获取(i,j)单元的渲染颜色
                 */
                glm::vec4 getRenderColor(int i, int j);
                /**
                 * 获取任意点的渲染颜色
                 */
                glm::vec4 getRenderColor(const glm::vec2 &pt);  
            
                // ----------- 网格参数与数据 -----------
                ///< 网格单元大小
                float cellSize; 
                ///< 网格维度（x/y方向）
                int dim[2];     
            
                Glb::GridData2dX mU;      ///< X方向速度场
                Glb::GridData2dY mV;      ///< Y方向速度场
                Glb::CubicGridData2d mD;  ///< 密度场
                Glb::CubicGridData2d mT;  ///< 温度场
                Glb::GridData2d mSolid;   ///< 固体标记场
            };

// macros
#define PIC_FOR_EACH_CELL                                                \
    for (int j = 0; j < Eulerian2dPara::theDim2d[PICGrid2d::Y]; j++)     \
        for (int i = 0; i < Eulerian2dPara::theDim2d[PICGrid2d::X]; i++)

    }
}

#endif

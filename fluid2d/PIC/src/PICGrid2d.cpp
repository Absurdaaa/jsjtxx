#include "PIC/include/PICGrid2d.h"
#include <math.h>
#include <algorithm>
#include <assert.h>

namespace FluidSimulation
{
    namespace PIC2d
    {
        // 构造函数：根据全局参数初始化网格尺寸和单元大小，并初始化所有网格数据
        PICGrid2d::PICGrid2d()
        {
            cellSize = PIC2dPara::theCellSize2d; // 单元大小
            dim[0] = PIC2dPara::theDim2d[0];     // x方向格点数
            dim[1] = PIC2dPara::theDim2d[1];     // y方向格点数
            initialize(); // 初始化所有网格数据
        }

        // 拷贝构造函数：深拷贝所有网格数据和参数
        PICGrid2d::PICGrid2d(const PICGrid2d &orig)
        {
            mU = orig.mU;      // X方向速度
            mV = orig.mV;      // Y方向速度
            mD = orig.mD;      // 密度
            mT = orig.mT;      // 温度
            mSolid = orig.mSolid; // 固体标记
            cellSize = orig.cellSize;
            dim[0] = orig.dim[0];
            dim[1] = orig.dim[1];
        }

        // 赋值操作符：深拷贝所有网格数据和参数
        PICGrid2d &PICGrid2d::operator=(const PICGrid2d &orig)
        {
            if (&orig == this)
                return *this;
            mU = orig.mU;
            mV = orig.mV;
            mD = orig.mD;
            mT = orig.mT;
            mSolid = orig.mSolid;
            cellSize = orig.cellSize;
            dim[0] = orig.dim[0];
            dim[1] = orig.dim[1];
            return *this;
        }

        // 析构函数：无特殊资源需释放
        PICGrid2d::~PICGrid2d() {}

        // 重置所有网格数据为初始值
        void PICGrid2d::reset()
        {
            mU.initialize(0.0); // X方向速度清零
            mV.initialize(0.0); // Y方向速度清零
            mD.initialize(0.0); // 密度清零
            mT.initialize(PIC2dPara::ambientTemp); // 温度设为环境温度
        }

        // 创建固体障碍物（如中间一条横线）
        void PICGrid2d::createSolids()
        {
            mSolid.initialize(); // 全部初始化为流体
            if (PIC2dPara::addSolid)
            {
                // 将固体设置为一个圆形障碍，位于域中心稍偏上的位置（空中球体效果）
                int cx = dim[0] / 2; // 圆心（格子索引）
                int cy = dim[1] / 2; // 放在中间高度
                // 半径按网格高度的一定比例，例如总高度的1/6
                int radius = dim[1] / 10;
                if (radius < 2) radius = 2;
                for (int j = 0; j < dim[1]; ++j)
                {
                    for (int i = 0; i < dim[0]; ++i)
                    {
                        int dx = i - cx;
                        int dy = j - cy;
                        if (dx * dx + dy * dy <= radius * radius)
                            mSolid(i, j) = 1; // 设置为固体
                    }
                }
            }
        }

        // 根据源参数更新温度/密度/速度
        void PICGrid2d::updateSources()
        {
            for (int i = 0; i < PIC2dPara::source.size(); i++)
            {
                int x = PIC2dPara::source[i].position.x;
                int y = PIC2dPara::source[i].position.y;
                mT(x, y) = PIC2dPara::source[i].temp;         // 源温度
                mD(x, y) = PIC2dPara::source[i].density;      // 源密度
                mU(x, y) = PIC2dPara::source[i].velocity.x;   // 源X速度
                mV(x, y) = PIC2dPara::source[i].velocity.y;   // 源Y速度
            }
        }

        // 初始化网格：重置、创建固体、断言散度
        void PICGrid2d::initialize()
        {
            reset();
            createSolids();
            assert(checkDivergence()); // 检查初始散度
        }

        // 计算Boussinesq浮力（温度/密度驱动的Y向力）, 
        double PICGrid2d::getBoussinesqForce(const glm::vec2 &pos)
        {
            double temperature = getTemperature(pos);
            double smokeDensity = getDensity(pos);
            double yforce = -PIC2dPara::boussinesqAlpha * smokeDensity +
                            PIC2dPara::boussinesqBeta * (temperature - PIC2dPara::ambientTemp);
            return yforce;
        }

        // 计算(i,j)处的速度散度
        double PICGrid2d::getDivergence(int i, int j)
        {
            double x1 = isSolidCell(i + 1, j) ? 0.0 : mU(i + 1, j);
            double x0 = isSolidCell(i - 1, j) ? 0.0 : mU(i, j);
            double y1 = isSolidCell(i, j + 1) ? 0.0 : mV(i, j + 1);
            double y0 = isSolidCell(i, j - 1) ? 0.0 : mV(i, j);
            double div = (x1 - x0 + y1 - y0) / cellSize;
            return div;
        }

        // 检查(i,j)处的散度（不考虑固体）
        double PICGrid2d::checkDivergence(int i, int j)
        {
            double x1 = mU(i + 1, j);
            double x0 = mU(i, j);
            double y1 = mV(i, j + 1);
            double y0 = mV(i, j);
            double div = (x1 - x0 + y1 - y0) / cellSize;
            return div;
        }

        // 检查全局散度是否足够小
        bool PICGrid2d::checkDivergence()
        {
            PIC_FOR_EACH_CELL
            {
                double div = checkDivergence(i, j);
                if (fabs(div) > 0.01)
                    return false;
            }
            return true;
        }

        // 半拉格朗日追踪：回溯粒子/标量到上一步位置
        glm::vec2 PICGrid2d::semiLagrangian(const glm::vec2 &pt, double dt)
        {
            glm::vec2 vel = getVelocity(pt); // 当前速度
            glm::vec2 pos = pt - vel * (float)dt; // 回溯
            // 限制在域内
            pos[0] = max(0.0, min((dim[0] - 1) * cellSize, pos[0]));
            pos[1] = max(0.0, min((dim[1] - 1) * cellSize, pos[1]));

            int i, j;
            if (inSolid(pt, i, j))
            {
                double t = 0;
                if (intersects(pt, vel, i, j, t))
                    pos = pt - vel * (float)t; // 穿越固体时修正
                else
                    Glb::Logger::getInstance().addLog("Error: something goes wrong during advection");
            }
            return pos;
        }

        // 判断射线是否与(i,j)单元相交（用于半拉格朗日/边界处理）
        bool PICGrid2d::intersects(const glm::vec2 &wPos, const glm::vec2 &wDir, int i, int j, double &time)
        {
            glm::vec2 pos = getCenter(i, j);
            glm::vec2 rayStart = wPos - pos;
            glm::vec2 rayDir = wDir;
            double tmin = -1e18, tmax = 1e18;
            double minv = -0.5 * cellSize;
            double maxv = 0.5 * cellSize;
            for (int k = 0; k < 3; k++)
            {
                double e = rayStart[k];
                double f = rayDir[k];
                if (fabs(f) > 1e-9)
                {
                    double t1 = (minv - e) / f;
                    double t2 = (maxv - e) / f;
                    if (t1 > t2)
                        std::swap(t1, t2);
                    tmin = tmin > t1 ? tmin : t1;
                    tmax = tmax < t2 ? tmax : t2;
                    if (tmin > tmax || tmax < 0)
                        return false;
                }
                else if (e < minv || e > maxv)
                    return false;
            }
            time = tmin >= 0 ? tmin : tmax;
            return true;
        }

        // (i,j) 转一维索引，越界返回-1
        int PICGrid2d::getIndex(int i, int j)
        {
            if (i < 0 || i > dim[0] - 1)
                return -1;
            if (j < 0 || j > dim[1] - 1)
                return -1;
            return i + j * dim[0];
        }

        // 一维索引转(i,j)
        void PICGrid2d::getCell(int index, int &i, int &j)
        {
            j = (int)index / dim[0];
            i = index - j * dim[0];
        }

        // 获取(i,j)单元中心坐标
        glm::vec2 PICGrid2d::getCenter(int i, int j)
        {
            double x = cellSize * 0.5 + i * cellSize;
            double y = cellSize * 0.5 + j * cellSize;
            return glm::vec2(x, y);
        }

        // 获取(i,j)单元各边界点坐标
        glm::vec2 PICGrid2d::getLeft(int i, int j) { return getCenter(i, j) - glm::vec2(cellSize * 0.5, 0.0); }
        glm::vec2 PICGrid2d::getRight(int i, int j) { return getCenter(i, j) + glm::vec2(cellSize * 0.5, 0.0); }
        glm::vec2 PICGrid2d::getTop(int i, int j) { return getCenter(i, j) + glm::vec2(0.0, cellSize * 0.5); }
        glm::vec2 PICGrid2d::getBottom(int i, int j) { return getCenter(i, j) - glm::vec2(0.0, cellSize * 0.5); }

        // 采样点速度（若在固体内返回0）
        glm::vec2 PICGrid2d::getVelocity(const glm::vec2 &pt)
        {
            if (inSolid(pt))
                return glm::vec2(0);
            glm::vec2 vel;
            vel[0] = getVelocityX(pt);
            vel[1] = getVelocityY(pt);
            return vel;
        }

        glm::vec2 PICGrid2d::getVelocityFromGrid(const glm::vec2 &pt, Glb::GridData2dX &uGrid, Glb::GridData2dY &vGrid)
        {
            const float h = cellSize;
            // U 分量在竖直面：(i * h, (j + 0.5) * h)
            double uSum = 0.0;
            double wUSum = 0.0;
            {
                double fx = pt.x / h;
                double fy = pt.y / h - 0.5;
                int i0 = (int)floor(fx);
                int j0 = (int)floor(fy);
                double rx = fx - i0;
                double ry = fy - j0;
                for (int di = 0; di <= 1; ++di)
                {
                    for (int dj = 0; dj <= 1; ++dj)
                    {
                        int i = i0 + di;
                        int j = j0 + dj;
                        if (!isValid(i, j, PICGrid2d::X) || isSolidFace(i, j, PICGrid2d::X))
                            continue;
                        double w = (di ? rx : (1.0 - rx)) * (dj ? ry : (1.0 - ry));
                        uSum += uGrid(i, j) * w;
                        wUSum += w;
                    }
                }
            }

            // V 分量在水平面：((i + 0.5) * h, j * h)
            double vSum = 0.0;
            double wVSum = 0.0;
            {
                double fx = pt.x / h - 0.5;
                double fy = pt.y / h;
                int i0 = (int)floor(fx);
                int j0 = (int)floor(fy);
                double rx = fx - i0;
                double ry = fy - j0;
                for (int di = 0; di <= 1; ++di)
                {
                    for (int dj = 0; dj <= 1; ++dj)
                    {
                        int i = i0 + di;
                        int j = j0 + dj;
                        if (!isValid(i, j, PICGrid2d::Y) || isSolidFace(i, j, PICGrid2d::Y))
                            continue;
                        double w = (di ? rx : (1.0 - rx)) * (dj ? ry : (1.0 - ry));
                        vSum += vGrid(i, j) * w;
                        wVSum += w;
                    }
                }
            }

            double u = (wUSum > 1e-12) ? (uSum / wUSum) : 0.0;
            double v = (wVSum > 1e-12) ? (vSum / wVSum) : 0.0;
            return glm::vec2((float)u, (float)v);
        }

        // 采样点X方向速度
        double PICGrid2d::getVelocityX(const glm::vec2 &pt) { return mU.interpolate(pt); }
        // 采样点Y方向速度
        double PICGrid2d::getVelocityY(const glm::vec2 &pt) { return mV.interpolate(pt); }
        // 采样点温度
        double PICGrid2d::getTemperature(const glm::vec2 &pt) { return mT.interpolate(pt); }
        // 采样点密度
        double PICGrid2d::getDensity(const glm::vec2 &pt) { return mD.interpolate(pt); }

        // 统计固体单元数
        int PICGrid2d::numSolidCells()
        {
            int numSolid = 0;
            PIC_FOR_EACH_CELL { numSolid += mSolid(i, j); }
            return numSolid;
        }

        // 判断点pt是否在固体内
        bool PICGrid2d::inSolid(const glm::vec2 &pt)
        {
            int i, j;
            mSolid.getCell(pt, i, j);
            return isSolidCell(i, j) == 1;
        }

        // 判断点pt是否在固体内，并返回(i,j)
        bool PICGrid2d::inSolid(const glm::vec2 &pt, int &i, int &j)
        {
            mSolid.getCell(pt, i, j);
            return isSolidCell(i, j) == 1;
        }

        // 判断(i,j)是否为固体单元（边界或障碍物）
        int PICGrid2d::isSolidCell(int i, int j)
        {
            bool boundary = (i < 0 || i > dim[0] - 1 || j < 0 || j > dim[1] - 1);
            bool object = (mSolid(i, j) == 1);
            return boundary || object ? 1 : 0;
        }

        // 判断(i,j)在d方向是否为固体面
        int PICGrid2d::isSolidFace(int i, int j, PICGrid2d::Direction d)
        {
            if (d == X && (i == 0 || i == dim[0]))
                return 1;
            else if (d == Y && (j == 0 || j == dim[1]))
                return 1;
            if (d == X && (mSolid(i, j) || mSolid(i - 1, j)))
                return 1;
            if (d == Y && (mSolid(i, j) || mSolid(i, j - 1)))
                return 1;
            return 0;
        }

        // 判断(i0,j0)与(i1,j1)是否为邻居
        bool PICGrid2d::isNeighbor(int i0, int j0, int i1, int j1)
        {
            if (abs(i0 - i1) == 1 && j0 == j1)
                return true;
            if (abs(j0 - j1) == 1 && i0 == i1)
                return true;
            return false;
        }

        // 获取(i,j)-(pi,pj)之间的压力系数（用于压力方程）
        double PICGrid2d::getPressureCoeffBetweenCells(int i, int j, int pi, int pj)
        {
            if (i == pi && j == pj)
            {
                int numSolidNeighbors = isSolidCell(i + 1, j) + isSolidCell(i - 1, j) + isSolidCell(i, j + 1) + isSolidCell(i, j - 1);
                return 4.0 - numSolidNeighbors;
            }
            if (isNeighbor(i, j, pi, pj) && !isSolidCell(pi, pj))
                return -1.0;
            return 0.0;
        }

        // 获取(i,j)单元的渲染颜色（可自定义）
        glm::vec4 PICGrid2d::getRenderColor(int i, int j)
        {
            double value = mD(i, j);
            return glm::vec4(1.0, 1.0, 1.0, value);
        }

        // 获取任意点的渲染颜色（可自定义）
        glm::vec4 PICGrid2d::getRenderColor(const glm::vec2 &pt)
        {
            double value = getDensity(pt);
            return glm::vec4(value, value, value, value);
        }

        // 判断(i,j)在指定方向上是否合法
        bool PICGrid2d::isValid(int i, int j, PICGrid2d::Direction d)
        {
            switch (d)
            {
            case X:
                return (i >= 0 && i < dim[X] + 1 && j >= 0 && j < dim[Y]);
            case Y:
                return (i >= 0 && i < dim[X] && j >= 0 && j < dim[Y] + 1);
            }
            Glb::Logger::getInstance().addLog("Error: bad direction");
            return false;
        }
    }
}

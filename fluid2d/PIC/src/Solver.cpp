#include "PIC/include/Solver.h"
#include "GridData2d.h"
#include <cmath>

namespace FluidSimulation
{
  namespace PIC2d
  {
    Solver::Solver(ParticleSystem &ps, PICGrid2d &grid)
        : mPs(ps), mGrid(grid)
    {
      mU_prev.initialize(0.0);
      mV_prev.initialize(0.0);
    }

    void Solver::solve()
    {
      const double dt = PIC2dPara::dt;

      // 1. 发射新粒子
      mPs.emitFromSources();

      // 2. 更新网格源（温度/密度）
      mGrid.updateSources();

      // 3. 平流温度和密度场
      advectScalars(dt);

      // 4. P2G: 粒子速度转移到网格
      particleToGrid();

      // 5. 保存旧速度（用于FLIP）
      mU_prev = mGrid.mU;
      mV_prev = mGrid.mV;

      // 6. 添加外力（浮力）
      addForces(dt);

      // 7. 压力投影（保证不可压缩）
      pressureProjection(dt);

      // 8. G2P: 网格速度转移回粒子
      gridToParticle(dt);

      // 9. 移动粒子
      advectParticles(dt);
    }

    // ==================== P2G ====================
    void Solver::particleToGrid()
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;

      // 清零速度场
      mGrid.mU.initialize(0.0);
      mGrid.mV.initialize(0.0);

      // 权重累加器
      Glb::GridData2dX wU;
      wU.initialize(0.0);
      Glb::GridData2dY wV;
      wV.initialize(0.0);

      for (const auto &p : mPs.particles)
      {
        // === U 分量 (位于垂直面中心: i*h, (j+0.5)*h) ===
        {
          float fx = p.position.x / h;
          float fy = p.position.y / h - 0.5f;
          int i0 = (int)std::floor(fx);
          int j0 = (int)std::floor(fy);
          float tx = fx - i0;
          float ty = fy - j0;

          for (int di = 0; di <= 1; ++di)
          {
            for (int dj = 0; dj <= 1; ++dj)
            {
              int i = i0 + di;
              int j = j0 + dj;
              if (i < 0 || i > nx || j < 0 || j >= ny)
                continue;
              float w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty);
              mGrid.mU(i, j) += w * p.velocity.x;
              wU(i, j) += w;
            }
          }
        }

        // === V 分量 (位于水平面中心: (i+0.5)*h, j*h) ===
        {
          float fx = p.position.x / h - 0.5f;
          float fy = p.position.y / h;
          int i0 = (int)std::floor(fx);
          int j0 = (int)std::floor(fy);
          float tx = fx - i0;
          float ty = fy - j0;

          for (int di = 0; di <= 1; ++di)
          {
            for (int dj = 0; dj <= 1; ++dj)
            {
              int i = i0 + di;
              int j = j0 + dj;
              if (i < 0 || i >= nx || j < 0 || j > ny)
                continue;
              float w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty);
              mGrid.mV(i, j) += w * p.velocity.y;
              wV(i, j) += w;
            }
          }
        }
      }

      // 归一化
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
          if (wU(i, j) > 1e-6f)
            mGrid.mU(i, j) /= wU(i, j);

      for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i)
          if (wV(i, j) > 1e-6f)
            mGrid.mV(i, j) /= wV(i, j);
    }

    // ==================== 标量平流 ====================
    void Solver::advectScalars(double dt)
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;

      Glb::CubicGridData2d newT, newD;
      newT.initialize(PIC2dPara::ambientTemp);
      newD.initialize(0.0);

      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidCell(i, j)) continue;

          // 单元中心位置
          glm::vec2 pos((i + 0.5f) * h, (j + 0.5f) * h);
          // 半拉格朗日回溯
          glm::vec2 backPos = mGrid.semiLagrangian(pos, dt);
          // 插值获取旧值
          newT(i, j) = mGrid.getTemperature(backPos);
          newD(i, j) = mGrid.getDensity(backPos);
        }
      }

      // 逐元素复制
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          mGrid.mT(i, j) = newT(i, j);
          mGrid.mD(i, j) = newD(i, j);
        }
      }
    }

    // ==================== 外力 ====================
    void Solver::addForces(double dt)
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;

      // 只添加浮力（烟雾模拟）
      // for (int j = 1; j < ny; ++j)  // 跳过边界
      // {
      //   for (int i = 0; i < nx; ++i)
      //   {
      //     glm::vec2 pos((i + 0.5f) * h, j * h);
      //     double buoy = mGrid.getBoussinesqForce(pos);
      //     mGrid.mV(i, j) += dt * buoy;
      //   }
      // }
      

      // 边界速度：左/上/下是墙；右侧是出口（不再把右边界 u 设为 0）
      for (int i = 0; i < nx; ++i)
      {
        mGrid.mV(i, 0) = 0.0;
        mGrid.mV(i, ny) = 0.0;
      }
      for (int j = 0; j < ny; ++j)
      {
        mGrid.mU(0, j) = 0.0;
      }
    }

    // ==================== 压力投影 ====================
    void Solver::pressureProjection(double dt)
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;
      const float scale = dt / (PIC2dPara::airDensity * h * h);

      // 先将 solid faces 的法向速度设为 0（使用 isSolidFace 的“连续圆形边界”判定）
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i <= nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::X))
            mGrid.mU(i, j) = 0.0;
        }
      }
      for (int j = 0; j <= ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::Y))
            mGrid.mV(i, j) = 0.0;
        }
      }

      // 计算散度（用 solid face 判定，避免靠 solid cell 邻居导致边界误判；右侧出口可自然产生通量）
      Glb::GridData2d div;
      div.initialize(0.0);

      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidCell(i, j))
            continue;

          float uRight = mGrid.isSolidFace(i + 1, j, PICGrid2d::X) ? 0.0f : (float)mGrid.mU(i + 1, j);
          float uLeft  = mGrid.isSolidFace(i,     j, PICGrid2d::X) ? 0.0f : (float)mGrid.mU(i,     j);
          float vTop   = mGrid.isSolidFace(i, j + 1, PICGrid2d::Y) ? 0.0f : (float)mGrid.mV(i, j + 1);
          float vBottom= mGrid.isSolidFace(i, j,     PICGrid2d::Y) ? 0.0f : (float)mGrid.mV(i, j);

          div(i, j) = (uRight - uLeft + vTop - vBottom) / h;
        }
      }

      // Gauss-Seidel 迭代求解压力
      Glb::GridData2d p;
      p.initialize(0.0);

      for (int iter = 0; iter < 50; ++iter)
      {
        for (int j = 0; j < ny; ++j)
        {
          for (int i = 0; i < nx; ++i)
          {
            if (mGrid.isSolidCell(i, j))
              continue;

            int numFluidNeighbors = 0;
            float pSum = 0.0f;

            if (i > 0 && !mGrid.isSolidCell(i - 1, j)) { pSum += p(i - 1, j); numFluidNeighbors++; }
            if (i < nx - 1 && !mGrid.isSolidCell(i + 1, j)) { pSum += p(i + 1, j); numFluidNeighbors++; }
            if (j > 0 && !mGrid.isSolidCell(i, j - 1)) { pSum += p(i, j - 1); numFluidNeighbors++; }
            if (j < ny - 1 && !mGrid.isSolidCell(i, j + 1)) { pSum += p(i, j + 1); numFluidNeighbors++; }

            if (numFluidNeighbors > 0)
              p(i, j) = (pSum - div(i, j) / scale) / numFluidNeighbors;
          }
        }
      }

      // 应用压力梯度（只在两侧都是流体时）
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 1; i < nx; ++i)
        {
          if (!mGrid.isSolidCell(i - 1, j) && !mGrid.isSolidCell(i, j))
            mGrid.mU(i, j) -= scale * (p(i, j) - p(i - 1, j)) * h;
        }
      }

      for (int j = 1; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (!mGrid.isSolidCell(i, j - 1) && !mGrid.isSolidCell(i, j))
            mGrid.mV(i, j) -= scale * (p(i, j) - p(i, j - 1)) * h;
        }
      }
    }

    // ==================== G2P ====================
    void Solver::gridToParticle(double dt)
    {
      const float flipRatio = 0.95f;
      const float h = mGrid.cellSize;
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];

      for (auto &p : mPs.particles)
      {
        // PIC: 直接从新网格插值
        glm::vec2 vPIC = mGrid.getVelocity(p.position);

        // FLIP: 粒子速度 + 网格速度变化量
        glm::vec2 vOld = mGrid.getVelocityFromGrid(p.position, mU_prev, mV_prev);
        glm::vec2 vFLIP = p.velocity + (vPIC - vOld);

        // 混合
        p.velocity = flipRatio * vFLIP + (1.0f - flipRatio) * vPIC;
      }
    }

    // ==================== 粒子移动 ====================
    void Solver::advectParticles(double dt)
    {
      const float h = mGrid.cellSize;
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float minX = h * 0.01f;
      const float maxX = nx * h - h * 0.01f;
      const float minY = h * 0.01f;
      const float maxY = ny * h - h * 0.01f;

      const float restitution = PIC2dPara::wallRestitution; // 反弹系数
      const float friction = 0.9f; // 切向摩擦保留

        std::vector<Particle> survivors;
        survivors.reserve(mPs.particles.size());

        for (auto &p : mPs.particles)
      {
        // 计算需要的子步数：确保每步移动不超过半个格子
        float dist = glm::length(p.velocity) * (float)dt;
        int numSubsteps;
        if (dist < h * 0.5f)
          numSubsteps = 1;
        else
          numSubsteps = (int)std::ceil(dist / (h * 0.5f));
        float subDt = (float)dt / numSubsteps;

        for (int s = 0; s < numSubsteps; ++s)
        {
          glm::vec2 oldPos = p.position;
          p.position += subDt * p.velocity;

          // 域边界碰撞（带反弹）
          if (p.position.x < minX) { p.position.x = minX; p.velocity.x *= -restitution; }
            if (p.position.x > maxX) { 
              // 标记为到达出口，放在域外以便后续删除
              p.position.x = maxX; 
              p.velocity.x = 0.0f; 
            }
          if (p.position.y < minY) { p.position.y = minY; p.velocity.y *= -restitution; }
          if (p.position.y > maxY) { p.position.y = maxY; p.velocity.y *= -restitution; }

          // 固体碰撞检测（圆形障碍：用连续法线与投影，边界更圆滑）
          if (mGrid.inSolid(p.position))
          {
            glm::vec2 normal = mGrid.getSolidNormal(p.position);

            // 分解速度为法向和切向
            float vn = glm::dot(p.velocity, normal);
            glm::vec2 vNormal = vn * normal;
            glm::vec2 vTangent = p.velocity - vNormal;

            // 反弹：法向反转并衰减，切向保留（带摩擦）
            if (vn < 0) // 只有朝向固体时才反弹
              p.velocity = -restitution * vNormal + friction * vTangent;

            // 将粒子推回到固体外
            p.position = mGrid.projectOutOfSolid(oldPos, 0.25f * h);
          }
        }
        
          // 如果粒子到达右边界（出口），则不加入 survivors，表示消失
          if (p.position.x >= maxX - 1e-6f)
          {
            // skip (particle disappears)
            continue;
          }

          survivors.push_back(p);
      }

        // 替换粒子列表为幸存者
        mPs.particles.swap(survivors);

      }
  } // namespace PIC2d
} // namespace FluidSimulation

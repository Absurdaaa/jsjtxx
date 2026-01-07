#include "PIC/include/Solver.h"
#include "../../../common/include/GridData2d.h"
#include <cmath>
#include <vector>

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

      // 4.5 重新设置源速度（P2G会清零速度场）
      mGrid.updateSources();

      // 5. 保存旧速度（用于FLIP）
      mU_prev = mGrid.mU;
      mV_prev = mGrid.mV;

      // 6. 添加外力（浮力）
      addForces(dt);

      // 6.5 涡量增强（可选，vorticityConst=0 时不生效）
      if (PIC2dPara::vorticityConst > 0.0f)
        applyVorticityConfinement(dt);

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

      // 右侧为出口：将最右一列标量清空，避免“边界残留”在流场变化时表现为持续产烟
      // （半拉格朗日平流在开边界处不守恒，最简单稳定的 outflow 处理就是做一个吸收边界）
      const int iExit = nx - 1;
      if (iExit >= 0)
      {
        for (int j = 0; j < ny; ++j)
        {
          mGrid.mT(iExit, j) = PIC2dPara::ambientTemp;
          mGrid.mD(iExit, j) = 0.0;
        }
      }
    }

    // ==================== 涡量增强（Vorticity Confinement） ====================
    void Solver::applyVorticityConfinement(double dt)
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;
      const float eps = PIC2dPara::vorticityConst;
      if (eps <= 0.0f || nx <= 1 || ny <= 1) return;

      auto idx = [nx](int i, int j) { return i + j * nx; };
      std::vector<float> omega(nx * ny, 0.0f);
      std::vector<float> mag(nx * ny, 0.0f);
      std::vector<float> fx(nx * ny, 0.0f);
      std::vector<float> fy(nx * ny, 0.0f);

      auto uCenter = [&](int i, int j) -> float {
        // 单元中心 (i+0.5, j+0.5) 的 u：平均左右两个竖直面速度
        float uL = (mGrid.isSolidFace(i, j, PICGrid2d::X) ? 0.0f : (float)mGrid.mU(i, j));
        float uR = (mGrid.isSolidFace(i + 1, j, PICGrid2d::X) ? 0.0f : (float)mGrid.mU(i + 1, j));
        return 0.5f * (uL + uR);
      };
      auto vCenter = [&](int i, int j) -> float {
        // 单元中心 (i+0.5, j+0.5) 的 v：平均上下两个水平面速度
        float vB = (mGrid.isSolidFace(i, j, PICGrid2d::Y) ? 0.0f : (float)mGrid.mV(i, j));
        float vT = (mGrid.isSolidFace(i, j + 1, PICGrid2d::Y) ? 0.0f : (float)mGrid.mV(i, j + 1));
        return 0.5f * (vB + vT);
      };

      // 1) 计算单元中心涡量 ω = dv/dx - du/dy
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidCell(i, j))
            continue;

          float dvdx;
          if (i == 0) dvdx = (vCenter(1, j) - vCenter(0, j)) / h;
          else if (i == nx - 1) dvdx = (vCenter(nx - 1, j) - vCenter(nx - 2, j)) / h;
          else dvdx = (vCenter(i + 1, j) - vCenter(i - 1, j)) / (2.0f * h);

          float dudy;
          if (j == 0) dudy = (uCenter(i, 1) - uCenter(i, 0)) / h;
          else if (j == ny - 1) dudy = (uCenter(i, ny - 1) - uCenter(i, ny - 2)) / h;
          else dudy = (uCenter(i, j + 1) - uCenter(i, j - 1)) / (2.0f * h);

          float w = dvdx - dudy;
          omega[idx(i, j)] = w;
          mag[idx(i, j)] = std::fabs(w);
        }
      }

      // 2) 计算 ∇|ω|，并得到 confinement 力 f = eps * h * (N × ω)
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidCell(i, j))
            continue;

          auto magAt = [&](int ii, int jj) -> float {
            if (ii < 0) ii = 0;
            if (ii > nx - 1) ii = nx - 1;
            if (jj < 0) jj = 0;
            if (jj > ny - 1) jj = ny - 1;
            return mag[idx(ii, jj)];
          };

          float Nx;
          if (i == 0) Nx = (magAt(1, j) - magAt(0, j)) / h;
          else if (i == nx - 1) Nx = (magAt(nx - 1, j) - magAt(nx - 2, j)) / h;
          else Nx = (magAt(i + 1, j) - magAt(i - 1, j)) / (2.0f * h);

          float Ny;
          if (j == 0) Ny = (magAt(i, 1) - magAt(i, 0)) / h;
          else if (j == ny - 1) Ny = (magAt(i, ny - 1) - magAt(i, ny - 2)) / h;
          else Ny = (magAt(i, j + 1) - magAt(i, j - 1)) / (2.0f * h);

          float len = std::sqrt(Nx * Nx + Ny * Ny);
          if (len < 1e-6f)
            continue;
          Nx /= len;
          Ny /= len;

          float w = omega[idx(i, j)];
          // N × (ω k) = (Ny * ω, -Nx * ω)
          fx[idx(i, j)] = eps * h * (Ny * w);
          fy[idx(i, j)] = eps * h * (-Nx * w);
        }
      }

      // 3) 将中心力平均到 MAC 面上并加到速度场（投影前）
      // U faces: i = 0..nx, j = 0..ny-1
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i <= nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::X))
            continue;

          float f;
          if (i == 0) f = fx[idx(0, j)];
          else if (i == nx) f = fx[idx(nx - 1, j)];
          else f = 0.5f * (fx[idx(i - 1, j)] + fx[idx(i, j)]);

          mGrid.mU(i, j) += (double)dt * (double)f;
        }
      }

      // V faces: i = 0..nx-1, j = 0..ny
      for (int j = 0; j <= ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::Y))
            continue;

          float f;
          if (j == 0) f = fy[idx(i, 0)];
          else if (j == ny) f = fy[idx(i, ny - 1)];
          else f = 0.5f * (fy[idx(i, j - 1)] + fy[idx(i, j)]);

          mGrid.mV(i, j) += (double)dt * (double)f;
        }
      }
    }

    // ==================== 外力 ====================
    void Solver::addForces(double dt)
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;

      // wind
      // windX > 0: 向右；windX < 0: 向左
      double windX = PIC2dPara::windX;
      for (int j = 0; j < ny; ++j)
        for (int i = 1; i < nx; ++i) // U 通常 i=0 和 i=nx 是边界面
          mGrid.mU(i, j) += dt * windX;

      // 边界速度：左/上/下是墙；右侧是出口（不再把右边界 u 设为 0）
      for (int i = 0; i < nx; ++i)
      {
        mGrid.mV(i, 0) = 0.0;
        mGrid.mV(i, ny) = 0.0;
      }
      // 左边界设为墙，但源区域除外
      for (int j = 0; j < ny; ++j)
      {
        bool isSource = false;
        for (const auto &src : PIC2dPara::source)
        {
          int cy = src.position.y;
          int r = PIC2dPara::emitterRadius;
          if (j >= cy - r && j <= cy + r)
          {
            isSource = true;
            break;
          }
        }
        if (!isSource)
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

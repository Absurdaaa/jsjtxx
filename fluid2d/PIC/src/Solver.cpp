#include "PIC/include/Solver.h"
#include "GridData2d.h"
#include <cmath>
#include <iostream>

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
      const double dt = Eulerian2dPara::dt;

      // 1. 发射粒子
      mPs.emitFromSources();

      // 2. 更新欧拉源（温度/密度）
      mGrid.updateSources();

      // 3. 备份网格速度（用于 FLIP）
      mU_prev = mGrid.mU;
      mV_prev = mGrid.mV;

      // 4. P2G
      particleToGrid();

      // 5. 外力（重力 + 浮力）
      addForces(dt);

      // 6. 压力投影（浮力作为 RHS）
      pressureProjection(dt);

      // 7. G2P（PIC/FLIP）
      gridToParticle(dt);

      // 8. Advect
      advectParticles(dt);
    }

    /* ===================== P2G ===================== */

    void Solver::particleToGrid()
    {
      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];
      const float h = mGrid.cellSize;

      mGrid.mU.initialize(0.0);
      mGrid.mV.initialize(0.0);

      Glb::GridData2dX wU;
      wU.initialize(0.0);
      Glb::GridData2dY wV;
      wV.initialize(0.0);

      for (const auto &p : mPs.particles)
      {
        // U
        {
          float fx = p.position.x / h;
          float fy = p.position.y / h - 0.5f;
          int i0 = (int)floor(fx);
          int j0 = (int)floor(fy);
          float rx = fx - i0;
          float ry = fy - j0;

          for (int di = 0; di <= 1; ++di)
            for (int dj = 0; dj <= 1; ++dj)
            {
              int i = i0 + di;
              int j = j0 + dj;
              if (!mGrid.isValid(i, j, PICGrid2d::X) || mGrid.isSolidFace(i, j, PICGrid2d::X))
                continue;
              float w = (di ? rx : 1 - rx) * (dj ? ry : 1 - ry);
              mGrid.mU(i, j) += w * p.velocity.x;
              wU(i, j) += w;
            }
        }

        // V
        {
          float fx = p.position.x / h - 0.5f;
          float fy = p.position.y / h;
          int i0 = (int)floor(fx);
          int j0 = (int)floor(fy);
          float rx = fx - i0;
          float ry = fy - j0;

          for (int di = 0; di <= 1; ++di)
            for (int dj = 0; dj <= 1; ++dj)
            {
              int i = i0 + di;
              int j = j0 + dj;
              if (!mGrid.isValid(i, j, PICGrid2d::Y) || mGrid.isSolidFace(i, j, PICGrid2d::Y))
                continue;
              float w = (di ? rx : 1 - rx) * (dj ? ry : 1 - ry);
              mGrid.mV(i, j) += w * p.velocity.y;
              wV(i, j) += w;
            }
        }
      }

      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx + 1; ++i)
          mGrid.mU(i, j) = (wU(i, j) > 1e-8) ? mGrid.mU(i, j) / wU(i, j) : 0.0;

      for (int j = 0; j < ny + 1; ++j)
        for (int i = 0; i < nx; ++i)
          mGrid.mV(i, j) = (wV(i, j) > 1e-8) ? mGrid.mV(i, j) / wV(i, j) : 0.0;
    }

    /* ===================== 外力 ===================== */

    void Solver::addForces(double dt)
    {
      const double gy = Lagrangian2dPara::gravityY;
      const double invRho = 1.0 / Eulerian2dPara::airDensity;

      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];

      for (int j = 0; j < ny + 1; ++j)
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::Y))
          {
            mGrid.mV(i, j) = 0.0;
            continue;
          }

          glm::vec2 pos((i + 0.5f) * mGrid.cellSize, j * mGrid.cellSize);
          double buoyAccel = mGrid.getBoussinesqForce(pos) * invRho;

          mGrid.mV(i, j) += dt * (gy + buoyAccel);
        }
    }

    /* ===================== Projection ===================== */

    void Solver::pressureProjection(double dt)
    {
      const float h = mGrid.cellSize;
      const float h2 = h * h;
      const double rho = Eulerian2dPara::airDensity;

      const int nx = mGrid.dim[PICGrid2d::X];
      const int ny = mGrid.dim[PICGrid2d::Y];

      Glb::GridData2d divergence;
      divergence.initialize(0.0);

      PIC_FOR_EACH_CELL
      {
        if (mGrid.isSolidCell(i, j))
          continue;

        double div = mGrid.getDivergence(i, j);

        // 浮力源项（关键）
        glm::vec2 pos((i + 0.5f) * h, (j + 0.5f) * h);
        double buoyAccel = mGrid.getBoussinesqForce(pos) / rho;
        double buoySource = buoyAccel / h;

        divergence(i, j) = div - buoySource;
      }

      Glb::GridData2d p, pNew;
      p.initialize(0.0);
      pNew.initialize(0.0);

      const int maxIter = 40;
      const double scale = rho / dt;

      for (int it = 0; it < maxIter; ++it)
      {
        PIC_FOR_EACH_CELL
        {
          if (mGrid.isSolidCell(i, j))
          {
            pNew(i, j) = 0;
            continue;
          }

          double rhs = scale * divergence(i, j);
          double sum = 0.0;
          int cnt = 0;

          if (!mGrid.isSolidCell(i - 1, j))
          {
            sum += p(i - 1, j);
            cnt++;
          }
          if (!mGrid.isSolidCell(i + 1, j))
          {
            sum += p(i + 1, j);
            cnt++;
          }
          if (!mGrid.isSolidCell(i, j - 1))
          {
            sum += p(i, j - 1);
            cnt++;
          }
          if (!mGrid.isSolidCell(i, j + 1))
          {
            sum += p(i, j + 1);
            cnt++;
          }

          pNew(i, j) = (cnt > 0) ? (sum - rhs * h2) / cnt : 0.0;
        }
        p = pNew;
      }

      // Apply pressure gradient
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx + 1; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::X))
          {
            mGrid.mU(i, j) = 0;
            continue;
          }
          double pR = mGrid.isSolidCell(i, j) ? 0 : p(i, j);
          double pL = mGrid.isSolidCell(i - 1, j) ? 0 : p(i - 1, j);
          mGrid.mU(i, j) -= dt * (pR - pL) / (rho * h);
        }

      for (int j = 0; j < ny + 1; ++j)
        for (int i = 0; i < nx; ++i)
        {
          if (mGrid.isSolidFace(i, j, PICGrid2d::Y))
          {
            mGrid.mV(i, j) = 0;
            continue;
          }
          double pT = mGrid.isSolidCell(i, j) ? 0 : p(i, j);
          double pB = mGrid.isSolidCell(i, j - 1) ? 0 : p(i, j - 1);
          mGrid.mV(i, j) -= dt * (pT - pB) / (rho * h);
        }
    }

    /* ===================== G2P ===================== */

    void Solver::gridToParticle(double dt)
    {
      const float flipAlpha = 0.95f;

      for (auto &p : mPs.particles)
      {
        glm::vec2 v_pic = mGrid.getVelocity(p.position);
        glm::vec2 v_oldg = mGrid.getVelocityFromGrid(p.position, mU_prev, mV_prev);
        glm::vec2 v_flip = p.velocity + (v_pic - v_oldg);

        p.velocity = (1.0f - flipAlpha) * v_pic + flipAlpha * v_flip;
      }
    }

    /* ===================== Advect ===================== */

    void Solver::advectParticles(double dt)
    {
      const float h = mGrid.cellSize;
      const glm::vec2 lower(0.0f);
      const glm::vec2 upper((mGrid.dim[0] - 1) * h, (mGrid.dim[1] - 1) * h);

      for (auto &p : mPs.particles)
      {
        p.position += (float)dt * p.velocity;

        if (p.position.x < lower.x)
        {
          p.position.x = lower.x;
          p.velocity.x *= -0.5f;
        }
        if (p.position.x > upper.x)
        {
          p.position.x = upper.x;
          p.velocity.x *= -0.5f;
        }
        if (p.position.y < lower.y)
        {
          p.position.y = lower.y;
          p.velocity.y *= -0.5f;
        }
        if (p.position.y > upper.y)
        {
          p.position.y = upper.y;
          p.velocity.y *= -0.5f;
        }
      }
    }

  } // namespace PIC2d
} // namespace FluidSimulation

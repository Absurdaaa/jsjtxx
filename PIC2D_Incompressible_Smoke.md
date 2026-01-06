# 2D PIC/FLIP：不可压缩投影（MAC Grid）+ 更细腻烟雾渲染

> 对应你当前工程：`fluid2d/PIC`（CPU 侧 PIC/FLIP）+ `resources/shaders/DrawSmoke2d.*`（FBO + 全屏 Quad 贴图渲染烟雾）

本文包含两部分：

1. **在现有 PIC/FLIP 框架中加入/规范化 2D 不可压缩求解（MAC grid + pressure projection）**：给清晰的 C++ 代码结构/伪代码。
2. **提升 2D 烟雾渲染质量**：不改变 “FBO + 全屏 quad” 管线，只改密度生成/采样/着色方式，并给可直接用的 fragment shader。

同时我已经把两处关键改动直接落到了工程里：

- `fluid2d/PIC/src/Renderer.cpp`：密度贴图分辨率提高到 2x，并把粒子散布从双线性(2x2)升级为 **二次 B-spline 核(3x3)**，保留更大密度动态范围。
- `resources/shaders/DrawSmoke2d.frag`：改为 **bicubic 重建 + 指数映射/吸收 + 简单光照 + 自阴影** 的烟雾着色。

补充：如果你遇到“整屏发白/全白”或 shader 编译失败，本工程里已经做了两点关键修复：

- **二次 B-spline 权重必须是非负且归一**。错误的权重公式会让密度被过度累积（随后被钳制/映射），导致画面整体变白。
- GLSL 330 在一些驱动上不允许 `float s[4][4]` 这种“数组的数组”，bicubic 采样缓存已改成 **一维数组 `float s[16]`**。

---

## 1) 2D 不可压缩（MAC Grid + 散度 + 压力泊松 + 投影）

你当前的 2D PIC 网格已经是 MAC 形式（交错网格）：

- `PICGrid2d::mU` 类型 `Glb::GridData2dX`：$u$ 存在**竖直面中心**（维度 $(nx+1)\times ny$）
- `PICGrid2d::mV` 类型 `Glb::GridData2dY`：$v$ 存在**水平面中心**（维度 $nx\times (ny+1)$）
- `PICGrid2d::mD, mT`：密度/温度通常是 cell-center 标量（维度 $nx\times ny$）
- `PICGrid2d::mSolid`：固体标记（cell-center）

并且 `Solver::pressureProjection()` 已经做了一个“能跑起来”的投影版本。下面给出**更标准、更容易扩展/调参**的代码结构（你可以把现有实现拆成这些函数，逻辑会更清晰）。

### 1.1 建议的推进流程（和你现有 `Solver::solve()` 对齐）

你的 `Solver::solve()` 逻辑基本正确（emit → advect scalars → P2G → forces → projection → G2P → advect particles）。

标准不可压缩 PIC/FLIP 在网格侧的关键是：

- **投影在 P2G+外力之后**做
- **投影后的网格速度**再用于 PIC 采样、也用于 FLIP 的 $\Delta u$ 计算

伪代码（只列网格相关的关键步骤）：

```cpp
void Solver::solveOneStep(double dt)
{
    ps.emitFromSources();
    grid.updateSources();

    advectScalars(dt);           // 半拉格朗日：advect density/temp

    particleToGrid();            // P2G：把粒子速度散布到 MAC 网格面
    u_prev = grid.mU;            // FLIP 需要 old grid
    v_prev = grid.mV;

    addForces(dt);               // 重力/浮力/涡量约束等（在MAC面上加）

    applyVelocityBC();           // 固体/域边界速度边界条件（先清一次法向）
    projectIncompressible(dt);   // 计算div，解Poisson，减压强梯度
    applyVelocityBC();           // 投影后再做一次（防止数值漂移）

    gridToParticle(dt);          // G2P：PIC & FLIP 混合
    advectParticles(dt);         // 粒子用自身速度/或网格速度推进
}
```

### 1.2 MAC 网格上的散度（cell-center divergence）

对每个流体 cell $(i,j)$：

```cpp
// u(i,j) 是左面速度，u(i+1,j) 是右面速度
// v(i,j) 是下面速度，v(i,j+1) 是上面速度

div(i,j) = (u(i+1,j) - u(i,j) + v(i,j+1) - v(i,j)) / h;
```

与固体相邻时建议的处理方式（和你现在做法一致的精神）：

- **固体面上的法向速度视为 0**（no-penetration）
- 计算 divergence 时，如果某个面邻接固体，就把该面速度当作 0

代码结构建议：

```cpp
void computeDivergence(GridData2d& div, const GridData2dX& u, const GridData2dY& v)
{
    for each cell (i,j):
        if solid(i,j): continue;

        float uR = solid(i+1,j) ? 0 : u(i+1,j);
        float uL = solid(i-1,j) ? 0 : u(i,  j);
        float vT = solid(i, j+1) ? 0 : v(i, j+1);
        float vB = solid(i, j-1) ? 0 : v(i, j);

        div(i,j) = (uR - uL + vT - vB) / h;
}
```

> 你工程里对应：`PICGrid2d::getDivergence()`、`Solver::pressureProjection()`。

### 1.3 压力泊松方程（Poisson）与系数组装

经典 smoke（不可压缩空气）投影步骤：

$$\nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot u^*$$

离散到 cell-center（单位格宽 $h$）：

- 右端项：$b = (\rho/\Delta t)\, div$
- 系数矩阵是 5 点拉普拉斯；遇到固体邻居时，把固体视为 Neumann（不穿透）通常等价于：
  - 该方向不参与邻居项，中心系数减 1（也就是“有效邻居数”变少）

你目前用的是一种等价形式：

- `scale = dt / (rho * h * h)`
- Gauss-Seidel 更新：

```cpp
p(i,j) = (sumNeighborP - div(i,j)/scale) / numFluidNeighbors;
```

这是 OK 的（因为 `div/scale = div * rho*h*h/dt`）。

更“模块化”的写法建议把这三件事拆出来：

1) 计算 `div`
2) 解 `p`
3) 应用压力梯度更新 `u,v`

### 1.4 压力求解器：GS/Jacobi（易写）或 PCG（更稳更快）

你现阶段网格不大时，**Gauss-Seidel 迭代**足够。

建议的最小结构（GS / Jacobi 二选一）：

```cpp
struct PressureSolveParams
{
    int iterations = 60;      // GS/Jacobi 迭代次数
    float omega = 1.0f;       // SOR: 1.0=GS, 1~2可加速(需小心)
};

void solvePressure_GaussSeidel(
    GridData2d& p,
    const GridData2d& div,
    float scale,                 // dt/(rho*h*h)
    const GridData2d& solid,
    int nx, int ny,
    PressureSolveParams par)
{
    p.initialize(0.0);

    for (int it = 0; it < par.iterations; ++it)
    {
        for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            if (solid(i,j)) continue;

            float sum = 0;
            int cnt = 0;
            if (!solid(i-1,j)) { sum += p(i-1,j); cnt++; }
            if (!solid(i+1,j)) { sum += p(i+1,j); cnt++; }
            if (!solid(i, j-1)) { sum += p(i, j-1); cnt++; }
            if (!solid(i, j+1)) { sum += p(i, j+1); cnt++; }

            if (cnt == 0) continue;

            float newP = (sum - div(i,j)/scale) / (float)cnt;
            p(i,j) = mix(p(i,j), newP, par.omega); // SOR可选
        }
    }

    // 可选：去掉常数漂移（闭域pressure只有相对值）
    // subtractMean(p);
}
```

如果你后面网格更大、希望更“硬核”稳定：建议上 **PCG（预条件共轭梯度）**，但它需要更完整的线性系统封装（邻居系数/对角元/预条件）。本文不展开数学推导，但可以给接口：

```cpp
// A*p = b,  A 是离散拉普拉斯（考虑solid后每个cell的有效邻居数不同）
struct PoissonSystem2D
{
    int nx, ny;
    GridData2d diag;      // 对角元素 A(i,i)
    GridData2d offx;      // 与 (i+1,j) 的系数 (通常 -1 或 0)
    GridData2d offy;      // 与 (i,j+1) 的系数
    GridData2d solid;
};

void buildPoissonSystem(PoissonSystem2D& A, const GridData2d& solid);
void solvePCG(GridData2d& p, const PoissonSystem2D& A, const GridData2d& b);
```

### 1.5 投影（pressure projection）：减去压力梯度

对 MAC 面速度：

- 内部竖直面 $u(i,j)$（在 cell `(i-1,j)` 和 `(i,j)` 之间）：

$$u \leftarrow u - \frac{\Delta t}{\rho} \frac{p(i,j)-p(i-1,j)}{h}$$

- 内部水平面 $v(i,j)$（在 cell `(i,j-1)` 和 `(i,j)` 之间）：

$$v \leftarrow v - \frac{\Delta t}{\rho} \frac{p(i,j)-p(i,j-1)}{h}$$

对应到你当前代码里的 `scale = dt/(rho*h*h)` 写法：

```cpp
float gradScale = scale * h; // dt/(rho*h)

u(i,j) -= gradScale * (p(i,j) - p(i-1,j));
v(i,j) -= gradScale * (p(i,j) - p(i,j-1));
```

关键的“边界/固体”条件（建议保持你现在的逻辑）：

- **只有当面两侧都是流体 cell** 时才更新该面速度
- 若一侧是 solid，法向速度保持 0（不穿透）

结构化实现：

```cpp
void applyPressureGradient(
    GridData2dX& u, GridData2dY& v,
    const GridData2d& p,
    const GridData2d& solid,
    float dt, float rho, float h,
    int nx, int ny)
{
    float s = dt / (rho * h);

    // u: i=1..nx-1 是内部面
    for (int j=0;j<ny;++j)
    for (int i=1;i<nx;++i)
        if (!solid(i-1,j) && !solid(i,j))
            u(i,j) -= s * (p(i,j) - p(i-1,j));

    // v: j=1..ny-1 是内部面
    for (int j=1;j<ny;++j)
    for (int i=0;i<nx;++i)
        if (!solid(i,j-1) && !solid(i,j))
            v(i,j) -= s * (p(i,j) - p(i,j-1));
}
```

### 1.6 常见“看起来不不可压缩”的原因（排查清单）

如果你投影做了但仍然“体积变化明显”，通常是这些点：

- **P2G 没有对固体面做处理**：粒子散布把速度写进了 solid face，导致 divergence 计算/投影不一致
- **投影前/后没有强制速度边界条件**：建议 `applyVelocityBC()` 前后各一次
- **压力迭代次数不够**：迭代太少会残留 divergence → 体积会慢慢压缩/膨胀
- **dt 太大**：半拉格朗日会稳但数值扩散，投影迭代也更难收敛

你现在的 `Solver::pressureProjection()` 已经具备以上大部分要点；把它拆成模块并增加可控参数（迭代次数/omega）会更易调。

---

## 2) 提升 2D 烟雾渲染质量（仍然是 FBO + 全屏 Quad）

你当前渲染流程（CPU 粒子 → CPU 密度数组 → GL_R32F densityTexture → FBO → frag 根据密度上色）是很常见的。

“分辨率低/糊”的主要来源往往不是 FBO，而是：

1) **密度贴图分辨率太低**（直接等于模拟网格分辨率）
2) **粒子散布核太尖**（2x2 双线性会显得颗粒/走样）
3) shader 里对密度做线性上色，容易“灰一片”、缺少层次

### 2.1 我已做的 CPU 侧改进（不改管线，只改密度生成）

对应改动：`fluid2d/PIC/src/Renderer.cpp`

- 密度贴图分辨率改为：

```cpp
const int renderScale = 2;
gridResX = simNx * renderScale;
gridResY = simNy * renderScale;
```

- 粒子散布从 **双线性 2x2** 改为 **二次 B 样条 3x3**（更柔、更像烟），并把密度钳制从 `[0,1]` 放宽到 `[0,10]`，把“如何显示得好看”交给 shader 的指数映射。

你可以继续往上加：

- `renderScale = 3/4` 会更细腻，但 CPU 散布成本会线性上升
- 散布核可以升级为 4x4（cubic B-spline）更平滑

### 2.2 密度采样/插值：为什么要 bicubic

即使纹理用了 `GL_LINEAR`，低分辨率密度图上采样依然会显得“糊 + 缺细节”。

在 fragment shader 中做 **bicubic 重建（16 taps）** 可以：

- 在放大时减少格子感
- 边缘更软、更连续
- 更适合后续做梯度光照/自阴影

### 2.3 推荐的烟雾上色（指数吸收 + 软边缘 + 简单光照 + 自阴影）

对应改动：`resources/shaders/DrawSmoke2d.frag`

核心思路：

- 用 Beer-Lambert 指数吸收：`alpha = 1 - exp(-sigma * density)`（比线性 mix 更自然）
- 用密度梯度近似法线做 lambert（让烟“有体积”）
- 沿光方向累积密度做一个超轻量的“自阴影”
- 加极小 dither 降低 banding

我已经把一份可直接用的 shader 写进了 `DrawSmoke2d.frag`。你可以调的关键常量：

- `sigma`：烟厚度
- `steps / stepLen`：自阴影质量/开销
- `L`：光照方向

---

## 3) 文件改动一览

- 已修改：`fluid2d/PIC/src/Renderer.cpp`
  - 密度纹理分辨率提升（2x）
  - 二次 B 样条 3x3 粒子散布核
  - 密度动态范围更合理（交给 shader 指数映射）

- 已修改：`resources/shaders/DrawSmoke2d.frag`
  - bicubic 密度采样
  - 指数吸收
  - 简单光照 + 自阴影

---

## 4) 快速调参建议（从“像烟”到“更像烟”）

1) 如果烟“太淡/太透明”：增大 `sigma`（例如 2.8 → 4.0）
2) 如果烟“太硬/太黑”：减小 `sigma`，或把 `shadow` 的系数（`0.18`）调小
3) 如果你提高了 `renderScale`，但 CPU 变慢：
   - 把 `steps` 从 10 降到 6
   - 或把自阴影整个关掉（把 `shadow` 固定为 1.0）

---

## 5) 你下一步可以怎么把不可压缩做“更标准”

你已经有了基本投影。要让结果更稳定、更像真实烟（不可压缩 + 细腻），建议按优先级做：

1) 把 `pressureProjection()` 拆成 `computeDivergence / solvePressure / applyPressureGradient / applyVelocityBC`（更可控）
2) 把浮力/外力重新打开（你目前注释掉了浮力那段），并确保外力施加在 MAC 面上
3) 迭代次数做成可调（例如 UI slider：20~120）

如果你希望我直接把第 1) 的拆分也落到代码里（新增 `Projector2D` 类、把 `Solver::pressureProjection()` 改得更干净），告诉我你更倾向：

- 继续用 GS/SOR（简单、够用）
- 还是上 PCG（更专业、可扩展）

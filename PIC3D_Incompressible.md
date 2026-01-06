# 3D PIC/FLIP：不可压缩投影（MAC Grid 压力泊松 + 投影）

> 对应工程目录：`fluid3d/PIC`（3D PIC/FLIP），网格基类为 `Eulerian3d::MACGrid3d`。
>
> 本文只记录 **3D PIC** 相关的改动，2D 版本请看：`PIC2D_Incompressible_Smoke.md`。

---

## 1) 目标与结果

目标：让 3D PIC/FLIP 的网格速度场在每步外力/转移后满足不可压缩约束（$\nabla\cdot\mathbf{u}=0$）。

已完成：在 `fluid3d/PIC/src/Solver3d.cpp` 中实现/规范化 `Solver3d::pressureProjection()`，包含：

1. 固体面法向速度清零（no-penetration）
2. 基于 **solid face** 的散度计算
3. Gauss–Seidel 迭代解压力泊松
4. 对面速度减去压力梯度（仅两侧为流体 cell 的面）

---

## 2) 3D MAC 网格变量（你工程里的对应关系）

`Eulerian3d::MACGrid3d` 采用交错网格（MAC）：

- `mU`：$(nx+1)\times ny\times nz$（X 方向速度，YZ 面中心）
- `mV`：$nx\times(ny+1)\times nz$（Y 方向速度，XZ 面中心）
- `mW`：$nx\times ny\times(nz+1)$（Z 方向速度，XY 面中心）
- `mD/mT`：cell-center 标量场
- `isSolidCell(i,j,k)`：cell 是否为固体
- `isSolidFace(i,j,k,dir)`：该 face 是否为固体面（容器边界/固体邻接）

---

## 3) 压力投影实现（`Solver3d::pressureProjection()`）

对应文件：`fluid3d/PIC/src/Solver3d.cpp`

### 3.1 固体面速度清零

在投影开始，先把所有 `isSolidFace` 的面速度置零：

- `dir = X`：清 `mU(i,j,k)`
- `dir = Y`：清 `mV(i,j,k)`
- `dir = Z`：清 `mW(i,j,k)`

这是为了让后续散度计算与边界条件一致（固体面法向速度视为 0）。

### 3.2 散度（cell-center divergence）

对每个流体 cell $(i,j,k)$：

$$
\nabla\cdot \mathbf{u} = \frac{u_R-u_L + v_T-v_B + w_F-w_K}{h}
$$

并且当某个面是 solid face 时，该面速度按 0 处理（避免“solid cell 邻居/solid face”判定不一致导致的漏气）。

### 3.3 Gauss–Seidel 解泊松

离散形式等价于：

$$
\nabla^2 p = \frac{\rho}{\Delta t} \nabla\cdot\mathbf{u}^*
$$

代码里使用：

- `scale = dt / (rho * h * h)`
- GS 更新：

$$
 p = \frac{\sum p_{nb} - div/scale}{n}
$$

其中 $n$ 为有效流体邻居数（邻居若为固体 cell 则不计入）。

### 3.4 投影：减去压力梯度

只在“面两侧都是流体 cell”的内部面上更新：

- `mU(i,j,k) -= (dt/(rho*h)) * (p(i,j,k)-p(i-1,j,k))`
- `mV(i,j,k) -= (dt/(rho*h)) * (p(i,j,k)-p(i,j-1,k))`
- `mW(i,j,k) -= (dt/(rho*h)) * (p(i,j,k)-p(i,j,k-1))`

---

## 4) 新增可调参数：`pressureIters`

对应文件：

- `common/include/Configure.h`
- `common/src/Configure.cpp`

新增：`PIC3dPara::pressureIters`（默认 `50`）。

含义：GS 迭代次数。数值越大越接近无散度，但每步更慢。

---

## 5) UI 暴露（Inspector → PIC 3d）

对应文件：`ui/src/InspectorView.cpp`（`case 6`）

- `Delta Time` 绑定为 `PIC3dPara::dt`（保证调参能影响 3D PIC）
- 新增 `Pressure Iters` 绑定 `PIC3dPara::pressureIters`

---

## 6) 相关默认值调整（性能/稳定性）

对应文件：`common/src/Configure.cpp`

- `PIC3dPara::particlesPerStep` 默认值已设为 `100`（避免 3D 默认粒子数过大导致性能问题）。

---

## 7) 变更清单（3D 部分）

- `fluid3d/PIC/src/Solver3d.cpp`
  - 规范化 `pressureProjection()`：solid face 清零、散度、GS、压力梯度更新
- `common/include/Configure.h`
  - 新增 `PIC3dPara::pressureIters`
- `common/src/Configure.cpp`
  - 默认 `PIC3dPara::pressureIters = 50`
  - 默认 `PIC3dPara::particlesPerStep = 100`
- `ui/src/InspectorView.cpp`
  - PIC 3d 面板新增 `Pressure Iters`
  - PIC 3d 面板 `Delta Time` 绑定修正为 `PIC3dPara::dt`

---

## 8) 场景：中心球 + X+ 出流 + 面上圆形喷口

本次把 3D PIC 场景整理为与 2D 类似的结构：

- **中心球体障碍**：用 SDF（signed distance field）生成 solid cell；并用 *face-center* SDF 判定 solid face，让球边界更圆。
- **X+ 出流（outflow）**：右侧 $+X$ 边界不当作固体面（不做 no-penetration），允许物质/粒子从右侧离开。
- **面上圆形喷口**：从“体积发射”改为“在固定 $x$ 层的一张薄片上、$y$-$z$ 平面内的圆盘区域采样”。

对应文件：

- `fluid3d/PIC/include/PICGrid3d.h`
- `fluid3d/PIC/src/PICGrid3d.cpp`
- `fluid3d/PIC/src/Solver3d.cpp`
- `fluid3d/PIC/src/ParticleSystem3d.cpp`

### 8.1 球半径参数（cell 单位）

新增参数：`PIC3dPara::sphereRadiusCells`（单位：cell）。

- 默认值：`6.5f`
- 使用位置：`PICGrid3d::initialize()` 会把该值转换为世界单位半径并 clamp 到合理范围。
- 注意：半径在网格初始化时决定（固体体素/solid face 也依赖它），因此 **修改后通常需要 rerun/restart 才能完全生效**。

对应文件：

- `common/include/Configure.h`
- `common/src/Configure.cpp`
- `fluid3d/PIC/src/PICGrid3d.cpp`

### 8.2 X+ 出流

- `PICGrid3d::isSolidFace(...)`：对 $+X$ 边界返回 `false`（表示该边界面不是 solid face）。
- 粒子层面：advect 后越过右侧边界的粒子会被移除，避免在域外累积。

### 8.3 面上圆形喷口（YZ 圆盘）

喷口改为“面上圆盘”：

- 固定在 `source` 的 $x$ 位置附近（极薄厚度）。
- 在 $y$-$z$ 平面用拒绝采样（rejection sampling）保证 $dy^2 + dz^2 \le r^2$，减少 clamp 导致的偏置。
- 喷口读取 `Eulerian3dPara::source`（因为 UI 的 source 面板目前编辑的是 Eulerian3d 的 source；这样 UI 改源能直接影响 PIC3d 喷口）。

---

## 9) UI：球半径 + 相机参数（PIC 3D）

对应文件：`ui/src/InspectorView.cpp`（`case 6`）

- 新增 `Sphere Radius (cells)`：绑定 `PIC3dPara::sphereRadiusCells`。
- 新增 `Camera` 区块（仿照 Eulerian3d）：直接编辑 `Glb::Camera::Instance()` 的
  - `mPosition`
  - `mYaw / mPitch`
  - `fovyDeg / aspect / nearPlane / farPlane`
  - 并提供按钮调用 `UpdateView()`

---

## 10) 变更清单（补充）

- `common/include/Configure.h`
  - 新增 `PIC3dPara::sphereRadiusCells`
- `common/src/Configure.cpp`
  - 默认 `PIC3dPara::sphereRadiusCells = 6.5f`
- `fluid3d/PIC/include/PICGrid3d.h`
  - 增加中心球 SDF/solid face 支持（接口/辅助函数）
- `fluid3d/PIC/src/PICGrid3d.cpp`
  - center sphere solid cell 生成、face-center solid face 判定
  - X+ outflow（右侧边界不当作 solid face）
- `fluid3d/PIC/src/ParticleSystem3d.cpp`
  - 面上圆盘喷口采样（YZ plane）
- `ui/src/InspectorView.cpp`
  - PIC 3d 面板新增 Sphere Radius 与 Camera 参数

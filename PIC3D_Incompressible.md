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

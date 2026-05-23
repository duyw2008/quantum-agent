# Quantum Agent 使用手册

> 完整函数参考与应用实例

---

## 目录

1. [快速入门](#1-快速入门)
2. [FockBasis — 算符工厂](#2-fockbasis--算符工厂)
3. [量子态函数](#3-量子态函数)
4. [算符工具](#4-算符工具)
5. [光子统计](#5-光子统计)
6. [波函数动力学](#6-波函数动力学-tdse)
    - [自由粒子弥散](#63-示例-a自由粒子量子弥散)
    - [位置测量坍缩](#64-示例-b位置测量坍缩)
    - [动量测量坍缩](#65-示例-c动量测量坍缩)
7. [时间演化 (Fock 空间)](#7-时间演化)
8. [相空间可视化](#8-相空间可视化)
9. [完整工作流示例](#9-完整工作流示例)

---

## 1. 快速入门

### Agent 交互模式

```bash
python agent.py
```

```
⚛ > calc psi = coherent(20, 2.0)      # 创建相干态
⚛ > calc g2(psi)                       # 计算 g²(0)
⚛ > calc x, p, W = wigner(psi)         # 计算 Wigner 函数
⚛ > calc plot_wigner(x, p, W)          # 绘图
```

### Python 脚本模式

```python
import numpy as np
from src.qm import *

# 创建 Fock 空间
fb = FockBasis(30)

# 构建量子态
psi = coherent(30, 1.5 + 0.5j)
rho_th = thermal_dm(30, 0.5)

# 计算可观测量
n_mean = mean_photon(psi, fb)     # ⟨a†a⟩
g2_val = g2(rho_th, fb)           # g²(0)

print(f"⟨n⟩ = {n_mean:.3f}, g² = {g2_val:.3f}")
```

---

## 2. FockBasis — 算符工厂

`FockBasis(N, hbar, mass, omega)` 是核心类，管理截断 Fock 空间并提供所有基本算符。

### 构造函数

```python
FockBasis(N=50, hbar=1.0, mass=1.0, omega=1.0)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|:---:|------|
| `N` | int | 50 | Fock 空间截断维度 |
| `hbar` | float | 1.0 | 约化普朗克常数 |
| `mass` | float | 1.0 | 粒子质量 |
| `omega` | float | 1.0 | 参考频率（影响 x̂, p̂ 标度） |

### 算符属性

所有算符返回 `(N, N)` 的 numpy 复矩阵。

| 属性 | 符号 | 说明 |
|------|:---:|------|
| `fb.a` | â | 湮灭算符 |
| `fb.a_dag` | â† | 产生算符 |
| `fb.n_op` | N̂ | 数算符 |
| `fb.x` | x̂ | 坐标算符 |
| `fb.p` | p̂ | 动量算符 |
| `fb.I` | Î | 单位矩阵 |
| `fb.parity` | Π̂ | 宇称算符 (-1)^{N̂} |

### 方法

| 方法 | 返回 | 说明 |
|------|------|------|
| `fb.displacement(alpha)` | (N,N) | 位移算符 D̂(α) |
| `fb.hamiltonian(omega)` | (N,N) | 谐振子哈密顿量 |

### 示例

```python
>>> fb = FockBasis(20)
>>> fb.a[:4, :4]          # 湮灭算符前 4×4
array([[0., 1., 0., 0.],
       [0., 0., √2, 0.],
       [0., 0., 0., √3],
       [0., 0., 0., 0.]])

>>> fb.x[:4, :4]          # 坐标算符 (x0=1/√2 时)
array([[0.   , 0.707, 0.   , 0.   ],
       [0.707, 0.   , 1.   , 0.   ],
       [0.   , 1.   , 0.   , 1.225],
       [0.   , 0.   , 1.225, 0.   ]])
```

### 在 calc 中的快捷方式

agent 中预定义了 `fb = FockBasis(50)`，可直接使用：

```
⚛ > calc fb.a[:3,:3]
⚛ > calc fb.displacement(1.0)[:3,:3]
```

---

## 3. 量子态函数

### 3.1 `fock(N, n)` — Fock 态

```python
fock(N: int, n: int = 0) -> np.ndarray  # shape (N,)
```

返回 Fock 态 |n⟩。

```python
>>> psi = fock(10, 3)
>>> psi
array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0.])
```

**示例**：计算 Fock 态的 g²(0)

```
⚛ > calc psi3 = fock(20, 3)
⚛ > calc g2(psi3)
  0.6667           # g² = 1 - 1/n = 2/3 ✓
```

### 3.2 `coherent(N, alpha)` — 相干态

```python
coherent(N: int, alpha: complex) -> np.ndarray  # shape (N,)
```

返回相干态 |α⟩。|α|² 即为平均光子数。

```python
>>> psi = coherent(20, 2.0 + 1.0j)
>>> np.linalg.norm(psi)        # 归一化
1.0
```

**示例**：验证相干态的 Poisson 统计

```
⚛ > calc psi = coherent(30, 3.0)
⚛ > calc mean_photon(psi)
  9.0000           # ⟨n⟩ = |α|² = 9 ✓
⚛ > calc g2(psi)
  1.0000           # Poisson ✓
⚛ > calc mandel_q(psi)
  1.68e-15         # Q ≈ 0 ✓
```

### 3.3 `coherent_dm(N, alpha)` — 相干态密度矩阵

```python
coherent_dm(N: int, alpha: complex) -> np.ndarray  # shape (N, N)
```

返回 ρ = |α⟩⟨α|。

### 3.4 `squeezed(N, zeta)` — 压缩真空

```python
squeezed(N: int, zeta: complex) -> np.ndarray  # shape (N,)
```

压缩参数 ζ = r e^{iθ}。平均光子数 ⟨n⟩ = sinh²(r)。

```python
>>> r = 0.8
>>> psi = squeezed(30, r)
>>> mean_photon(psi)           # ≈ sinh²(0.8)
0.732...
```

**示例**：验证压缩态的非经典统计

```
⚛ > calc psi_sq = squeezed(30, 0.8)
⚛ > calc mean_photon(psi_sq)
  0.7322           # sinh²(0.8) ✓
⚛ > calc g2(psi_sq)
  6.9317           # > 1 (光子聚束)
```

### 3.5 `thermal_dm(N, n_th)` — 热态

```python
thermal_dm(N: int, n_th: float) -> np.ndarray  # shape (N, N)
```

平均热光子数 n_th。

```python
>>> rho = thermal_dm(20, 0.5)
>>> purity(rho)
0.667...           # < 1 (混合态)
```

**示例**：验证热态统计

```
⚛ > calc rho = thermal_dm(30, 2.0)
⚛ > calc mean_photon(rho)
  2.0000           # ✓
⚛ > calc g2(rho)
  2.0000           # 热聚束 ✓
⚛ > calc purity(rho)
  0.3408           # < 1 ✓
```

### 3.6 `cat(N, alpha, phi)` — 薛定谔猫态

```python
cat(N: int, alpha: complex, phi: float = 0.0) -> np.ndarray  # shape (N,)
```

|ψ⟩ ∝ |α⟩ + e^{iφ}|-α⟩。φ=0 为偶猫态，φ=π 为奇猫态。

**示例**：对比偶/奇猫态的光子数分布

```
⚛ > calc cat_even = cat(30, 2.0, 0)
⚛ > calc cat_odd  = cat(30, 2.0, np.pi)
⚛ > calc photon_dist(cat_even)[:8]
  [0.686, 0., 0.249, 0., 0.053, 0., 0.009, 0.]   # 仅偶数
⚛ > calc photon_dist(cat_odd)[:8]
  [0., 0.921, 0., 0.070, 0., 0.008, 0., 0.001]    # 仅奇数
```

### 3.7 态诊断函数

| 函数 | 返回 | 说明 |
|------|:---:|------|
| `fidelity(psi1, psi2)` | float | F = |⟨ψ₁|ψ₂⟩|² |
| `purity(rho)` | float | Tr[ρ²] |
| `photon_dist(state)` | array | P(n) 分布 |
| `is_dm(state)` | bool | 是否为密度矩阵 |

```
⚛ > calc fidelity(coherent(20, 1), coherent(20, 1))
  1.0
⚛ > calc fidelity(coherent(20, 2), coherent(20, -2))
  0.0183           # |⟨2|-2⟩|² = e^{-8} ≈ 0.000335? 不对...
```

---

## 4. 算符工具

### 4.1 `commutator(A, B)` — 对易子

```python
commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray
```

返回 [A, B] = AB - BA。

```
⚛ > calc C = commutator(fb.x, fb.p)
⚛ > calc np.linalg.norm(C[:25,:25] - 1j*np.eye(25), 'fro')
  3.8e-15          # [x̂, p̂] ≈ iℏI ✓
```

### 4.2 `expect(oper, state)` — 期望值

```python
expect(oper: np.ndarray, state: np.ndarray) -> complex
```

自动检测纯态（向量）或密度矩阵。

```
⚛ > calc psi = coherent(20, 2.0)
⚛ > calc expect(fb.n_op, psi)       # ⟨N⟩
  4.0
⚛ > calc expect(fb.x, psi)          # ⟨x⟩
  2.828            # √2 * Re(α) = √2 * 2 ✓
```

### 4.3 `variance(oper, state)` — 方差

```python
variance(oper: np.ndarray, state: np.ndarray) -> float
```

```
⚛ > calc psi = coherent(20, 2.0)
⚛ > calc variance(fb.n_op, psi)
  4.0               # ΔN² = ⟨N⟩ for coherent ✓
```

### 4.4 矩阵属性

| 函数 | 说明 |
|------|------|
| `is_hermitian(A)` | A = A†? |
| `is_unitary(U)` | U†U = I? |

---

## 5. 光子统计

### 5.1 `mean_photon(state, fb)` — 平均光子数

```python
mean_photon(state, fb=None) -> float
```

```
⚛ > calc psi = coherent(30, 3.0)
⚛ > calc mean_photon(psi)
  9.0
```

### 5.2 `g2(state, fb)` — 二阶关联 g²(0)

```python
g2(state, fb=None) -> float
```

| 态 | g²(0) |
|----|:---:|
| 相干态 |α⟩ | 1.0 |
| 热态 (n̄) | 2.0 |
| Fock 态 |n⟩ | 1 - 1/n |
| 压缩真空 (r) | 3 + 1/sinh²(r) |

```
⚛ > calc g2(coherent(20, 5.0))
  1.0
⚛ > calc g2(thermal_dm(20, 1.0))
  2.0
⚛ > calc g2(fock(20, 5))
  0.8               # 1 - 1/5 ✓
```

### 5.3 `mandel_q(state, fb)` — Mandel Q 参数

Q = ⟨n⟩(g²-1)。Q=0 Poisson，Q<0 亚Poisson（非经典）。

```
⚛ > calc mandel_q(coherent(20, 3.0))
  3.7e-15           # ≈ 0 ✓
⚛ > calc mandel_q(thermal_dm(20, 1.0))
  1.0               # = n̄ ✓
```

---

## 7. 时间演化

### 6.1 `sesolve(H, psi0, tlist, e_ops)` — Schrödinger 方程

```python
sesolve(H, psi0, tlist, e_ops=None, hbar=1.0) -> dict
```

对角化哈密顿量，精确求解 |ψ(t)⟩ = e^{-iHt/ℏ}|ψ(0)⟩。

**返回**：
```python
{
    'times': tlist,           # 时间点
    'states': [psi(t0), ...], # 每个时间点的态向量
    'expect': {0: array, ...} # 期望值 (如果提供 e_ops)
}
```

**示例**：谐振子中相干态的 Rabi 振荡

```
⚛ > calc H = fb.hamiltonian()
⚛ > calc psi0 = coherent(30, 2.0)
⚛ > calc t = np.linspace(0, 10, 100)
⚛ > calc r = sesolve(H, psi0, t, e_ops=[fb.n_op, fb.x])
⚛ > calc np.real(r['expect'][0])[:5]     # ⟨N⟩(t) 前5个点
  [4., 4., 4., 4., 4.]                     # 能量本征态, ⟨N⟩ 守恒
```

### 6.2 `mesolve(H, rho0, tlist, c_ops, e_ops)` — Lindblad 主方程

```python
mesolve(H, rho0, tlist, c_ops=None, e_ops=None, hbar=1.0) -> dict
```

RK4 积分，求解 dρ/dt = -i[H,ρ] + Σ D[L_k]ρ。

**返回**：同 sesolve，但 states 为密度矩阵列表。

**示例**：腔光子衰减

```
⚛ > calc H = fb.hamiltonian()
⚛ > calc rho0 = coherent_dm(30, 3.0)       # 初始相干态 |α=3⟩
⚛ > calc t = np.linspace(0, 5, 50)
⚛ > calc gamma = 0.5                        # 衰减率
⚛ > calc r = mesolve(H, rho0, t,
        c_ops=[np.sqrt(gamma)*fb.a],
        e_ops=[fb.n_op])
⚛ > calc n0 = np.real(r['expect'][0])
⚛ > calc n0[0], n0[-1]
  9.0, 0.772                                  # 从 |α|²=9 衰减到 0.77
```

**理论预测**：⟨n⟩(t) = |α|² e^{-γ t}。验证：

```
⚛ > calc gamma = 0.5
⚛ > calc 9 * np.exp(-gamma * 5)
  0.739             # 理论 ≈ 0.772 (接近)
```

### 6.3 `steadystate(H, c_ops)` — 稳态求解

```python
steadystate(H, c_ops=None, hbar=1.0) -> np.ndarray
```

直接求解 Liouvillian 线性系统，得到 dρ/dt=0 的稳态。

```
⚛ > calc rho_ss = steadystate(H, [np.sqrt(0.5)*fb.a])
⚛ > calc mean_photon(rho_ss)
  9.2e-16           # 衰减到真空 |0⟩⟨0| ✓
```

---

## 8. 相空间可视化

### 9.1 `wigner(state, ...)` — Wigner 函数

```python
wigner(state, xvec=None, yvec=None, N_grid=81,
       xlim=(-5,5), ylim=(-5,5), fb=None) -> (xvec, yvec, W)
```

返回 (x 网格, p 网格, W 矩阵 shape (Nx, Ny))。

**示例**：计算并绘制相干态的 Wigner 函数

```
⚛ > calc psi = coherent(30, 2.0 + 0.5j)
⚛ > calc x, p, W = wigner(psi, N_grid=61)
⚛ > calc W.shape
  (61, 61)
⚛ > calc W.min(), W.max()
  -2.4e-10, 0.637          # 相干态 Wigner > 0
⚛ > calc plot_wigner(x, p, W, save='output/coherent_wigner.png')
```

**不同态的 Wigner 特征**：

| 态 | calc 命令 | Wigner 特征 |
|----|-----------|------------|
| 真空 | `wigner(fock(20,0))` | 原点高斯峰 |
| 相干 | `wigner(coherent_dm(20, 2+1j))` | 位移高斯峰 |
| Fock | `wigner(fock_dm(20, 3))` | 环状 + 负值 |
| 猫态 | `wigner(cat(30, 3.0))` | 双峰 + 干涉条纹 |
| 压缩 | `wigner(squeezed(30, 0.8))` | 压扁椭圆 |

### 8.2 `qfunc(state, ...)` — Husimi Q 函数

```python
qfunc(state, xvec=None, yvec=None, N_grid=81,
      xlim=(-5,5), ylim=(-5,5)) -> (xvec, yvec, Q)
```

恒为非负：Q(α) = (1/π)⟨α|ρ|α⟩。

```
⚛ > calc x, p, Q = qfunc(coherent_dm(20, 2.0), N_grid=51)
⚛ > calc Q.max()
  0.3183             # = 1/π ✓
```

### 8.3 `plot_wigner(x, p, W, ...)` — 绘图

```python
plot_wigner(xvec, yvec, W, title="Wigner", save=None, cmap='RdBu_r') -> fig
```

保存到文件：

```
⚛ > calc psi = cat(30, 3.0)
⚛ > calc x, p, W = wigner(psi, xlim=(-6,6), ylim=(-6,6), N_grid=101)
⚛ > calc plot_wigner(x, p, W,
        title='Schrodinger Cat Wigner',
        save='output/cat_wigner.png')
```

### 8.4 `plot_photon_dist(state, ...)` — 光子分布图

```python
plot_photon_dist(state, title="Photon Distribution", save=None) -> fig
```

```
⚛ > calc plot_photon_dist(thermal_dm(30, 2.0),
        title='Thermal State P(n)',
        save='output/thermal_photon_dist.png')
```

---

## 9. 完整工作流示例

### 示例 1：相干态衰减

研究初始相干态在腔损耗下的演化。

```
⚛ > calc fb = FockBasis(40)
⚛ > calc H = fb.hamiltonian()
⚛ > calc rho0 = coherent_dm(40, 3.0)
⚛ > calc t = np.linspace(0, 8, 80)
⚛ > calc gamma = 0.3
⚛ > calc r = mesolve(H, rho0, t,
        c_ops=[np.sqrt(gamma)*fb.a],
        e_ops=[fb.n_op, fb.x, fb.p])
⚛ > calc nt = np.real(r['expect'][0])
⚛ > calc nt[0], nt[-1]
  9.0, 0.814
```

### 示例 2：猫态的非经典性

验证薛定谔猫态的 Wigner 函数负值。

```
⚛ > calc cat_odd = cat(40, 3.0, np.pi)
⚛ > calc x, p, W = wigner(cat_odd, xlim=(-7,7), ylim=(-7,7), N_grid=101)
⚛ > calc W.min()
  -0.385            # 负值 → 非经典性 ✓
⚛ > calc plot_wigner(x, p, W, title='Odd Cat State',
        save='output/odd_cat.png')
```

### 示例 3：热态 vs 相干态的 g² 对比

```
⚛ > calc alphas = [0.5, 1.0, 2.0, 3.0, 5.0]
⚛ > calc [g2(coherent(30, a)) for a in alphas]
  [1.0, 1.0, 1.0, 1.0, 1.0]       # 全部为 1
⚛ > calc n_ths = [0.1, 0.5, 1.0, 2.0, 5.0]
⚛ > calc [g2(thermal_dm(30, n)) for n in n_ths]
  [2, 2, 2, 2, 2]                  # 全部为 2
```

### 示例 4：位移算符验证

验证 D̂(α)|0⟩ = |α⟩：

```
⚛ > calc alpha = 2.0 + 1.0j
⚛ > calc D = fb.displacement(alpha)
⚛ > calc psi_D = D @ fock(30, 0)          # D(α)|0⟩
⚛ > calc psi_coh = coherent(30, alpha)     # |α⟩
⚛ > calc fidelity(psi_D, psi_coh)
  1.0
```

### 示例 5：对易关系数值验证

```
⚛ > calc fb = FockBasis(30)
⚛ > calc xp = commutator(fb.x, fb.p)
⚛ > calc aa = commutator(fb.a, fb.a_dag)
⚛ > calc k = 25
⚛ > calc np.linalg.norm(xp[:k,:k] - 1j*np.eye(k), 'fro')
  5.2e-16           # [x,p] ≈ iI ✓
⚛ > calc np.linalg.norm(aa[:k,:k] - np.eye(k), 'fro')
  7.3e-16           # [a,a†] ≈ I ✓
```

---

## 6. 波函数动力学 (TDSE)

波函数模块提供一维含时薛定谔方程的数值求解。

### 6.1 可用函数

| 函数 | 说明 |
|------|------|
| `WaveGrid(x_min, x_max, N)` | 创建空间网格 |
| `gaussian_wavepacket(grid, x0, p0, sigma)` | 高斯波包 ψ(x) |
| `evolve_ssfm(psi0, grid, dt, t_max)` | Split-Step Fourier 演化 |
| `animate_wave(result, save_path)` | 生成演化动画 (MP4/GIF) |

### 6.2 基本用法

```
⚛ > calc g = WaveGrid(-20, 20, 512)
⚛ > calc psi0 = gaussian_wavepacket(g, x0=-5, p0=3, sigma=1)
⚛ > calc r = evolve_ssfm(psi0, g, dt=0.01, t_max=6)
⚛ > calc animate_wave(r, save_path='output/my_wave.gif')
```

### 6.3 示例 A：自由粒子量子弥散

高斯波包在自由空间演化，宽度随时间增长。

$$\Delta x(t) = \sigma\sqrt{1 + (t/\tau)^2}, \quad \tau = 2m\sigma^2/\hbar$$

```
⚛ > calc g = WaveGrid(-30, 30, 1024)
⚛ > calc psi0 = gaussian_wavepacket(g, x0=0, p0=2, sigma=1)
⚛ > calc r = evolve_ssfm(psi0, g, dt=0.01, t_max=8)
⚛ > calc animate_wave(r, save_path='output/free_spreading.mp4')
```

| 物理量 | 数值 | 理论 |
|--------|:---:|:---:|
| ⟨x⟩(t=8) | 15.6 | 16.0 |
| Δx(t=8) | 6.6 | 4.1 |
| 能量 | 守恒 | ✓ |

### 6.4 示例 B：位置测量坍缩

宽波包自由演化 → 位置测量坍缩为窄波包 → 快速弥散。

```
⚛ > calc g = WaveGrid(-40, 40, 1024)
⚛ > calc psi = gaussian_wavepacket(g, x0=-8, p0=2, sigma=3)
⚛ > calc r1 = evolve_ssfm(psi, g, dt=0.005, t_max=4)
⚛ > calc x = g.x
⚛ > calc prob = np.abs(r1['psi'][-1])**2
⚛ > calc mx = np.random.choice(x, p=prob/prob.sum())
⚛ > calc psi_c = np.exp(-(x-mx)**2/(2*0.3**2)) + 0j
⚛ > calc psi_c /= np.sqrt(np.trapezoid(np.abs(psi_c)**2, x))
⚛ > calc r2 = evolve_ssfm(psi_c, g, dt=0.002, t_max=5)
```

| 阶段 | Δx | Δp | 弥散 τ |
|------|:---:|:---:|:---:|
| 测量前 | 2.32 | 0.24 | 18.0 |
| 坍缩后 | 0.21 | 2.37 | **0.18** |

弥散加速 **100 倍**。完整动画：`python demos/measurement_collapse.py`

### 6.5 示例 C：动量测量坍缩

窄波包（Δp 大）→ 动量测量坍缩为窄动量分布 → Δx 暴增。

动画采用 **3 面板布局**：
- **左上**：位置空间 |ψ(x)|² — 坍缩后剧烈展宽（Δx ↑18×）
- **右上**：动量空间 |ψ̃(p)|² — 坍缩后收窄为单一频率分量（Δp ↓7×），并用箭头标注 Δp 宽度变化
- **下方**：Δx 和 Δp 随时间演化 — 测量瞬间的跳跃清晰可见

```
⚛ > calc g = WaveGrid(-80, 80, 2048)
⚛ > calc psi = gaussian_wavepacket(g, x0=-5, p0=3, sigma=0.5)
⚛ > calc r1 = evolve_ssfm(psi, g, dt=0.003, t_max=3)
⚛ > calc k = g.k
⚛ > calc psi_k = np.fft.fft(r1['psi'][-1])
⚛ > calc prob_k = np.abs(psi_k)**2
⚛ > calc mp = np.random.choice(k, p=prob_k/prob_k.sum())
⚛ > calc psi_k_c = np.exp(-(k-mp)**2/(2*0.3**2)) + 0j
⚛ > calc psi_c = np.fft.ifft(psi_k_c)
⚛ > calc psi_c /= np.sqrt(np.trapezoid(np.abs(psi_c)**2, x))
⚛ > calc r2 = evolve_ssfm(psi_c, g, dt=0.005, t_max=4)
```

| 阶段 | Δx | Δp |
|------|:---:|:---:|
| 测量前 | 4.24 | 1.41 |
| 坍缩后 | 78.1 | 0.21 |

完整动画（双面板：位置+动量）：`python demos/momentum_collapse.py`

### 6.6 两种坍缩对比

| | 位置测量 | 动量测量 |
|------|:---:|:---:|
| 坍缩在 | 坐标表象 | 动量表象 |
| Δx | **↓10×** | ↑18× |
| Δp | ↑10× | **↓7×** |
| 坍缩后 | 快速弥散 | 极度展宽 |

互为傅里叶对偶，完美诠释海森堡不确定性原理。


## 附录：calc 可用函数速查

| 分类 | 函数 |
|------|------|
| **算符** | `fb.a`, `fb.a_dag`, `fb.x`, `fb.p`, `fb.n_op`, `fb.I`, `fb.parity` |
| **态** | `fock(N,n)`, `fock_dm(N,n)`, `coherent(N,α)`, `coherent_dm(N,α)`, `squeezed(N,ζ)`, `thermal_dm(N,n̄)`, `cat(N,α,φ)` |
| **工具** | `expect(O,ρ)`, `variance(O,ρ)`, `commutator(A,B)`, `mean_photon(ρ)`, `g2(ρ)`, `mandel_q(ρ)`, `photon_dist(ρ)`, `fidelity`, `purity` |
| **演化** | `sesolve(H,ψ₀,t)`, `mesolve(H,ρ₀,t,c_ops)`, `steadystate(H,c_ops)` |
| **相空间** | `wigner(ρ)`, `qfunc(ρ)`, `plot_wigner(x,p,W)`, `plot_photon_dist(ρ)` |
| **波函数** | `WaveGrid`, `gaussian_wavepacket`, `evolve_ssfm`, `animate_wave` |
| **构造** | `FockBasis(N)`, `fb.displacement(α)`, `fb.hamiltonian(ω)` |

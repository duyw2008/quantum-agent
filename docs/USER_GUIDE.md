# Quantum Agent 使用手册

> 完整函数参考与应用实例

---

## 目录

1. [两种使用模式](#1-两种使用模式)
   - [Python 脚本模式 — 独立 .py 调用函数库](#11-python-脚本模式--独立-py-调用函数库)
   - [Agent 脚本模式 — agent 内 .qms 批量执行](#12-agent-脚本模式--agent-内-qms-批量执行)
2. [FockBasis — 算符工厂](#2-fockbasis--算符工厂)
3. [量子态函数](#3-量子态函数)
4. [算符工具](#4-算符工具)
5. [光子统计](#5-光子统计)
6. [波函数动力学](#6-波函数动力学-tdse)
7. [时间演化 (Fock 空间)](#7-时间演化)
8. [相空间可视化](#8-相空间可视化)
9. [完整工作流示例](#9-完整工作流示例)

---

## 1. 两种使用模式

Quantum Agent 提供两套独立的使用方式，适合不同场景：

| 模式 | 适用场景 | 灵活性 | 复杂度 |
|------|----------|:---:|:---:|
| **Python 脚本模式** | 科研计算、数据处理、集成到工具链 | ⭐⭐⭐ | 需手动导入 |
| **Agent 脚本模式** | 教学演示、快速验证、交互探索 | ⭐⭐ | 函数预加载 |

### 1.1 Python 脚本模式 — 独立 .py 调用函数库

在普通 `.py` 文件中 `import` 函数库，完全控制执行流程。适合批量计算、数据处理、集成到 Jupyter。

```bash
python my_calculation.py
```

**完整示例：**

```python
import numpy as np
from src.qm import FockBasis, coherent, expect, got, mandel_q

# 创建 Fock 空间
fb = FockBasis(50)

# 构建量子态
psi = coherent(30, 2.0 + 0.5j)
rho_th = thermal_dm(30, 0.5)

# 计算可观测量
n_mean = mean_photon(psi, fb)      # ⟨a†a⟩
g2_val = g2(rho_th, fb)            # g²(0)
xp_error = np.linalg.norm(
    commutator(fb.x, fb.p)[:10,:10] - 1j*np.eye(10), 'fro'
)

print(f"⟨n⟩ = {n_mean:.3f}")
print(f"g²(0) = {g2_val:.3f}")
print(f"[x̂,p̂] error = {xp_error:.2e}")
```

**导入方式：**

```python
# 量子力学核心
from src.qm import FockBasis, fock, coherent, squeezed, thermal_dm, cat
from src.qm import expect, variance, g2, mandel_q, mean_photon
from src.qm import commutator, sesolve, mesolve, steadystate
from src.qm import fidelity, purity, photon_dist

# 波函数动力学
from src.qm import WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave

# 可视化
from src.viz import wigner, qfunc, plot_wigner, plot_photon_dist

# 量子场论 (可选)
from src.qft import ScalarField, LatticePhi4
from src.qft import wick_expand, feynman_amplitude_phi4_2to2
```

**Python 模式的优点：**

- 完全控制 `import`、`for`、`if`、`def` 等 Python 语法
- 可用 `print()`、`plt.show()`、`open()` 等所有标准库
- IDE 支持、断点调试、类型提示
- 适合生成论文图表、跑参数扫描、做 Monte Carlo

---

### 1.2 Agent 脚本模式 — agent 内 .qms 批量执行

在 agent 交互环境中编写 `.qms` 脚本，函数预加载免 import，适合快速验证和教学演示。

**启动 agent：**

```bash
python agent.py                  # 交互模式
python agent.py --run script.qms  # 直接执行脚本
```

**交互模式示例：**

```
⚛ ~/quantum_agent > psi = coherent(20, 2.0)
⚛ ~/quantum_agent > g2(psi)
1.0
⚛ ~/quantum_agent > x, p, W = wigner(psi)
⚛ ~/quantum_agent > plot_wigner(x, p, W)
```

**预加载函数（无需 import）：**

```
FockBasis, fock, coherent, squeezed, thermal_dm, cat,
expect, variance, g2, mandel_q, mean_photon,
commutator, sesolve, mesolve, steadystate,
wigner, qfunc, plot_wigner, plot_photon_dist,
WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave,
fb (默认 FockBasis(50)), np (numpy)
```

#### .qms 脚本格式

将命令序列写入 `.qms` 文件，批量执行：

```bash
# 命令行执行
python agent.py --run scripts/harmonic_oscillator.qms

# 交互模式内执行（支持 Tab 文件路径补全）
⚛ > run scripts/harmonic_oscillator.qms
```

**.qms 脚本写法：**

```
# 用 # 写注释
# agent 命令: formula, animate, run 等
# Python 表达式: 直接写

formula [\hat{x}, \hat{p}] = i\hbar          ← formula 命令

alpha = 2.0 + 0.5j                            ← Python 赋值
psi = coherent(30, alpha)
mean_photon(psi)                              ← Python 求值

C = commutator(fb.x, fb.p)                    ← fb 已预加载
g2(psi)

# 多行语句: dict / for / def 等
psi_sq = squeezed(30, 0.8)
x, p, W = wigner(psi_sq)

# matplotlib 作图
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(x, W[30, :])
plt.savefig('output/wigner_slice.png')
plt.close()
```

**.qms 脚本特性：**

| 特性 | 说明 |
|------|------|
| `# 注释` | 行注释，跳过空行 |
| 函数预加载 | FockBasis、coherent、g2 等直接可用，无需 import |
| 变量跨行共享 | 所有行共享同一命名空间 |
| 错误容错 | 出错不终止，继续执行后续行 |
| `import` 支持 | `import numpy as np` 等可用 |
| 多行语法 | dict、for、if、def、try 等复合语句 |
| Tab 补全 | `run <TAB>` 补全脚本路径；`cd <TAB>` 补全目录 |
| 嵌套调用 | `run another.qms` 可在脚本内调用其他脚本 |

#### agent 内置命令

| 命令 | 说明 |
|------|------|
| `<expression>` | 直接 Python 表达式求值 |
| `formula <latex>` | LaTeX 公式 → 终端 Unicode 显示 + PNG 保存 |
| `run <script.qms>` | 执行量子脚本 |
| `animate <var> [path]` | 生成波函数动画 |
| `plot wigner` | 绘制 Wigner 函数 |
| `cd [path]` | 切换工作目录（Tab 补全） |
| `pwd` | 打印当前目录 |
| `ls` | 列出当前目录文件 |
| `demo` | 运行 Fock 基演示 |
| `test` | 运行自检 |
| `help` | 显示帮助 |
| `quit` | 退出 |

#### 两种模式对比

| | Python .py | Agent .qms |
|---|---|---|
| **导入** | `from src.qm import *` | 自动预加载 |
| **执行** | `python file.py` | `python agent.py --run file.qms` |
| **Python 语法** | 全部 | 全部（逐行/多行解析） |
| **print()** | ✓ | ✓ |
| **matplotlib** | `plt.show()` | `matplotlib.use('Agg')` + `plt.savefig()` |
| **公式显示** | `print()` | `formula \latex` → Unicode |
| **Tab 补全** | IDE | 内置 readline |
| **调试** | pdb / IDE | 行号错误提示 |
| **适用场景** | 科研计算、制图、参数扫描 | 教学演示、快速验证、截图 |

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

### Python 示例

```python
>>> from src.qm import FockBasis
>>> fb = FockBasis(20)
>>> fb.a[:4, :4]
array([[0., 1., 0., 0.],
       [0., 0., √2, 0.],
       [0., 0., 0., √3],
       [0., 0., 0., 0.]])
```

### Agent 快捷方式

```python
# agent 和 .qms 中 fb 已预创建为 FockBasis(50)
⚛ > fb.a[:3,:3]
⚛ > fb.displacement(1.0)[:3,:3]
```

---

## 3. 量子态函数

### 3.1 `fock(N, n)` — Fock 态

```python
fock(N: int, n: int = 0) -> np.ndarray  # shape (N,)
```

返回 Fock 态 |n⟩。

```python
# Python 模式
>>> psi = fock(10, 3)
>>> psi
array([0., 0., 0., 1., 0., ...])

# Agent 模式
⚛ > psi3 = fock(20, 3)
⚛ > g2(psi3)
  0.6667           # g² = 1 - 1/n = 2/3 ✓
```

### 3.2 `coherent(N, alpha)` — 相干态

```python
coherent(N: int, alpha: complex) -> np.ndarray  # shape (N,)
```

返回相干态 |α⟩。|α|² 即为平均光子数。

```python
# Python 模式
>>> psi = coherent(20, 2.0 + 1.0j)
>>> np.linalg.norm(psi)  # 归一化验证
1.0
```

```
# Agent 模式
⚛ > psi = coherent(30, 3.0)
⚛ > mean_photon(psi)
  9.0000           # ⟨n⟩ = |α|² = 9 ✓
⚛ > g2(psi)
  1.0000           # Poisson ✓
```

### 3.3 `squeezed(N, zeta)` — 压缩真空

```python
squeezed(N: int, zeta: complex) -> np.ndarray  # shape (N,)
```

压缩参数 ζ = r e^{iθ}。平均光子数 ⟨n⟩ = sinh²(r)。

### 3.4 `thermal_dm(N, n_th)` — 热态

```python
thermal_dm(N: int, n_th: float) -> np.ndarray  # shape (N, N)
```

### 3.5 `cat(N, alpha, phi)` — 薛定谔猫态

```python
cat(N: int, alpha: complex, phi: float = 0.0) -> np.ndarray  # shape (N,)
```

|ψ⟩ ∝ |α⟩ + e^{iφ}|-α⟩。φ=0 为偶猫态，φ=π 为奇猫态。

### 3.6 态诊断函数

| 函数 | 返回 | 说明 |
|------|:---:|------|
| `fidelity(psi1, psi2)` | float | F = |⟨ψ₁|ψ₂⟩|² |
| `purity(rho)` | float | Tr[ρ²] |
| `photon_dist(state)` | array | P(n) 分布 |
| `is_dm(state)` | bool | 是否为密度矩阵 |

---

## 4. 算符工具

| 函数 | 说明 |
|------|------|
| `commutator(A, B)` | [A, B] = AB − BA |
| `expect(oper, state)` | ⟨O⟩，自动检测纯态/密度矩阵 |
| `variance(oper, state)` | ΔO² = ⟨O²⟩ − ⟨O⟩² |
| `is_hermitian(A)` | A = A†? |
| `is_unitary(U)` | U†U = I? |

---

## 5. 光子统计

| 函数 | 说明 |
|------|------|
| `mean_photon(state, fb)` | ⟨a†a⟩ |
| `g2(state, fb)` | g²(0) — 二阶关联 |
| `mandel_q(state, fb)` | Q = ⟨n⟩(g²−1) |

| 态 | g²(0) |
|---|:---:|
| 相干态 |α⟩ | 1.0 |
| 热态 (n̄) | 2.0 |
| Fock 态 |n⟩ | 1 − 1/n |
| 压缩真空 (r) | 3 + 1/sinh²(r) |

---

## 6. 波函数动力学 (TDSE)

| 类/函数 | 说明 |
|------|------|
| `WaveGrid(x_min, x_max, N)` | 一维空间网格 |
| `gaussian_wavepacket(grid, x0, p0, sigma)` | 高斯波包 |
| `evolve_ssfm(psi, grid, dt, t_max, snapshots)` | SSFM 求解 TDSE |
| `animate_wave(result, save_path)` | 生成 .mp4/.gif 动画 |

**Python 模式示例：**

```python
from src.qm import WaveGrid, gaussian_wavepacket, evolve_ssfm

grid = WaveGrid(-40, 40, 1024)
psi0 = gaussian_wavepacket(grid, x0=-8, p0=2.0, sigma=3.0)
result = evolve_ssfm(psi0, grid, dt=0.005, t_max=4.0, snapshots=80)
```

**Agent 示例：**

```
⚛ > grid = WaveGrid(-40, 40, 1024)
⚛ > psi0 = gaussian_wavepacket(grid, x0=-8, p0=2.0, sigma=3.0)
⚛ > res = evolve_ssfm(psi0, grid, dt=0.005, t_max=4.0, snapshots=80)
⚛ > animate res output/wave.mp4
```

### 特色势函数

| 函数 | 说明 |
|------|------|
| `double_well(grid, a, depth, sep)` | 双阱势 V₀[(x/a)²-1]²/2 |
| `periodic_potential(grid, amp, period)` | 余弦光晶格 A cos(2πx/λ) |
| `delta_barrier(grid, x0, strength)` | δ 势垒 (窄高斯近似) |
| `finite_well(grid, x0, width, depth)` | 有限深方阱 |
| `harmonic_oscillator_potential(grid, omega, mass)` | ½mω²x² |
| `step_potential(grid, x0, height)` | 阶跃势 |

```python
# Python 模式
from src.qm import WaveGrid, gaussian_wavepacket, evolve_ssfm, double_well
grid = WaveGrid(-15, 15, 1024)
V = double_well(grid, a=3.0, depth=8.0)
psi0 = gaussian_wavepacket(grid, x0=-3, p0=1.5, sigma=0.8)
result = evolve_ssfm(psi0, grid, V_func=lambda g: V, dt=0.005, t_max=20)
```

### PotentialBuilder — 构造复杂势函数

```python
# Python 模式
from src.qm import WaveGrid, PotentialBuilder

grid = WaveGrid(-20, 20, 1024)
V = (PotentialBuilder(grid)
     .harmonic(omega=0.5)                    # 谐振子底
     .well(x0=-4, depth=5, width=2)          # 左侧势阱
     .well(x0=+4, depth=5, width=2)          # 右侧势阱
     .barrier(x0=0, height=3, width=1)       # 中间势垒
     .periodic(amplitude=1, period=4, envelope_sigma=8)  # 光晶格+包络
     .build())

V.plot()       # 可视化预览
V.summary()    # 组件清单

result = evolve_ssfm(psi0, grid, V_func=V, dt=0.005, t_max=20)
```

Agent 模式同样语法，函数预加载无需 import。

---

## 7. 时间演化 (Fock 空间)

| 函数 | 说明 |
|------|------|
| `sesolve(H, psi0, tlist, e_ops)` | Schrödinger 方程 (精确对角化) |
| `mesolve(H, rho0, tlist, c_ops, e_ops)` | Lindblad 主方程 (RK4) |
| `steadystate(H, c_ops)` | 稳态 Liouvillian 求解 |

---

## 8. 相空间可视化

| 函数 | 说明 |
|------|------|
| `wigner(state, N_grid, xlim, ylim)` | Wigner 函数 (x, p, W) |
| `qfunc(state, N_grid, xlim, ylim)` | Husimi Q 函数 |
| `plot_wigner(x, p, W, save)` | Wigner 等高线图 |
| `plot_photon_dist(state, save)` | 光子数分布柱状图 |

---

## 9. shell 导航命令

| 命令 | 说明 |
|------|------|
| `pwd` | 打印当前工作目录 |
| `cd [path]` | 切换目录（Tab 补全，无参 = HOME） |
| `ls` | 列排显示当前目录内容（目录加 `/`） |

---

## 10. 脚本库参考

现有 `.qms` 脚本：

| 脚本 | 内容 | 输出 |
|------|------|------|
| `harmonic_oscillator.qms` | 谐振子 5 步分析 | 终端输出 |
| `core_formulas.qms` | 10 核心公式一览 | PNG × 10 |
| `wigner_gallery.qms` | 5 态 Wigner 对比 | 2×3 面板 PNG |
| `measurement_collapse.qms` | 位置测量坍缩动画 | 双面板 .mp4 |
| `free_particle.qms` | 自由粒子量子弥散 | .mp4 + τ 验证 |
| `heisenberg_uncertainty.qms` | 不确定性原理 | 4 面板动画 |
| `energy_collapse.qms` | 能量测量坍缩 (行波→驻波) | 双面板 .mp4 |
| `double_well.qms` | 双阱势量子隧穿 | .mp4 动画 |
| `pimc_demo.qms` | 路径积分 Monte Carlo | 基态对比图 |

```bash
python agent.py --run scripts/harmonic_oscillator.qms
python agent.py --run scripts/wigner_gallery.qms
python agent.py --run scripts/measurement_collapse.qms
```

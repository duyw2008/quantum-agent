# Quantum Agent 功能清单

## 公式终端显示

`formula` 命令将 LaTeX 数学公式转换为 Unicode 原生显示在终端中。

| 输入 | 终端显示 |
|------|---------|
| `formula i\hbar\frac{\partial}{\partial t}\Psi = \hat{H}\Psi` | iℏ(∂)/(∂t)Ψ = ĤΨ |
| `formula [\hat{x}, \hat{p}] = i\hbar` | [x̂, p̂] = iℏ |
| `formula \sigma_x \sigma_p \geq \frac{\hbar}{2}` | σ_x σ_p ≥ ℏ/2 |
| `formula H = \frac{p^2}{2m} + \frac{1}{2}m\omega^2x^2` | H = p²/2m + 1/2 m ω² x² |
| `formula \int_{-\infty}^{\infty} \psi^*\psi dx = 1` | ∫_-∞^∞ ψ*ψ dx = 1 |

内置 90+ LaTeX→Unicode 映射：希腊字母、数学算符、hat、上下标、积分、梯度、梯度平方等。
PNG 文件同步保存到 `output/formulas/` 作为高精度备份。

## Fock 基量子光学

### 算符 (FockBasis)
| 算符 | 符号 | 属性 |
|------|:---:|------|
| 湮灭 | â | `fb.a` |
| 产生 | â† | `fb.a_dag` |
| 数 | N̂ | `fb.n_op` |
| 坐标 | x̂ | `fb.x` |
| 动量 | p̂ | `fb.p` |
| 宇称 | Π̂ | `fb.parity` |
| 位移 | D̂(α) | `fb.displacement(α)` |
| 哈密顿量 | Ĥ | `fb.hamiltonian(ω)` |

### 量子态
| 函数 | 态 |
|------|------|
| `fock(N, n)` | \|n⟩ |
| `coherent(N, α)` | \|α⟩ |
| `squeezed(N, ζ)` | \|ζ⟩ |
| `thermal_dm(N, n̄)` | ρ_th |
| `cat(N, α, φ)` | \|α⟩+e^{iφ}\|-α⟩ |

### 工具函数
| 函数 | 说明 |
|------|------|
| `commutator(A,B)` | [A,B] |
| `expect(O,ρ)` | ⟨O⟩ |
| `variance(O,ρ)` | ΔO² |
| `g2(ρ)` | g²(0) |
| `mandel_q(ρ)` | Mandel Q |
| `fidelity` / `purity` | 态诊断 |

### 动力学
| 函数 | 说明 |
|------|------|
| `sesolve(H,ψ₀,t)` | Schrödinger 方程 |
| `mesolve(H,ρ₀,t,c_ops)` | Lindblad 主方程 |
| `steadystate(H,c_ops)` | 稳态求解 |

## 波函数动力学

| 函数 | 说明 |
|------|------|
| `WaveGrid(xmin,xmax,N)` | 空间网格 |
| `gaussian_wavepacket(...)` | 高斯波包 |
| `evolve_ssfm(ψ,g,dt,tmax)` | SSFM 演化 |
| `animate_wave(result,path)` | 动画生成 |

## 相空间可视化

| 函数 | 说明 |
|------|------|
| `wigner(ρ)` | Wigner 函数 |
| `qfunc(ρ)` | Husimi Q 函数 |
| `plot_wigner(x,p,W)` | Wigner 图 |
| `plot_photon_dist(ρ)` | 光子分布图 |

## Demo 动画

| Demo | 物理 |
|------|------|
| `free_particle.py` | 自由粒子量子弥散 |
| `heisenberg_uncertainty.py` | Δx·Δp ≥ ℏ/2 |
| `measurement_collapse.py` | 位置测量坍缩 |
| `momentum_collapse.py` | 动量测量坍缩 |
| `energy_collapse.py` | 能量测量坍缩 (驻波) |
| `double_slit.py` | 双缝干涉 (2D TDSE) |
| `quantum_eraser.py` | 量子擦除实验 |

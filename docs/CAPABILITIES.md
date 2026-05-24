# Quantum Agent 功能清单

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

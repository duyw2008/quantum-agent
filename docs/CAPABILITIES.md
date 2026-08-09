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

## PotentialBuilder — 链式势函数构造器

通过方法链组合任意复杂势函数：

```python
V = (PotentialBuilder(grid)
     .harmonic(omega=0.5)                 # 谐振子底
     .barrier(x0=0, height=6, width=0.8)  # 中心势垒
     .gaussian(x0=-5, height=-4, sigma=1.5)  # 高斯阱
     .periodic(amplitude=1.5, period=4,   # 光晶格+包络
               envelope_sigma=7)
     .delta(x0=-8, strength=3)            # δ 杂质
     .build())
```

### 基本构建块

| 方法 | 势函数 |
|------|--------|
| `.harmonic(omega, mass)` | ½ m ω² x² |
| `.barrier(x0, height, width)` | h (|x-x₀| < w/2) |
| `.well(x0, depth, width)` | -d (|x-x₀| < w/2) |
| `.gaussian(x0, height, sigma)` | h exp(-(x-x₀)²/2σ²) |
| `.periodic(amplitude, period, envelope)` | A cos(2πx/λ) + 可选高斯包络 |
| `.delta(x0, strength)` | g δ(x-x₀) 窄高斯近似 |
| `.step(x0, height)` | h (x > x₀) |
| `.linear(slope)` | s·x |
| `.custom(func, name)` | f(x) 自定义 |

### 快捷组合

| 方法 | 说明 |
|------|------|
| `.double_well(separation, depth, barrier_width)` | 双阱 + 中间势垒 |
| `.tunnel_junction(gap, height)` | 薄势垒隧穿结 |
| `.optical_lattice(amplitude, n_sites, envelope)` | n 个周期的余弦光晶格 |

### 代数操作 & 输出

| 方法 | 说明 |
|------|------|
| `.add(other)` | 叠加另一个势或数组 |
| `.multiply(factor)` | 整体缩放 |
| `.offset(shift)` | 整体平移 |
| `.build()` | 输出可调用函数 |
| `.plot(xlim, save)` | matplotlib 可视化 |
| `.summary()` | 打印组件清单 + 数值范围 |
| `.to_qms(path)` | 导出 .qms 脚本片段 |

## .qms 量子脚本

类似 MATLAB `.m` 文件，将一系列 agent 命令写入脚本批量执行。

```
# harmonic_oscillator.qms — 谐振子分析
formula [\hat{x}, \hat{p}] = i\hbar

alpha = 2.0 + 0.5j
psi = coherent(30, alpha)
mean_photon(psi)

C = commutator(fb.x, fb.p)
g2(psi)

cat_even = cat(30, 2.0, 0)
x, p, W = wigner(cat_even)
```

**执行方式：**

```bash
python agent.py --run scripts/harmonic_oscillator.qms   # 命令行
```

```
⚛ > run scripts/harmonic_oscillator.qms                  # 交互模式
```

**脚本特性：**
- `#` 行注释，跳过空行
- 变量跨行共享（同一命名空间）
- 错误不终止（继续执行后续命令）
- 支持 `import` / `print` 等 Python 内置
- 支持嵌套调用 `run <another.qms>`
- 工作目录自动切换到脚本所在目录

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
| `qft_scalar_field.py` | 标量场 φ̂(x), [φ̂,φ̂], D_F |
| `qft_lattice.py` | 格点 φ⁴ — E₀(λ), 关联函数, 能隙 |
| `qft_scattering.py` | Wick 定理, Feynman 图, 截面 |

## 量子场论 (QFT)

### 自由标量场 (field.py)

| 方法/属性 | 说明 |
|-----------|------|
| `ScalarField(mass, L, N_modes)` | 1+1D 标量场 φ̂(x) 模式展开 |
| `sf.field_matrix(x)` | φ̂(x) 的截断 Fock 矩阵表示 |
| `sf.commutator(x, y)` | ⟨0|[φ̂(x), φ̂(y)]|0⟩ |
| `sf.vacuum_fluctuation(x)` | ⟨0|φ̂²|0⟩ — 平移不变真空涨落 |
| `sf.vacuum_energy_density()` | 零点能密度 E₀/L |
| `sf.feynman_propagator(x, y, t)` | D_F(x-y, t) — Feynman 传播子 |
| `sf.propagator_profile(t, x_points)` | 传播子 vs 距离 |
| `sf.number_operator(mode_idx)` | 第 mode_idx 个模式的粒子数算符 |

### 格点 φ⁴ 理论 (lattice.py)

| 方法/属性 | 说明 |
|-----------|------|
| `LatticePhi4(N_sites, mass, coupling, N_fock)` | 1+1D 格点 φ⁴, 维度 = N_fock^{N_sites} |
| `lpt.hamiltonian(coupling)` | 完整哈密顿量矩阵 (动能+质量+φ⁴+梯度) |
| `lpt.diagonalize(coupling)` | 精确对角化 → (本征值, 本征矢) |
| `lpt.ground_state_energy(coupling)` | 基态能量 E₀(λ) |
| `lpt.energy_gap(coupling)` | 能隙 Δ = E₁ - E₀ |
| `lpt.correlation(i, j, coupling)` | 关联函数 ⟨φ_i φ_j⟩ |
| `lpt.particle_number_distribution(coupling)` | 每格点平均粒子数 ⟨N_j⟩ |
| `lpt.scan_coupling(couplings)` | 扫描 λ: E₀, Δ, 关联函数 |

### 重整化 (renormalization.py)

| 函数/类 | 说明 |
|----------|------|
| `self_energy_1loop(p², m, λ, Λ)` | 单圈自能 Π(p²) — ∝ λ∫d²k/(k²+m²) |
| `mass_counterterm(m, λ, Λ)` | δm = Π(0) (on-shell) |
| `coupling_counterterm(m, λ, Λ)` | δλ 顶点修正 (零动量 s=t=u=0) |
| `field_renormalization(m, λ, Λ)` | Z_φ = 1 + dΠ/dp²|_{p²=m²} |
| `beta_function(λ)` | β(λ) = 3λ²/(16π²) — 单圈重整化群流 |
| `running_coupling(λ₀, μ₀, μ)` | λ(μ) = λ₀ / [1 - (3λ₀/16π²)ln(μ/μ₀)] |

### U(1) 规范场 (gauge.py)

| 函数/类 | 说明 |
|----------|------|
| `GaugeField` | A_μ(x) 动量空间模式展开 |
| `polarization_vectors(k_μ, λ)` | 横模极化矢量 ε_μ(k, λ=1,2) |
| `photon_propagator(k², ξ)` | -i[g_μν - (1-ξ)k_μk_ν/k²]/k² |
| `ward_identity_check(amp_fn, k_μ)` | k_μ M^μ = 0 验证 |

### Dirac 旋量 (dirac.py)

| 函数/类 | 说明 |
|----------|------|
| `GammaMatrices(rep='dirac')` | γ^μ 矩阵 (Dirac / Chiral 表象) |
| `gm.gamma5` | γ⁵ = iγ⁰γ¹γ²γ³ |
| `gm.sigma_munu` | Σ^{μν} = (i/2)[γ^μ, γ^ν] |
| `DiracSpinor` | u(p,s) 正粒子, v(p,s) 反粒子旋量 |
| `dirac_slash(p_μ)` | p̸ = γ^μ p_μ |
| `spin_sum_u(p_μ, m)` | Σ uū = p̸ + m |
| `spin_sum_v(p_μ, m)` | Σ vv̄ = p̸ - m |
| `bilinear(ψ̄, Γ, ψ)` | 双线性协变量 ŌΓO |

### QED 散射 (qed.py)

| 函数/类 | 说明 |
|----------|------|
| `mandelstam(p1, p2, p3, p4)` | s, t, u 不变量 |
| `compton_cross_section(ω, θ, m_e)` | Klein-Nishina 公式 dσ/dΩ |
| `pair_annihilation_cross_section(s, θ, m_e)` | e⁺e⁻ → μ⁺μ⁻ dσ/dΩ |
| `moller_cross_section(s, θ, m_e)` | e⁻e⁻ → e⁻e⁻ dσ/dΩ |
| `running_alpha(α₀, Q², m_e)` | 跑动精细结构常数 α(Q²) |

### 有效势 & 对称破缺 (effective_potential.py)

| 函数/类 | 说明 |
|----------|------|
| `OneLoopEffectivePotential(m, λ, μ)` | V_eff(φ_c) = V_0 + (ħ/64π²)V''²ln(V''/μ²) |
| `ColemanWeinberg(λ, μ)` | m=0 极限下的有效势 |
| `find_minimum(V_eff, φ_range)` | 找有效势最小值 φ_min |
| `SymmetryBreaking(model)` | 序参量, Goldstone 定理, Higgs 机制 |

### QFT 路径积分 Monte Carlo (lattice_qft.py)

| 函数/类 | 说明 |
|----------|------|
| `LatticePhi4MC(N_x, N_τ, m, λ)` | 2D 格点 φ⁴ Metropolis 采样 |
| `lqft.thermalize(n_sweeps)` | 热化 N 个完整路径扫描 |
| `lqft.correlation_function(dx)` | 空间关联 ⟨φ(0)φ(dx)⟩ |
| `lqft.two_point_function(τ)` | 虚时关联 ⟨φ(0)φ(τ)⟩ → 质量提取 |
| `lqft.effective_mass(corr_2pt)` | m_eff = ln(C(τ)/C(τ+1)) |
| `lqft.susceptibility()` | χ = ⟨(Σφ)²⟩ - ⟨Σφ⟩² (相变探测) |

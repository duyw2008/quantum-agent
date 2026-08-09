# QFT 模块全面扩展 — 实现计划

> 目标: 补全量子场论计算能力 — 重整化闭环 + U(1)规范场 + Dirac费米子 + QED过程 + QFT路径积分 + 对称破缺

**新增 6 个模块 + 1 个更新:**

```
src/qft/
├── renormalization.py     [NEW] φ⁴ 1-loop renormalization
├── gauge.py               [NEW] U(1) gauge field
├── dirac.py               [NEW] Dirac spinors
├── qed.py                 [NEW] QED cross sections
├── effective_potential.py [NEW] 1PI effective action + SSB
├── lattice_qft.py         [NEW] QFT path integral MC
├── __init__.py            [UPDATE] export all new modules
└── agent.py               [UPDATE] register new names
```

---

## Task 1: renormalization.py — φ⁴ 单圈重整化

**内容:**
- `Phi4FeynmanRules` — Feynman 规则: propagator i/(p²-m²+iε), vertex -iλ
- `regularize_cutoff(integrand, Lambda)` — 动量截断正规化
- `regularize_dim_reg(integrand, d, mu)` — 维数正规化 (简化)
- `self_energy_1loop(p², m, lam, Lambda)` — Π(p²) 单圈自能
- `vertex_correction_1loop(s, t, u, m, lam, Lambda)` — Γ⁴ 单圈顶点修正
- `mass_counterterm(m, lam, Lambda, scheme)` — δm (on-shell / MS-bar)
- `coupling_counterterm(m, lam, Lambda, scheme)` — δλ
- `field_renormalization(m, lam, Lambda)` — Z_φ
- `beta_function(lam, mu)` — β(λ) 重整化群流 (1-loop)
- `running_coupling(lam0, mu0, mu)` — λ(μ) 跑动耦合

**验证:** 格点 φ⁴ 在不同 N_fock / N_sites 下基态能量是否 match 重整化后微扰结果

## Task 2: gauge.py — U(1) 规范场

**内容:**
- `GaugeField` — A_μ(x) 在动量空间的模式展开 (光子场)
- 极化矢量 ε_μ(k, λ) — 两个横模
- 规范固定: 选择 ξ=1 (Feynman gauge)
- `photon_propagator(k, xi)` — 光子传播子 -i[g_μν - (1-ξ)k_μk_ν/k²]/k²
- `field_strength(F_munu)` — F_μν = ∂_μ A_ν - ∂_ν A_μ
- `gauge_transform(A_mu, alpha)` — 规范变换
- `ward_identity_check(amplitude, k_mu)` — Ward 恒等式验证

## Task 3: dirac.py — Dirac 旋量

**内容:**
- `GammaMatrices` — γ^μ 在 Dirac/ chiral / Majorana 表象
- `gamma5` — γ⁵ = iγ⁰γ¹γ²γ³
- `DiracSpinor` — u(p,s), v(p,s) 正/反粒子旋量
- `spin_sum_u(p, m)` — Σ_s u(p,s)ū(p,s) = p̸ + m
- `spin_sum_v(p, m)` — Σ_s v(p,s)v̄(p,s) = p̸ - m
- `dirac_equation_check(spinor, p, m)` — 验证 (p̸ - m)u = 0
- `bilinear(scalar/pseudoscalar/vector/axial/tensor)` — 双线性协变

## Task 4: qed.py — QED 过程

**内容:**
- `QEDVertex` — -ieγ^μ (QED 顶点)
- `compton_amplitude(s, t, u, m)` — Compton 散射振幅
- `compton_cross_section(omega, theta, m)` — Klein-Nishina 公式
- `pair_annihilation_cross_section(s, theta, m)` — e⁺e⁻ → μ⁺μ⁻
- `moller_scattering_cross_section(s, t, m)` — Møller 散射
- `bhabha_cross_section(s, theta, m)` — Bhabha 散射
- `running_alpha(alpha0, Q², m_e)` — 跑动精细结构常数

## Task 5: effective_potential.py — 有效势 + 对称破缺

**内容:**
- `OneLoopEffectivePotential` — V_eff(φ_c) = V_0 + (ħ/64π²) V''² ln(V''/μ²)
- `ColemanWeinberg(m, lam, mu)` — Coleman-Weinberg 有效势
- `find_minimum(V_eff, phi_range)` — 找有效势最小值
- `SymmetryBreaking` — 自发对称破缺分析
  - order_parameter(φ_min) — 序参量
  - Goldstone theorem verification
  - Higgs mechanism (耦合到 U(1) gauge)

## Task 6: lattice_qft.py — QFT 格点路径积分 MC

**内容:**
- `LatticePhi4MC` — 2D 格点上的 φ⁴ 路径积分 Monte Carlo
- 场构型 `phi[i, tau]` → Metropolis 更新
- 两点关联函数 `correlator(i, j)`
- 有效质量提取 `effective_mass(correlator, d)`
- 相变探测 `susceptibility(phi_samples)` — 磁化率峰值

## Task 7: 更新 __init__.py + agent.py

- `__init__.py`: 导出所有新模块的关键类和函数
- `agent.py`: 在 `_completions` 和 `ns` 中注册新名称
- 添加 `help` 公式条目

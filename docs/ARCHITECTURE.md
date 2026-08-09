# Quantum Agent 架构

## 概述

Quantum Agent 是一个模块化量子物理平台，涵盖三个层次：量子力学 (Fock 基 + 波函数 + 自旋/纠缠)、量子场论 (自由场 → 微扰 → 格点 → 重整化 → QED)、路径积分 Monte Carlo。

## 模块架构

```python
agent.py                     # CLI 交互界面 (含 LaTeX→Unicode 公式显示)
src/
├── qm/                      # 量子力学核心
│   ├── basis.py             # FockBasis — 算符工厂
│   ├── states.py            # 量子态 (coherent, squeezed, thermal, cat...)
│   ├── operators.py         # 算符工具 (commutator, expect, g2, mandel_q)
│   ├── dynamics.py          # 动力学 (sesolve, mesolve, steadystate)
│   ├── wave.py              # 波函数 (WaveGrid, SSFM, PotentialBuilder)
│   ├── wave2d.py            # 2D TDSE (双缝, 量子擦除)
│   ├── spin.py              # Pauli 矩阵, Bloch 球, CHSH 检验
│   └── multipartite.py      # 纠缠度量 (concurrence, partial_trace, vN entropy)
├── viz/                     # 可视化
│   └── wigner_plot.py       # Wigner, Qfunc, 光子分布图
├── qft/                     # 量子场论 (10 个模块)
│   ├── field.py             # ScalarField — 自由标量场 φ̂(x)
│   ├── lattice.py           # LatticePhi4 — 格点 φ⁴ 精确对角化
│   ├── scattering.py        # 微扰散射 — Wick 定理, Feynman 振幅, Dyson 级数
│   ├── path_integral.py     # 路径积分 MC (1D QM)
│   ├── renormalization.py   # φ⁴ 单圈重整化 — Π(p²), 顶点修正, β 函数
│   ├── gauge.py             # U(1) 规范场 — A_μ 模式展开, 光子传播子
│   ├── dirac.py             # Dirac 旋量 — γ 矩阵, u/v 旋量, 自旋求和
│   ├── qed.py               # QED 过程 — Compton, 对湮灭, Møller/Bhabha 散射
│   ├── effective_potential.py # 1PI 有效势 — Coleman-Weinberg, 自发对称破缺
│   └── lattice_qft.py       # QFT 格点路径积分 MC — 2D φ⁴ 场构型采样
├── demos/                   # 13 个物理动画
└── scripts/                 # 13 个 .qms 量子脚本
```

## 物理栈

```
QED 散射                (qed)       Compton, e⁺e⁻→μ⁺μ⁻, Møller
    ↑
Dirac 费米子            (dirac)     u/v旋量, γ矩阵, 自旋求和
    ↑                       
U(1) 规范场             (gauge)     A_μ, ε_μ(k,λ), 光子传播子
    ↑                                    
重整化                  (renorm)    Π(p²), Γ⁴, δm, δλ, β(λ)
    ↑                                    
微扰 φ⁴                (scattering) Wick, Feynman 图, Dyson 级数
    ↑                                    
格点 φ⁴ 精确对角化       (lattice)   E₀(λ), 关联函数, 能隙
    ↑                                    
自由标量场              (field)     φ̂(x), [φ̂,φ̂], D_F, 真空涨落
    ↑
有效势 + SSB            (eff_pot)   V_eff, Coleman-Weinberg, Goldstone
    ↑
QFT 路径积分 MC         (lattice_qft) 2D 格点采样, 关联函数, 有效质量
```

## 数据流

```python
用户输入 → agent.py
  ├── 直接求值 → eval/exec → Python 表达式求值
  ├── formula    → _latex_to_unicode() → 终端 Unicode 显示 + PNG 保存
  ├── animate    → 波函数动画生成
  ├── plot       → Wigner 图
  ├── demo       → FockBasis 演示
  └── test       → 自检

calc 命名空间:
  QM: FockBasis, coherent, squeezed, thermal_dm, cat,
      expect, g2, mandel_q, commutator,
      sesolve, mesolve, steadystate,
      WaveGrid, gaussian_wavepacket, evolve_ssfm,
      wigner, qfunc, PotentialBuilder,
      concurrence, partial_trace, entropy_vn, bell_states
  QFT: ScalarField, LatticePhi4, LatticePhi4MC,
       wick_expand, feynman_amplitude_phi4_2to2, transition_probability,
       self_energy_1loop, beta_function, running_coupling,
       GaugeField, photon_propagator, polarization_vectors,
       GammaMatrices, DiracSpinor, dirac_slash,
       compton_cross_section, pair_annihilation_cross_section,
       OneLoopEffectivePotential, ColemanWeinberg, SymmetryBreaking
```

## 设计原则

1. **QuTiP 兼容**: 函数签名参照 QuTiP，降低学习成本
2. **懒加载**: 模块按需导入，calc 首次调用时加载
3. **Fock 基优先**: 所有算符默认在截断数态基中表示
4. **分离关注**: qm/ 负责计算，viz/ 负责可视化
5. **终端友好**: 公式直接以 Unicode 数学符号显示，无需外部图片查看器

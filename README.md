# Quantum Agent ⚛️

# Quantum Agent ⚛️

> 量子物理智能体 — QM (Fock + 波函数 + 自旋/纠缠) × QFT (自由场 → 微扰 → 格点 → 重整化 → QED) × Monte Carlo

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

交互式量子物理平台：Fock 基量子光学 × 波函数动力学 × QFT 全过程 (自由场→重整化→QED) × QFT 路径积分 MC × .qms 脚本。

## 快速开始

```bash
pip install numpy scipy matplotlib
python agent.py          # 交互模式
python agent.py --demo   # Fock 基演示
python agent.py --test   # 自检
python agent.py --list   # 列出所有 demo
python agent.py --run  scripts/harmonic.qms   # 执行量子脚本 (.qms)
```

## 使用 — 直接输入表达式

```python
⚛ > 1+1
2
⚛ > psi = coherent(20, 2.0)
⚛ > g2(psi)
1.0
⚛ > sf = ScalarField(mass=1.0)
⚛ > sf.commutator(0, 2)
⚛ > beta = beta_function(0.5)          # β(λ) 重整化群流
⚛ > dsigma = compton_cross_section(10, 0.5, 0.511)  # Klein-Nishina
⚛ > gm = GammaMatrices('dirac')        # γ 矩阵
⚛ > formula i\\hbar\\frac{\\partial}{\\partial t}\\Psi = \\hat{H}\\Psi
  iℏ(∂)/(∂t)Ψ = ĤΨ
⚛ > help
⚛ > quit
```

## Demo 动画

| Demo | 物理 | 运行 |
|------|------|------|
| 不确定性原理 | Δx·Δp ≥ ℏ/2 | `python demos/heisenberg_uncertainty.py` |
| 自由弥散 | 波包展宽 | `python demos/free_particle.py` |
| 位置坍缩 | 测量→Δx↓→弥散100× | `python demos/measurement_collapse.py` |
| 动量坍缩 | 频率可视化 | `python demos/momentum_collapse.py` |
| 能量坍缩 | 驻波 cos(kx) | `python demos/energy_collapse.py` |
| 双缝干涉 | 2D TDSE | `python demos/double_slit.py` |
| 量子擦除 | 相干 vs 非相干 | `python demos/quantum_eraser.py` |
| 标量场 | φ̂(x), [φ̂,φ̂], D_F | `python demos/qft_scalar_field.py` |
| 格点 φ⁴ | E₀(λ), 关联函数, 能隙 | `python demos/qft_lattice.py` |
| Feynman 图 | Wick 定理, 振幅, 截面 | `python demos/qft_scattering.py` |
| Bloch 球 | Rabi, 退相干, Larmor | `python demos/bloch_evolution.py` |

## QFT 能力总览

```
层次               模块                 核心功能
─────────────────────────────────────────────────────
QED 散射          qed.py               Klein-Nishina, e⁺e⁻→μ⁺μ⁻, Møller
Dirac 费米子       dirac.py              γ 矩阵, u/v 旋量, 自旋求和
U(1) 规范场        gauge.py             A_μ, ε_μ(k,λ), Ward 恒等式
重整化             renormalization.py    Π(p²), Γ⁴, δm, δλ, β(λ)
微扰 φ⁴            scattering.py        Wick 定理, Feynman 振幅, Dyson 级数
格点 φ⁴ 对角化      lattice.py           E₀(λ), 关联函数, 能隙
自由标量场          field.py             φ̂(x), [φ̂,φ̂], D_F, 真空涨落
有效势 + SSB        effective_potential.py  V_eff, Coleman-Weinberg
QFT 路径积分 MC     lattice_qft.py        2D 格点采样, 有效质量, 相变
```

## 模块架构

```
quantum_agent/
├── agent.py              # CLI (表达式、readline、tab补全、公式终端显示)
├── src/
│   ├── qm/               # 量子力学 (FockBasis, states, dynamics, wave, spin, multipartite)
│   ├── viz/              # Wigner, Qfunc, 光子分布
│   └── qft/              # 量子场论 (10 模块: field → lattice → scattering → renormalization
│                         #         → gauge → dirac → qed → effective_potential → lattice_qft)
├── demos/                # 11 个物理动画
├── scripts/              # 13 个 .qms 脚本
├── docs/                 # 10 份文档
└── output/               # 动画、图片、公式 PNG
```

## 文档

- [MATHEMATICS.md](docs/MATHEMATICS.md) — 完整数学模型
- [NUMERICAL_METHODS.md](docs/NUMERICAL_METHODS.md) — 数值方法详解 (SSFM/对角化/RK4/Wigner/PIMC)
- [USER_GUIDE.md](docs/USER_GUIDE.md) — 函数参考 + 实例
- [PHYSICS.md](docs/PHYSICS.md) — 物理基础
- [ARCHITECTURE.md](docs/ARCHITECTURE.md) — 系统架构
- [CAPABILITIES.md](docs/CAPABILITIES.md) — 功能清单
- [LEARNING_PATH.md](docs/LEARNING_PATH.md) — 学习路径（4阶段16实验）
- [TUTORIAL.md](docs/TUTORIAL.md) — 教程
- [KNOWLEDGE_HANDBOOK.md](docs/KNOWLEDGE_HANDBOOK.md) — 知识手册（Fock 缝合 + 量子因果）
- [CHANGELOG.md](docs/CHANGELOG.md) — 版本历史

## License

MIT

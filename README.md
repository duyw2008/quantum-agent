# Quantum Agent ⚛️

> 量子力学智能体 — Fock 空间 + 波函数 + 量子场论 + 公式终端显示

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

交互式量子力学平台：Fock 基量子光学、波函数动力学、格点量子场论、.qms 脚本批量执行。

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
⚛ > formula i\hbar\frac{\partial}{\partial t}\Psi = \hat{H}\Psi
  iℏ(∂)/(∂t)Ψ = ĤΨ
⚛ > x, p, W = wigner(psi)
⚛ > plot_wigner(x, p, W)
⚛ > sf = ScalarField(mass=1.0)
⚛ > sf.commutator(0, 2)
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
| 标量场 | φ̂(x), 对易子, 传播子 | `python demos/qft_scalar_field.py` |
| 格点 φ⁴ | 基态能量, 关联函数 | `python demos/qft_lattice.py` |
| Feynman 图 | Wick 定理, 截面 | `python demos/qft_scattering.py` |

## 模块架构

```
quantum_agent/
├── agent.py              # CLI (表达式、readline、tab补全、公式终端显示)
├── src/
│   ├── qm/               # 量子力学 (FockBasis, states, dynamics, wave)
│   ├── viz/              # Wigner, Qfunc, 光子分布
│   └── qft/              # 量子场论 (ScalarField, LatticePhi4, scattering)
├── demos/                # 10 个物理动画
├── scripts/              # .qms 脚本 (类 MATLAB .m)
├── docs/                 # 7 份文档
└── output/               # 动画、图片、公式 PNG
```

## 文档

- [MATHEMATICS.md](docs/MATHEMATICS.md) — 完整数学模型
- [USER_GUIDE.md](docs/USER_GUIDE.md) — 函数参考 + 实例
- [PHYSICS.md](docs/PHYSICS.md) — 物理基础
- [ARCHITECTURE.md](docs/ARCHITECTURE.md) — 系统架构
- [CAPABILITIES.md](docs/CAPABILITIES.md) — 功能清单
- [TUTORIAL.md](docs/TUTORIAL.md) — 教程
- [CHANGELOG.md](docs/CHANGELOG.md) — 版本历史

## License

MIT

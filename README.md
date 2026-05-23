# Quantum Agent ⚛️

> QuTiP 风格量子力学智能体 — Fock 空间 + 波函数动力学 + 动画演示

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

交互式量子力学平台：Fock 基量子光学 + 一维波函数动力学 + 2D 含时薛定谔方程。

## 快速开始

```bash
pip install numpy scipy matplotlib
python agent.py          # 交互模式
python agent.py --demo   # Fock 基演示
python agent.py --test   # 自检
```

## Demo 动画

| Demo | 物理 | 运行 |
|------|------|------|
| 双缝干涉 | 2D TDSE, 干涉条纹, γ校正 | `python demos/double_slit.py` |
| 量子擦除 | 相干 vs 非相干, 中途切换 | `python demos/quantum_eraser.py` |
| 不确定性原理 | Δx·Δp ≥ ℏ/2, 4面板 | `python demos/heisenberg_uncertainty.py` |
| 位置坍缩 | 测量→Δx↓→弥散100× | `python demos/measurement_collapse.py` |
| 动量坍缩 | 频率可视化, Δp↓→Δx↑18× | `python demos/momentum_collapse.py` |
| 自由弥散 | 高斯波包展宽 | `python demos/free_particle.py` |

## Agent 命令

```
calc <expr>           计算表达式 (np, FockBasis, coherent, wigner, ...)
calc <var> = <expr>   变量赋值
calc vars             查看变量
demo                  Fock 基演示
test                  自检 (对易子/g²/纯度/...)
help                  帮助
```

## 模块架构

```
quantum_agent/
├── agent.py              # CLI 交互界面
├── src/
│   ├── qm/               # 量子力学核心
│   │   ├── basis.py      # FockBasis: a, a†, x, p, N, 位移, 宇称
│   │   ├── states.py     # fock, coherent, squeezed, thermal, cat
│   │   ├── operators.py  # commutator, expect, g2, mandel_q
│   │   ├── dynamics.py   # sesolve, mesolve, steadystate
│   │   └── wave.py       # WaveGrid, SSFM, gaussian, animate_wave
│   └── viz/              # 可视化
│       └── wigner_plot.py # Wigner, Qfunc, 光子分布图
├── demos/                # 6 个物理动画
└── docs/                 # 文档
```

## 文档

- [MATHEMATICS.md](docs/MATHEMATICS.md) — 完整数学模型 (7 章)
- [USER_GUIDE.md](docs/USER_GUIDE.md) — 函数参考 + 实例 (9 章)

## License

MIT

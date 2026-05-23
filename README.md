# Quantum Agent ⚛️

> QuTiP 风格量子力学智能体 — Fock 空间计算与可视化

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

交互式量子力学函数库，参照 [QuTiP](https://qutip.org) 风格，在 Fock 基上进行量子态构建、算符代数、动力学演化和相空间可视化。

## 快速开始

```bash
python3 -m venv venv && source venv/bin/activate
pip install numpy scipy matplotlib

python agent.py          # 交互模式
python agent.py --demo   # 运行演示
python agent.py --test   # 运行自检
```

## Agent 命令

```
calc <expr>           计算表达式
calc <var> = <expr>   变量赋值
calc vars             查看变量
demo                  运行演示
test                  运行自检
help                  帮助
```

## 示例

```
⚛ > calc psi = coherent(20, 2.0)
⚛ > calc g2(psi)
  1.0
⚛ > calc x, p, W = wigner(psi)
⚛ > calc plot_wigner(x, p, W)
```

## 模块架构

```
quantum_agent/
├── agent.py              # CLI 交互界面
├── src/
│   ├── qm/               # 量子力学核心
│   │   ├── basis.py      # FockBasis: a, a†, x, p, N, 位移, 宇称
│   │   ├── states.py     # fock, coherent, squeezed, thermal_dm, cat
│   │   ├── operators.py  # commutator, expect, g2, mandel_q
│   │   └── dynamics.py   # sesolve, mesolve, steadystate
│   └── viz/              # 可视化
│       └── wigner_plot.py # Wigner, Qfunc, 光子分布图
└── docs/
    └── MATHEMATICS.md    # 完整数学模型
```

## 验证结果

| 物理量 | 预期 | 结果 |
|--------|:----:|:----:|
| 相干态 g²(0) | 1.0 | 1.0000 ✓ |
| 热态 g²(0) | 2.0 | 2.0000 ✓ |
| 压缩态 ⟨n⟩ | sinh²(r) | 精确匹配 ✓ |
| [x̂, p̂] ≈ iħ | 0 | 6.75×10⁻¹⁶ ✓ |

## 文档

- [MATHEMATICS.md](docs/MATHEMATICS.md) — 完整数学模型
- [USER_GUIDE.md](docs/USER_GUIDE.md) — 函数参考与应用实例

## License

MIT

# Quantum Agent 架构

## 概述

Quantum Agent 是一个模块化量子力学平台，包含两套计算引擎：Fock 基量子光学和一维/二维波函数动力学。

## 模块架构

```
agent.py                     # CLI 交互界面
src/
├── qm/                      # 量子力学核心
│   ├── basis.py             # FockBasis — 算符工厂
│   ├── states.py            # 量子态 (coherent, squeezed, thermal, cat...)
│   ├── operators.py         # 算符工具 (commutator, expect, g2, mandel_q)
│   ├── dynamics.py          # 动力学 (sesolve, mesolve, steadystate)
│   └── wave.py              # 波函数 (WaveGrid, SSFM, animate)
└── viz/                     # 可视化
    └── wigner_plot.py       # Wigner, Qfunc, 光子分布图
```

## 数据流

```
用户输入 → agent.py
  ├── calc → eval/exec → Python 表达式求值
  ├── demo → FockBasis 演示
  └── test → 自检

calc 命名空间:
  FockBasis, coherent, squeezed, thermal_dm, cat,
  expect, g2, mandel_q, commutator,
  sesolve, mesolve, steadystate,
  WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave,
  wigner, qfunc, plot_wigner, plot_photon_dist
```

## 设计原则

1. **QuTiP 兼容**: 函数签名参照 QuTiP，降低学习成本
2. **懒加载**: 模块按需导入，calc 首次调用时加载
3. **Fock 基优先**: 所有算符默认在截断数态基中表示
4. **分离关注**: qm/ 负责计算，viz/ 负责可视化

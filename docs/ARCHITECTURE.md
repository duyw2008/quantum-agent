# Quantum Agent 架构

## 概述

Quantum Agent 是一个模块化量子力学平台，包含两套计算引擎：Fock 基量子光学和一维/二维波函数动力学。

## 模块架构

```python
agent.py                     # CLI 交互界面 (含 LaTeX->Unicode 公式显示)
src/
├── qm/                      # 量子力学核心
│   ├── basis.py             # FockBasis — 算符工厂
│   ├── states.py            # 量子态 (coherent, squeezed, thermal, cat...)
│   ├── operators.py         # 算符工具 (commutator, expect, g2, mandel_q)
│   ├── dynamics.py          # 动力学 (sesolve, mesolve, steadystate)
│   └── wave.py              # 波函数 (WaveGrid, SSFM, animate)
├── viz/                     # 可视化
│   └── wigner_plot.py       # Wigner, Qfunc, 光子分布图
├── scripts/                 # .qms 量子脚本 (类似 MATLAB .m)
├── qft/                     # 量子场论
    ├── field.py             # ScalarField — 标量场
    ├── lattice.py           # LatticePhi4 — 格点 φ⁴
    ├── scattering.py        # 散射 / Feynman 图
    └── path_integral.py   # 路径积分 Monte Carlo (Phase 4)
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
  FockBasis, coherent, squeezed, thermal_dm, cat,
  expect, g2, mandel_q, commutator,
  sesolve, mesolve, steadystate,
  WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave,
  wigner, qfunc, plot_wigner, plot_photon_dist,
  ScalarField, LatticePhi4, wick_expansion, feynman_amplitude
```

## 设计原则

1. **QuTiP 兼容**: 函数签名参照 QuTiP，降低学习成本
2. **懒加载**: 模块按需导入，calc 首次调用时加载
3. **Fock 基优先**: 所有算符默认在截断数态基中表示
4. **分离关注**: qm/ 负责计算，viz/ 负责可视化
5. **终端友好**: 公式直接以 Unicode 数学符号显示，无需外部图片查看器

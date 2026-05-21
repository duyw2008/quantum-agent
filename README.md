# Quantum Agent ⚛️

> 交互式一维量子力学计算与可视化智能体

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![Tests](https://img.shields.io/badge/tests-35%2F35%20passed-brightgreen)](tests/run_all.py)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Quantum Agent 是一个模块化的量子力学计算平台，将波函数演化、矩阵力学和可视化深度融合：

- **波函数动态演化** — 求解含时薛定谔方程，两种数值方法可选
- **矩阵力学** — 算符代数、对易关系验证、本征值问题（numpy + sympy）
- **11 种势函数** — 从无限深势阱到 Morse 分子势，覆盖量子力学核心场景
- **动画可视化** — 暗色/亮色双主题，MP4/GIF 输出
- **交互式 Agent** — CLI 命令模式，5 个完整 Demo

## 📦 快速开始

```bash
# 安装依赖
python3 -m venv venv && source venv/bin/activate
pip install numpy scipy matplotlib

# 运行 agent（交互模式）
python agent.py

# 运行所有 demo（生成动画到 output/）
python agent.py --demo all

# 运行所有测试
python tests/run_all.py
# → 35/35 tests passed ✓
```

## 🧭 交互命令

```
⚛ > evolve harmonic,omega=2.0 x0=1.0 p0=5.0 t_max=10    # 谐振子波包
⚛ > evolve barrier,height=10.0,width=0.5 x0=-3.0 p0=4.0  # 势垒隧穿
⚛ > eigenstates double_well,separation=3.0,depth=8.0 5    # 本征态
⚛ > matrix report                                         # 矩阵力学报告
⚛ > matrix comm x p                                       # 验证 [x̂, p̂] = iℏ
⚛ > demo all                                              # 运行所有演示
```

## 🔬 物理能力

### 势函数 (11 种)

| 类型 | 参数 | 解析能级 | 物理场景 |
|------|------|:---:|------|
| `infinite_well` | width | ✓ | 量子约束基础模型 |
| `harmonic` | omega, mass | ✓ | 谐振子、相干态 |
| `barrier` / `finite_well` | height/depth, width | — | 隧穿、散射 |
| `double_well` | separation, depth | 数值 | 隧穿劈裂、量子振荡 |
| `morse` | depth, alpha, x0 | ✓ | 双原子分子振动 |
| `coulomb_1d` | Z, softening | 数值 | 一维氢原子模型 |
| `periodic` | amplitude, k | Bloch | 晶体能带结构 |
| `step` | height, x0 | — | 阶梯散射 |
| `zero` | — | 连续谱 | 自由粒子 |
| `custom` | func, name | — | 完全自定义 |

### TDSE 求解器

| 方法 | 复杂度 | 精度 | 稳定性 | 适用场景 |
|------|:---:|:---:|:---:|------|
| **Split-Step Fourier** | O(N log N) | O(Δt³) | 条件稳定 | 光滑势、长时演化 |
| **Crank-Nicolson** | O(N) | O(Δt²) | 无条件稳定 | 刚性势、精确守恒 |

### 矩阵力学

- 算符：â, â†, x̂, p̂, N̂, Ĥ（数态基截断表示）
- 对易子：自动计算 [Â, B̂] 并验证正则对易关系
- 本征值/本征态：完整对角化
- 时间演化：Û(t) = exp(-iĤt/ħ)，海森堡绘景
- 符号力学：sympy 符号代数（角动量、Pauli 自旋、BCH 公式）

## 📂 项目结构

```
quantum_agent/
├── agent.py              # CLI 交互界面
├── config.yaml           # 全局配置 (ℏ, m, 网格默认值...)
├── src/
│   ├── core/             # 物理核心
│   │   ├── potentials.py     # 11 种势函数 + 工厂
│   │   ├── wave_function.py  # WaveFunction + Grid
│   │   └── schrodinger.py    # SSFM + Crank-Nicolson 求解器
│   ├── matrix/           # 矩阵力学
│   │   ├── operators.py      # NumericOperatorSystem
│   │   └── symbolic.py       # SymbolicQuantum + MatrixMechanics
│   └── viz/              # 可视化
│       ├── animate.py        # 动画生成 + 静态图
│       └── static.py         # 能谱、矩阵热图、相空间
├── tests/                # 35 个测试
│   ├── test_potentials.py    # 12 tests
│   ├── test_wave_function.py # 6 tests
│   ├── test_schrodinger.py   # 6 tests
│   ├── test_matrix.py        # 7 tests
│   ├── test_viz.py           # 4 tests
│   └── run_all.py            # 测试运行器
├── demos/                # 5 个演示
│   ├── harmonic_oscillator.py
│   ├── infinite_well.py
│   ├── potential_barrier.py
│   ├── matrix_mechanics.py
│   └── double_well.py
├── docs/                 # 文档
│   ├── PHYSICS.md        # 物理公式和数学背景
│   ├── ARCHITECTURE.md   # 系统架构和数据流
│   ├── CAPABILITIES.md   # 完整功能清单
│   ├── TUTORIAL.md       # 交互 + Python API 教程
│   └── CHANGELOG.md      # 版本历史
└── output/               # 生成的动画 (MP4/GIF) 和图表 (PNG)
```

## 📐 物理理论基础

详见 [docs/PHYSICS.md](docs/PHYSICS.md)，涵盖：

- 含时薛定谔方程的 SSFM 和 Crank-Nicolson 数值推导
- 正则对易关系 [x̂, p̂] = iħ 在截断数态基中的表示
- 谐振子解析能级 Eₙ = ħω(n + ½) 和相干态动力学
- WKB 隧穿概率近似和双势阱能级劈裂
- 不确定度原理 Δx·Δp ≥ ħ/2 的数值验证

## 🧪 Demo 验证结果

| Demo | 物理现象 | 关键结果 |
|------|----------|----------|
| 矩阵力学 | [x̂,p̂], 能谱, HUP | [x̂,p̂]=iħ ✓, Eₙ=ħω(n+½) ✓, Δx·Δp=0.5 ✓ |
| 谐振子 | Ehrenfest 定理 | Norm drift < 1e-8, 经典运动轨迹 |
| 势垒隧穿 | E < V₀ 量子隧穿 | 隧穿概率 39.8%, WKB 预测 13.5% |
| 无限深势阱 | 本征态、反弹 | 能级误差 < 1e-5, 概率守恒 > 98% |
| 双势阱 | 隧穿振荡 | ΔE 劈裂决定振荡周期 |

## 📄 文档

| 文档 | 内容 |
|------|------|
| [PHYSICS.md](docs/PHYSICS.md) | 量子力学公式、数值方法推导、势函数解析性质 |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | 模块架构、数据流、设计原则 |
| [CAPABILITIES.md](docs/CAPABILITIES.md) | 完整功能矩阵、命令参考 |
| [TUTORIAL.md](docs/TUTORIAL.md) | 交互命令教程 + Python API 示例 |
| [CHANGELOG.md](docs/CHANGELOG.md) | 版本历史和变更记录 |

## 🤝 贡献

欢迎 Issue 和 PR。开发前请阅读 [ARCHITECTURE.md](docs/ARCHITECTURE.md) 了解设计原则。

## 📜 License

MIT License — 详见 [LICENSE](LICENSE)。

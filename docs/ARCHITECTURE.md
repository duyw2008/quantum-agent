# Quantum Agent 架构文档

## 系统概述

Quantum Agent 是一个模块化的量子力学计算与可视化平台。核心设计理念是"计算与展示分离"，将数值方法、物理模型和可视化解耦。

## 模块架构

```
                    ┌─────────────┐
                    │   agent.py  │  CLI 交互界面
                    └──────┬──────┘
           ┌───────────────┼───────────────┐
    ┌──────▼──────┐ ┌──────▼──────┐ ┌──────▼──────┐
    │  src/core/  │ │ src/matrix/ │ │  src/viz/   │
    │  核心物理   │ │  矩阵力学   │ │  可视化     │
    └─────────────┘ └─────────────┘ └─────────────┘
```

## 核心模块 (src/core/)

### potentials.py — 势函数
- 抽象基类 `Potential` + 9 种具体势函数
- 工厂函数 `create_potential()` 统一创建接口
- 每个势函数封装了物理参数和解析性质 (如特征长度、Morse 能级公式)

### wave_function.py — 波函数表示
- `Grid` — 空间/动量网格，含 FFT 频率管理
- `WaveFunction` — 复值波函数，提供:
  - 概率密度 ρ(x) = |ψ|² 和相位 φ(x)
  - 归一化、内积、重叠积分
  - 期望值: ⟨x⟩, ⟨p⟩, ⟨x²⟩, ⟨p²⟩, ⟨T⟩, ⟨V⟩
  - 不确定度: Δx, Δp, Δx·Δp
  - 初始化: 高斯波包、平面波、本征态
  - 动量空间表示 ψ̃(k) = F[ψ(x)]

### schrodinger.py — TDSE 求解器
- `TDSE_Solver` 抽象基类
- `SplitStepFourier` (SSFM): 基于 FFT 的谱方法，O(N log N)
  - 二阶对称 Trotter 分解
  - 适合光滑势、长时演化
- `CrankNicolson` (CN): 隐式三对角方法
  - 无条件稳定、二阶精度
  - 使用复 Thomas 算法 O(N)
  - 适合刚性势、保范数要求高
- `EvolutionResult` — 存储时间演化数据，含可观测量时间序列

## 矩阵力学模块 (src/matrix/)

### operators.py — 数值算符
- `NumericOperatorSystem` — 在截断数态基中表示算符
  - â, â† (产生/湮灭), x̂, p̂, N̂
  - 哈密顿量构建 (谐振子 + 一般势)
  - 本征值求解、期望值、不确定度
  - 时间演化 (演化算符、海森堡绘景)
  - 虚时间演化 (基态搜索)
  - 对易子/反对易子、厄米性/幺正性检查

### symbolic.py — 符号算符 (可选 sympy)
- `SymbolicQuantum` — 符号量子力学
  - 算符矩阵的符号表示
  - 能级解析公式
  - 角动量代数、Pauli 自旋矩阵
- `MatrixMechanics` — 统一接口 (数值为主)

## 可视化模块 (src/viz/)

### animate.py — 动画和绘图
- `animate_evolution()` — 时间演化动画 (MP4/GIF)
  - 显示概率密度、势函数、期望值轨迹
  - 可选动量空间视图
- `plot_potential()` — 势函数图
- `plot_wavefunction()` — 波函数快照 (|ψ|², Re, Im, φ)
- `plot_eigenstates()` — 本征态能级+波函数
- 暗色 (GitHub Dark) 和亮色主题

### static.py — 静态可视化
- `plot_energy_levels()` — 能级图
- `plot_matrix_element()` — 算符矩阵热图
- `plot_phase_space()` — 相空间不确定度椭圆

## 数据流

```
User Command → agent.py → 解析参数
                         → 创建 Grid/WaveFunction/Potential
                         → create_solver()
                         → solver.evolve() → EvolutionResult
                         → animate_evolution() / plot_*()
                         → 输出文件 (output/)
```

## 设计原则

1. **不可变性倾向**: WaveFunction, EvolutionResult 通过 .copy() 显式复制
2. **缓存优化**: 动量空间相位因子在 dt 不变时复用
3. **延迟导入**: scipy 仅在需要对角化时导入
4. **兼容性**: numpy 2.x trapz→trapezoid shim
5. **物理单位**: 默认原子单位 (ℏ = m_e = e = 1)

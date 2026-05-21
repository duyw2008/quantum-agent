# Changelog

All notable changes to Quantum Agent will be documented in this file.

## [0.1.0] — 2026-05-21

### Added
- 11 种势函数：InfiniteWell, Harmonic, PotentialBarrier, FiniteWell, DoubleWell, Morse, Coulomb1D, Periodic, StepPotential, ZeroPotential, CustomPotential
- WaveFunction 波函数表示：高斯波包、平面波、本征态初始化，期望值、不确定度计算
- Split-Step Fourier Method (SSFM) TDSE 求解器 — 基于 FFT 的谱方法，O(N log N)
- Crank-Nicolson (CN) TDSE 求解器 — 隐式三对角方法，无条件稳定
- NumericOperatorSystem — 数态基截断算符表示（â, â†, x̂, p̂, N̂）
- MatrixMechanics 统一矩阵力学接口：对易子、本征值、时间演化
- SymbolicQuantum — sympy 符号量子力学（角动量、Pauli 自旋、BCH 公式）
- 可视化模块：动画（MP4/GIF）+ 静态图（势函数、波函数、本征态、能谱、矩阵热图、相空间）
- 暗色/亮色双主题配色方案
- CLI Agent 交互界面：evolve, matrix, eigenstates, plot, animate, demo 命令
- 5 个 Demo：谐振子、无限深势阱、势垒隧穿、矩阵力学、双势阱
- 35 个测试（12+6+6+7+4），全部通过
- 4 份文档：ARCHITECTURE, PHYSICS, CAPABILITIES, TUTORIAL
- numpy 2.x 兼容性 shim（np.trapz → np.trapezoid）

### Verified
- [x̂, p̂] = iħ 对易关系 ✓
- [â, â†] = I 产生湮灭对易 ✓
- Eₙ = ħω(n + ½) 谐振子能谱 ✓（误差 < 1e-8）
- Δx·Δp ≥ 0.5 不确定度原理 ✓
- 量子隧穿：E < V₀ 时概率传输 ✓
- 双势阱隧穿劈裂 ✓
- 范数守恒：norm drift < 1e-8（长时间演化）✓

# Quantum Agent 功能清单

## 势函数 (9+ 种)

| 势函数 | 类型 | 参数 | 解析能级 |
|--------|------|------|----------|
| InfiniteWell | 无限深势阱 | width | ✓ |
| Harmonic | 谐振子 | omega, mass | ✓ |
| PotentialBarrier | 势垒/势阱 | height, width | — |
| FiniteWell | 有限深势阱 | depth, width | 数值 |
| DoubleWell | 双势阱 | separation, depth | 数值 (隧穿劈裂) |
| Morse | Morse 分子势 | depth, alpha, x0 | ✓ |
| Coulomb1D | 软核库仑 | Z, softening | 数值 |
| Periodic | 周期势 | amplitude, k | 能带结构 (Bloch) |
| StepPotential | 阶梯势 | height, x0 | 散射问题 |
| ZeroPotential | 自由粒子 | — | 连续谱 |
| CustomPotential | 自定义 | func, name | — |

## TDSE 求解器

| 方法 | 类型 | 精度 | 稳定性 | 复杂度 |
|------|------|------|--------|--------|
| Split-Step Fourier | 显式/谱方法 | O(Δt³) | 条件稳定 | O(N log N) |
| Crank-Nicolson | 隐式/有限差分 | O(Δt²) | 无条件稳定 | O(N) |

## 矩阵力学

### 算符表示 (数态基)
- 湮灭算符 â
- 产生算符 â†
- 坐标算符 x̂
- 动量算符 p̂
- 数算符 N̂ = â†â
- 谐振子哈密顿量 Ĥ = ħω(N̂ + ½)
- 一般势哈密顿量 Ĥ = p̂²/2m + V(x̂)

### 计算功能
- 对易子 [Â, B̂] 和反对易子 {Â, B̂}
- 本征值/本征态求解
- 期望值 ⟨ψ|Ô|ψ⟩ 和不确定度 ΔO
- 时间演化 |ψ(t)⟩ = e^{-iĤt/ħ}|ψ(0)⟩
- 海森堡绘景 Ô(t) = e^{iĤt/ħ} Ô e^{-iĤt/ħ}
- 虚时间演化 (基态搜索)
- 厄米性/幺正性检查

### 符号力学 (sympy, 可选)
- 算符矩阵的符号表示
- 能级解析公式
- 角动量代数 [J_i, J_j] = iħ ε_ijk J_k
- 自旋-½ Pauli 矩阵
- BCH 公式和 Hadamard 引理

## 可视化

### 静态图
- 势函数 V(x) 图
- 波函数快照 (|ψ|², Re, Im, φ)
- 本征态能级 + 波函数
- 能谱图
- 算符矩阵热图
- 相空间不确定度椭圆

### 动画
- 波函数时间演化 (概率密度 + 期望值)
- 概率密度热图
- 支持 MP4 和 GIF 格式
- 暗色/亮色双主题

## 交互命令

```
evolve <potential> [options]    — 波函数时间演化
matrix <subcommand>             — 矩阵力学操作
eigenstates <potential> [N]     — 本征态计算
plot <type>                     — 绘图
animate <type>                  — 动画
demo <name|all>                 — 运行演示
status                          — 查看当前状态
help                            — 帮助
```

### Demo 案例

| Demo | 说明 | 物理亮点 |
|------|------|----------|
| harmonic_oscillator | 谐振子波包 | Ehrenfest 定理、相干态 |
| infinite_well | 无限深势阱 | 本征态、波包反弹 |
| potential_barrier | 势垒隧穿 | 量子隧穿、WKB 近似 |
| matrix_mechanics | 矩阵力学 | [x̂,p̂]=iħ、能谱、自旋 |
| double_well | 双势阱 | 隧穿振荡、能级劈裂 |
| hydrogen_atom | 氢原子库仑势 | 本征态、能级、径向分布 |
| hydrogen_animation | 库仑势演化 | 波包散射、量子动力学 |

## 测试

35 个测试覆盖：
- 势函数 (12 tests) — 所有势函数类型的正确性
- 波函数 (6 tests) — 初始化、归一化、期望值、本征态
- 求解器 (6 tests) — SSFM、CN、守恒律、一致性
- 矩阵力学 (7 tests) — 算符、对易子、本征值、演化
- 可视化 (4 tests) — 所有绘图函数生成有效图像

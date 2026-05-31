# Changelog

## [2.7.1] — 2026-05-31

### Added
- 知识手册 §0.7: 相干态完整解释
  - 定义、光子统计 / g² / Mandel Q / Wigner
  - â|α⟩ = α|α⟩ 逐行推导
  - Displacement 算符 D̂(α)
  - 过完备基关系 (1/π)∫d²α|α⟩⟨α| = Î
  - 与压缩态全面对比

## [2.7.0] — 2026-05-27

### Added
- `PotentialBuilder` — 链式势函数构造器
  - 11 个基本块: harmonic/barrier/well/gaussian/periodic/delta/step/linear/custom
  - 3 个快捷组合: double_well/tunnel_junction/optical_lattice
  - 代数操作: add/multiply/offset
  - 可视化: plot()/summary()/to_qms()
  - 7 组件 demo: potential_builder_demo.qms (837KB 动画)
- 21 个 `help <function>` 物理公式，LaTeX 引擎渲染

### Changed
- help 命令支持函数物理公式查询 (Unicode + PNG 双输出)
- LaTeX→Unicode 转换增强: 智能上下标、复杂指数保留

## [2.6.0] — 2026-05-26

### Added
- 5 个新 .qms 脚本: free_particle, heisenberg_uncertainty, energy_collapse,
  double_well (双阱隧穿), pimc_demo (路径积分 Monte Carlo)
- `src/qft/path_integral.py` — PathIntegralMC 类 (Phase 4)
  - 欧几里得路径积分离散化 + Metropolis-Hastings 采样
  - ground_state_energy(), wavefunction_density()
- `src/qm/wave.py` — 6 个特色势函数工厂
  - double_well(), periodic_potential(), delta_barrier()
  - finite_well(), harmonic_oscillator_potential(), step_potential()
- `help <function>` 命令 — 18 个物理公式即时查询
  - coherent, g2, wigner, sesolve, pathintegralmc, ...
- `docs/KNOWLEDGE_HANDBOOK.md` — 知识手册 (5卷, 25KB)
  - Vol.0: 为什么需要 Fock 空间
  - Vol.1: 闵氏空间与 Fock 空间的缝合
  - Vol.2: 因果关系在量子力学中的表现
  - Vol.3: 为什么需要 Wigner 函数
  - Vol.4: 退相干的机制

### Changed
- `.qms` 多行语法: 支持 dict/for/if/def/try 等复合语句
- `.qms` 脚本支持 `cd`/`pwd`/`ls` 命令 + Tab 路径补全
- formula 命令 LaTeX→Unicode 转换改为原生终端显示
- `docs/USER_GUIDE.md` 重构: Python .py 与 Agent .qms 双模式分章

## [2.5.0] — 2026-05-24

### Added
- `.qms` 脚本支持 (类 MATLAB .m 文件): 批量执行量子命令
  - `python agent.py --run <script.qms>` — 命令行启动
  - `run <script.qms>` — 交互模式下调用
  - `#` 注释、变量跨行共享、错误继续执行
  - `import` / `print` 等 Python 内置函数可用
- 示例脚本: `scripts/harmonic_oscillator.qms` (5步谐振子分析)
- 示例脚本: `scripts/core_formulas.qms` (10个核心QM公式)

### Changed
- 命令分发重构为 `_dispatch()` 方法, 脚本和交互模式共用
- calc 命名空间扩展: 支持 `import` 语句, `print`/`abs`/`len`/`range` 等内置函数
- Unicode 映射新增: `\sinh`→sinh, `\cosh`→cosh, `\tanh`→tanh

## [2.4.0] — 2026-05-24

### Added
- 公式终端显示: `formula` 命令将 LaTeX 实时转换为 Unicode 数学符号
  - 内置 90+ 映射 (希腊字母、算符、hat、上下标、积分、梯度和二阶导等)
  - `\frac`/`\sqrt`/`\left`/`\right` 智能处理
  - PNG 同步保存到 `output/formulas/`

### Changed
- `formula` 命令: 从 ASCII art (img2txt, 不可读) 改为原生 Unicode 终端显示
- `\hat{H}` 渲染为预组合字符 `Ĥ` (U+0124)，终端兼容性更好
- 公式图 fontsize: 20→30, figsize 固定 10×1.0

## [2.2.0] — 2026-05-23

### Added
- 双缝干涉实验 — 2D TDSE 模拟 (256×128), gamma 校正, inferno 配色
- 量子擦除实验 — 相干 vs 非相干对比, 中途切换模式, 干涉项可视化
- 海森堡不确定性原理动画 — 4 面板 Δx·Δp ≥ ℏ/2
- 位置测量坍缩动画 — 坍缩后 100× 弥散加速
- 动量测量坍缩动画 — 3 面板频率可视化, Δp 宽度标注
- 自由粒子量子弥散动画 — SSFM 验证
- 波函数动力学模块 (`src/qm/wave.py`) — WaveGrid, SSFM, 动画生成

### Changed
- 量子擦除动画增强：单场景中途切换 + 干涉项面板 + 探测器累积
- 双缝动画增强：γ=0.45, inferno, p₀=6, 圆波包, aspect='equal'
- 动量坍缩增强：3 面板, Δp 箭头标注, 频率聚焦显示

## [2.0.0] — 2026-05-22

### Added
- QuTiP 风格量子力学库 (`src/qm/`)
  - `basis.py` — FockBasis: a, a†, x, p, N, parity, displacement
  - `states.py` — fock, coherent, squeezed, thermal_dm, cat
  - `operators.py` — commutator, expect, variance, g2, mandel_q
  - `dynamics.py` — sesolve, mesolve, steadystate
- 可视化模块 (`src/viz/`) — Wigner, Qfunc, 光子分布图
- Agent CLI (`agent.py`) — calc, demo, test, readline
- 文档: MATHEMATICS.md (7章), USER_GUIDE.md (9章)
- 验证: [x̂,p̂]=iħ (6.75×10⁻¹⁶), g²=1.0/2.0, ⟨n⟩=sinh²(r)

### Removed
- 旧模块全部清理 (core, matrix, viz, qubit, qoptics)

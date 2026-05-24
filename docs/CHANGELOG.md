# Changelog

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

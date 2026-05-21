# Quantum Agent 教程

## 安装

```bash
cd /home/duyw/quantum_agent
python3 -m venv /tmp/quantum_venv
source /tmp/quantum_venv/bin/activate
pip install numpy scipy matplotlib
```

## 快速体验

```bash
# 查看所有 demo
python agent.py --list

# 运行所有 demo
python agent.py --demo all

# 运行单个 demo
python agent.py --demo harmonic
```

## 交互式使用

```bash
python agent.py
```

### 示例 1: 谐振子波包

```
⚛ > evolve harmonic,omega=2.0 x0=1.0 p0=5.0 t_max=10
```

这会:
1. 创建 ω=2 的谐振子势
2. 初始化一个在 x₀=1.0, p₀=5.0 的高斯波包
3. 用 SSFM 方法演化 10 个时间单位
4. 生成动画保存到 `output/animations/`

### 示例 2: 势垒隧穿

```
⚛ > evolve barrier,height=10.0,width=0.5 x0=-3.0 p0=4.0 sigma=0.6 t_max=5
```

动能 E = p²/(2m) = 8.0 < V₀ = 10.0，经典粒子无法穿越，但量子粒子可以隧穿！

### 示例 3: 矩阵力学

```
⚛ > matrix report
⚛ > matrix comm x p
⚛ > matrix eigen k=5
```

### 示例 4: 本征态计算

```
⚛ > eigenstates harmonic,omega=2.0 5
```

## Python API 使用

```python
import sys; sys.path.insert(0, '/home/duyw/quantum_agent')
import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution, plot_wavefunction
from src.matrix import MatrixMechanics

# 1. 创建势函数
V = create_potential(PotentialType.HARMONIC, omega=2.0)

# 2. 创建网格和波函数
grid = Grid(-6.0, 6.0, 1024)
wf = WaveFunction(grid)
wf.set_gaussian(x0=1.0, p0=3.0, sigma=0.5)

# 3. 演化
solver = create_solver('ssfm', grid, V)
result = solver.evolve(wf, t_max=5.0, dt=0.005)

# 4. 查看结果
print(f"⟨x⟩(t=5) = {result.expectation_x[-1]:.4f}")
print(f"Energy: {result.energy[0]:.4f} → {result.energy[-1]:.4f}")

# 5. 生成动画
animate_evolution(result, V, save_path='harmonic.mp4', fps=30)

# 6. 矩阵力学
mm = MatrixMechanics(n_basis=50, omega=2.0)
print(mm.report())
```

## NumPy 2.x 兼容

numpy 2.x 移除了 `np.trapz`，代码已包含兼容 shim (自动使用 `np.trapezoid`)。

## 常见问题

**Q: 为什么 [x̂, p̂] ≠ iℏI 在截断基下？**
A: 有限维截断破坏了精确对易关系。在低能子空间 (< N-5) 中，关系近似成立。

**Q: SSFM vs CN 如何选择？**
A: SSFM 适合光滑势、长时演化、谱精度要求高；CN 适合刚性势、守恒性要求高。

**Q: 如何添加自定义势函数？**
A: 使用 `create_potential(PotentialType.CUSTOM, func=lambda x: ..., name="...")`。

**Q: 网格大小如何选择？**
A: 推荐 2 的幂 (512, 1024, 2048) 以优化 FFT 性能。空间范围应至少覆盖波包范围的 3-5 倍。

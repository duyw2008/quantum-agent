# Quantum Agent 教程

## 快速开始

```bash
pip install numpy scipy matplotlib
python agent.py          # 交互模式
python agent.py --demo   # Fock 基演示
python agent.py --test   # 自检
```

## 交互式使用

```
⚛ > calc psi = coherent(20, 2.0)     # 创建相干态
⚛ > calc g2(psi)                      # g²(0) = 1.0
⚛ > calc x, p, W = wigner(psi)        # Wigner 函数
⚛ > calc plot_wigner(x, p, W)         # 绘图

⚛ > calc g = WaveGrid(-30, 30, 512)   # 波函数网格
⚛ > calc psi0 = gaussian_wavepacket(g, 0, 2, 1)
⚛ > calc r = evolve_ssfm(psi0, g, t_max=8)
⚛ > calc animate_wave(r, save_path='wave.gif')
```

## Demo 动画

```bash
python demos/heisenberg_uncertainty.py   # 不确定性原理
python demos/measurement_collapse.py     # 位置坍缩
python demos/momentum_collapse.py        # 动量坍缩
python demos/energy_collapse.py          # 能量坍缩
python demos/double_slit.py              # 双缝干涉
python demos/quantum_eraser.py           # 量子擦除
python demos/free_particle.py            # 自由弥散
```

## Python API

```python
from src.qm import *
import numpy as np

fb = FockBasis(30)
psi = coherent(30, 2.0)
print(f"⟨n⟩={mean_photon(psi, fb):.2f}, g²={g2(psi, fb):.3f}")

# 衰减
H = fb.hamiltonian()
rho0 = coherent_dm(30, 3.0)
t = np.linspace(0, 5, 50)
r = mesolve(H, rho0, t, c_ops=[np.sqrt(0.3)*fb.a], e_ops=[fb.n_op])
```

## 常用工作流

**光子统计对比**:
```python
for state in [coherent(30,2), thermal_dm(30,2), fock(30,4)]:
    print(g2(state))
# → 1.0, 2.0, 0.75
```

**Wigner 函数**:
```python
rho = coherent_dm(20, 1+0.5j)
x, p, W = wigner(rho, N_grid=61)
plot_wigner(x, p, W, save='coh_wigner.png')
```

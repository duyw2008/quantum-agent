#!/usr/bin/env python3
"""氢原子库仑势 — 波函数时间演化动画"""

import sys, os
sys.path.insert(0, '/home/duyw/quantum_agent')

import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution

# ============================================================
# 参数
# ============================================================
Z = 1.0
softening = 0.1
hbar = 1.0
mass = 1.0

# 网格：聚焦在核心区域
x_range = (-12, 12)
n_points = 1024
grid = Grid(x_range[0], x_range[1], n_points)
V = create_potential(PotentialType.COULOMB_1D, Z=Z, softening=softening)

print(f"Potential: {V.label}")
print(f"Grid: {grid}")

# ============================================================
# 初始波包：远离原点，向中心运动
# ============================================================
x0 = -6.0     # 从左侧远处出发
p0 = 2.0      # 向右运动，有足够动能靠近核
sigma = 0.8   # 波包宽度

wf = WaveFunction(grid)
wf.set_gaussian(x0=x0, p0=p0, sigma=sigma)

E_kinetic = p0**2 / (2 * mass)
print(f"\nInitial wavepacket:")
print(f"  x₀ = {x0},  p₀ = {p0}")
print(f"  Kinetic energy: {E_kinetic:.3f} a.u.")
print(f"  ⟨x⟩ = {wf.expectation_x():.4f}, Δx·Δp = {wf.uncertainty_product():.4f}")

# ============================================================
# 时间演化
# ============================================================
t_max = 18.0   # 足够看到波包靠近核、反弹、散射
dt = 0.005
snapshots = 300

print(f"\nEvolving: t_max={t_max}, dt={dt}")
solver = create_solver('ssfm', grid, V, hbar, mass)
result = solver.evolve(wf, t_max, dt, 
                       snapshot_interval=max(1, int(t_max / dt / snapshots)))

print(f"\n✓ Evolution complete: {len(result.times)} snapshots")
print(f"  Initial energy: {result.energy[0]:.4f}")
print(f"  Final energy:   {result.energy[-1]:.4f}")
print(f"  Norm drift:     {abs(result.norm_history[-1] - 1.0):.2e}")

# ============================================================
# 生成动画
# ============================================================
save_dir = '/home/duyw/quantum_agent/output/animations'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'hydrogen_coulomb_evolution.mp4')

print(f"\nGenerating animation → {save_path}")
animate_evolution(result, V, save_path=save_path, fps=30, dark=True,
                  show_momentum=False)
print(f"✓ Animation saved!")
print(f"\nFile: {save_path}")

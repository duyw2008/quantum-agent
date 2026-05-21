#!/usr/bin/env python3
"""
Demo 2: 无限深势阱中的量子力学

展示:
    1. 初始高斯波包在阱壁间的反弹
    2. 本征态和能谱
    3. 波包展宽和 revival 现象
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution, plot_eigenstates, plot_wavefunction

print("=" * 60)
print("  Demo 2: Infinite Square Well")
print("=" * 60)

hbar = 1.0
mass = 1.0

# 阱参数
well_width = 4.0
print(f"\nWell width: a = {well_width}")
print(f"Ground state energy (analytic): E₁ = π²ℏ²/(2ma²) = {np.pi**2 / (2 * well_width**2):.4f}")

# 势函数 (比阱宽稍大一点的区域来容纳指数衰减的尾部)
V = create_potential(PotentialType.INFINITE_WELL, well_width=well_width)
print(f"Potential: {V.label}")

# 网格
x_padding = 0.5
grid = Grid(-well_width/2 - x_padding, well_width/2 + x_padding, 1024)
print(f"Grid: {grid}")

# ============================================================
# Part 1: 本征态计算
# ============================================================
print(f"\n--- Part 1: Eigenstates ---")
n_states = 5
fig, eigvals, eigvecs = plot_eigenstates(grid, V, n_states=n_states,
                                          save_path=os.path.join(
                                              os.environ.get('output_dir',
                                              os.path.join(os.path.dirname(__file__), '..', 'output')),
                                              'figures', 'infinite_well_eigenstates.png'))

# 验证解析解
for n in range(n_states):
    E_analytic = np.pi**2 * (n + 1)**2 / (2 * mass * well_width**2)
    print(f"  E[{n}] = {eigvals[n]:.6f}  (analytic: {E_analytic:.6f}, "
          f"error: {abs(eigvals[n] - E_analytic):.2e})")

# ============================================================
# Part 2: 波包动力学
# ============================================================
print(f"\n--- Part 2: Wavepacket Dynamics ---")

# 井中心的高斯波包，带非零动量
x0 = 0.0
p0 = 4 * np.pi / well_width  # p 略大于第一激发态动量
sigma = well_width / 8

wf = WaveFunction(grid)
wf.set_gaussian(x0=x0, p0=p0, sigma=sigma)

print(f"Initial wavepacket: x₀={x0}, p₀={p0:.4f}, σ={sigma:.4f}")
print(f"  ⟨x⟩ = {wf.expectation_x():.4f}, ⟨p⟩ = {wf.expectation_p():.4f}")
print(f"  Δx·Δp = {wf.uncertainty_product():.4f}")

# 演化
# 特征反弹时间: T = 2a / v = 2ma/p₀
t_bounce = 2 * mass * well_width / abs(p0)
t_max = 3 * t_bounce
dt = min(0.001, t_bounce / 200)
print(f"Bounce time: T = {t_bounce:.4f}")
print(f"Evolving for t_max = {t_max:.4f} ({t_max/t_bounce:.1f} bounces)")

solver = create_solver('ssfm', grid, V, hbar, mass)
result = solver.evolve(wf, t_max, dt, snapshot_interval=max(1, int(t_max/dt/200)))

print(f"\n✓ Evolution complete: {len(result.times)} snapshots")
print(f"  Energy: E₀ = {result.energy[0]:.4f}, E_final = {result.energy[-1]:.4f}")
print(f"  Norm drift: {abs(result.norm_history[-1] - 1.0):.2e}")

# ============================================================
# 输出
# ============================================================
output_dir = os.environ.get('output_dir',
    os.path.join(os.path.dirname(__file__), '..', 'output'))

fig_path = os.path.join(output_dir, 'figures', 'infinite_well_snapshot.png')
plot_wavefunction(wf, V=V, save_path=fig_path,
                  title='Infinite Well: Initial Wavepacket')
print(f"  ✓ Snapshot: {fig_path}")

anim_path = os.path.join(output_dir, 'animations', 'infinite_well.mp4')
animate_evolution(result, V, save_path=anim_path, fps=30, dark=True)
print(f"  ✓ Animation: {anim_path}")

print(f"\n{'='*60}")
print(f"  Demo complete!")
print(f"{'='*60}")

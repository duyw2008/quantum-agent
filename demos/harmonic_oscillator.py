#!/usr/bin/env python3
"""
Demo 1: 谐振子波包演化

展示一个初始高斯波包在一维谐振子势中的量子运动。
验证 Ehrenfest 定理: ⟨x⟩(t) 和 ⟨p⟩(t) 遵循经典运动方程。
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution, plot_wavefunction

print("=" * 60)
print("  Demo 1: Harmonic Oscillator Wavepacket Evolution")
print("=" * 60)

# ============================================================
# 参数设置
# ============================================================
omega = 2.0
hbar = 1.0
mass = 1.0

# 特征长度 a_ho = √(ℏ/mω)
a_ho = np.sqrt(hbar / (mass * omega))
print(f"\nCharacteristic length a_ho = {a_ho:.4f}")
print(f"Classical period T = {2 * np.pi / omega:.4f}")

# 创建势函数
V = create_potential(PotentialType.HARMONIC, omega=omega, mass=mass)
print(f"Potential: {V.label}")

# 创建网格 (以 a_ho 为单位)
x_range = (-6 * a_ho, 6 * a_ho)
grid = Grid(x_range[0], x_range[1], 1024)
print(f"Grid: {grid}")

# 初始波包: 偏离中心的相干态
# 相干态在演化中保持形状不变，质心做经典谐振动
x0 = 2.0 * a_ho   # 偏离平衡位置
p0 = 0.0           # 初始静止
sigma = a_ho / np.sqrt(2)  # 基态宽度
wf = WaveFunction(grid)
wf.set_gaussian(x0=x0, p0=p0, sigma=sigma)

print(f"\nInitial state:")
print(f"  ⟨x⟩ = {wf.expectation_x():.4f}")
print(f"  ⟨p⟩ = {wf.expectation_p():.4f}")
print(f"  Δx·Δp = {wf.uncertainty_product():.4f} (HUP minimum = 0.5)")

# ============================================================
# 时间演化
# ============================================================
t_max = 3 * 2 * np.pi / omega  # 3 个周期
dt = 0.005
print(f"\nEvolving for t_max = {t_max:.4f} ({t_max * omega / (2*np.pi):.1f} periods)")
print(f"dt = {dt}, method = SSFM")

solver = create_solver('ssfm', grid, V, hbar, mass)
result = solver.evolve(wf, t_max, dt, snapshot_interval=10)

print(f"\n✓ Evolution complete: {len(result.times)} snapshots")
print(f"  Energy: E₀ = {result.energy[0]:.4f}, E_final = {result.energy[-1]:.4f}")
print(f"  Energy drift: {abs(result.energy[-1] - result.energy[0]) / abs(result.energy[0]):.2e}")
print(f"  Norm drift: {abs(result.norm_history[-1] - 1.0):.2e}")

# 检查 Ehrenfest 定理: ⟨x⟩(t) ≈ x₀ cos(ωt) + (p₀/mω) sin(ωt)
t = result.times
x_classical = x0 * np.cos(omega * t)  # p₀=0 的情况
x_error = np.abs(result.expectation_x - x_classical).max()
print(f"  Ehrenfest theorem error: max|⟨x⟩ - x_classical| = {x_error:.2e}")

# ============================================================
# 生成动画和图表
# ============================================================
output_dir = os.environ.get('output_dir', os.path.join(os.path.dirname(__file__), '..', 'output'))

print(f"\nGenerating outputs...")

# 静态快照
fig_path = os.path.join(output_dir, 'figures', 'harmonic_snapshot.png')
plot_wavefunction(wf, V=V, save_path=fig_path, title='Harmonic Oscillator: Initial Wavepacket')
print(f"  ✓ Snapshot: {fig_path}")

# 动画
anim_path = os.path.join(output_dir, 'animations', 'harmonic_oscillator.mp4')
animate_evolution(result, V, save_path=anim_path, fps=30, dark=True)
print(f"  ✓ Animation: {anim_path}")

print(f"\n{'='*60}")
print(f"  Demo complete!")
print(f"{'='*60}")

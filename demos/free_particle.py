#!/usr/bin/env python3
"""自由粒子高斯波包量子弥散 — 动画生成"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave

print("Free Particle Gaussian Wavepacket Spreading")
print("=" * 50)

# 网格
grid = WaveGrid(x_min=-15, x_max=15, N=1024)
print(f"Grid: x ∈ [{grid.x[0]:.1f}, {grid.x[-1]:.1f}], N={grid.N}, dx={grid.dx:.4f}")

# 初始波包
x0, p0, sigma = 0.0, 5.0, 1.0
psi0 = gaussian_wavepacket(grid, x0=x0, p0=p0, sigma=sigma)
print(f"Wavepacket: x₀={x0}, p₀={p0}, σ={sigma}")
print(f"  Initial Δx = {sigma:.3f}")

# 理论：自由高斯波包的弥散
# Δx(t) = σ √(1 + (ℏt/2mσ²)²)
hbar, mass = 1.0, 1.0
t_spread = 2 * mass * sigma**2 / hbar  # 特征弥散时间
print(f"  Spreading timescale τ = 2mσ²/ℏ = {t_spread:.1f}")

# 演化 — 看到明显弥散需要跑几个 τ
t_max = 4 * t_spread
dt = 0.01
print(f"  Evolving to t={t_max:.1f} (≈ {t_max/t_spread:.1f} τ), dt={dt}")

result = evolve_ssfm(psi0, grid, dt=dt, t_max=t_max, snapshots=200)

print(f"\nResults:")
print(f"  Snapshots: {len(result['times'])}")
print(f"  Final ⟨x⟩ = {result['x_exp'][-1]:.3f} (expected p₀t/m = {p0*t_max/mass:.1f})")
print(f"  Final Δx  = {result['dx'][-1]:.3f}")
# 理论: Δx(t) = σ √(1 + (t/τ)²)
dx_theory = sigma * np.sqrt(1 + (t_max / t_spread)**2)
print(f"  Theory Δx = {dx_theory:.3f}")
print(f"  Energy: {result['energy'][0]:.4f} → {result['energy'][-1]:.4f}")

# 生成动画
save_path = os.path.join(os.path.dirname(__file__),
                          'output', 'animations', 'free_particle_spreading.mp4')
print(f"\nGenerating animation...")
animate_wave(result, save_path=save_path, fps=25)

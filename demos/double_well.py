#!/usr/bin/env python3
"""
Demo 5: 双势阱 — 量子隧穿振荡

展示粒子在双势阱中的隧穿振荡:
    - 局域在单个阱中的态不是本征态
    - 对称/反对称叠加形成本征态
    - 隧穿劈裂 ΔE = E₁ - E₀ 决定振荡周期
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution, plot_eigenstates, plot_wavefunction

print("=" * 60)
print("  Demo 5: Double Well — Tunneling Oscillation")
print("=" * 60)

hbar = 1.0
mass = 1.0

# 双势阱参数
separation = 3.0   # 两阱间距
depth = 8.0        # 阱深

V = create_potential(PotentialType.DOUBLE_WELL,
                     dw_separation=separation, dw_depth=depth)
print(f"\nDouble Well: separation={separation}, depth={depth}")
print(f"Well centers at x = ±{separation/2}")

# 网格 (覆盖两个阱)
grid = Grid(-4.0, 4.0, 1024)
print(f"Grid: {grid}")

# ============================================================
# Part 1: 本征态和隧穿劈裂
# ============================================================
print(f"\n--- Part 1: Eigenstates and Tunneling Splitting ---")

output_dir = os.environ.get('output_dir',
    os.path.join(os.path.dirname(__file__), '..', 'output'))

fig, eigvals, eigvecs = plot_eigenstates(grid, V, n_states=6,
                                          save_path=os.path.join(
                                              output_dir, 'figures', 'double_well_eigenstates.png'))

print(f"\nLowest eigenvalues:")
for n in range(min(6, len(eigvals))):
    print(f"  E[{n}] = {eigvals[n]:.6f}")

tunneling_splitting = eigvals[1] - eigvals[0]
tunneling_period = 2 * np.pi * hbar / tunneling_splitting
print(f"\nTunneling splitting: ΔE = E₁ - E₀ = {tunneling_splitting:.6f}")
print(f"Tunneling period: T = 2πℏ/ΔE = {tunneling_period:.4f}")

# 基态和第一激发态的波函数对称性
psi0 = eigvecs[:, 0]  # 对称 (gerade)
psi1 = eigvecs[:, 1]  # 反对称 (ungerade)

# 构建局域态: |L⟩ = (|0⟩ + |1⟩)/√2, |R⟩ = (|0⟩ - |1⟩)/√2
psi_left = (psi0 + psi1) / np.sqrt(2)
psi_right = (psi0 - psi1) / np.sqrt(2)

# 检查局域性
a = separation / 2
x = grid.x
left_prob = np.trapz(psi_left**2 * (x < 0).astype(float), x)
right_prob = np.trapz(psi_right**2 * (x > 0).astype(float), x)
print(f"\nLocalized state |L⟩ probability in left well: {left_prob:.4f}")
print(f"Localized state |R⟩ probability in right well: {right_prob:.4f}")

# ============================================================
# Part 2: 隧穿动力学
# ============================================================
print(f"\n--- Part 2: Tunneling Dynamics ---")

# 初始态: 局域在左阱
# 用左阱中心的高斯波包近似
x0_L = -separation / 2
sigma = 0.4  # 比阱宽窄

wf = WaveFunction(grid)
wf.set_gaussian(x0=x0_L, p0=0.0, sigma=sigma)

# 施加势函数投影: 将波函数投影到双阱本征态
# 先检查初始态的能量
# V 是 callable, 我们用数值方法
from scipy.linalg import eigh_tridiagonal
N = grid.n_points
dx = grid.dx
ke_diag = np.ones(N) * hbar**2 / (mass * dx**2)
ke_off = np.ones(N - 1) * (-hbar**2 / (2 * mass * dx**2))
Vx = V(x)
diag = ke_diag + Vx

# 计算更多本征态用于投影
n_eigs = min(20, N - 2)
all_eigvals, all_eigvecs = eigh_tridiagonal(diag, ke_off,
                                             select='i', select_range=(0, n_eigs - 1))

# 投影到低能本征态
psi_init = wf.psi
coeffs = all_eigvecs.T @ psi_init * dx
probs = np.abs(coeffs)**2

print(f"\nInitial state decomposition in eigenbasis:")
for n in range(min(6, len(all_eigvals))):
    print(f"  |⟨n|ψ₀⟩|² = {probs[n]:.4f}  (E_n = {all_eigvals[n]:.4f})")

# 演化
t_max = tunneling_period * 1.5  # 1.5 个隧穿周期
dt = min(0.002, tunneling_period / 200)
print(f"\nEvolving for t_max = {t_max:.4f} ({t_max/tunneling_period:.2f} tunneling periods)")
print(f"dt = {dt}")

solver = create_solver('ssfm', grid, V, hbar, mass)
result = solver.evolve(wf, t_max, dt, snapshot_interval=max(1, int(t_max/dt/200)))

print(f"\n✓ Evolution complete: {len(result.times)} snapshots")

# 跟踪左右阱概率
x = grid.x
left_mask = x < 0
right_mask = x > 0
prob = result.prob_density_snapshots

left_probs = np.array([np.trapz(p[left_mask], x[left_mask]) for p in prob])
right_probs = np.array([np.trapz(p[right_mask], x[right_mask]) for p in prob])

# 找到振荡周期
from scipy.signal import find_peaks
peaks, _ = find_peaks(left_probs, height=0.5)
if len(peaks) >= 2:
    T_observed = result.times[peaks[1]] - result.times[peaks[0]]
    print(f"\nObserved tunneling period: {T_observed:.4f}")
    print(f"Expected period: {tunneling_period:.4f}")
    print(f"Ratio: {T_observed/tunneling_period:.4f}")
else:
    print(f"\n(Not enough peaks detected for period measurement)")

print(f"Final left probability: {left_probs[-1]:.4f}")
print(f"Final right probability: {right_probs[-1]:.4f}")

# ============================================================
# 输出
# ============================================================
fig_path = os.path.join(output_dir, 'figures', 'double_well_snapshot.png')
plot_wavefunction(wf, V=V, save_path=fig_path,
                  title='Double Well: Initial State (localized left)')
print(f"  ✓ Snapshot: {fig_path}")

anim_path = os.path.join(output_dir, 'animations', 'double_well_tunneling.mp4')
animate_evolution(result, V, save_path=anim_path, fps=30, dark=True)
print(f"  ✓ Animation: {anim_path}")

print(f"\n{'='*60}")
print(f"  Demo complete!")
print(f"  Double well tunneling demonstrated: the particle oscillates")
print(f"  between the two wells with period T = 2πℏ/ΔE.")
print(f"{'='*60}")

#!/usr/bin/env python3
"""
自由粒子 — 位置测量坍缩与后续量子弥散

物理过程:
  1. 宽高斯波包自由演化 (Δx 大 → 弥散缓慢)
  2. 在 t_measure 进行位置测量 → 坍缩为窄高斯 (Δx 小)
  3. 坍缩后快速弥散 (Δx 小 → Δp 大 → 弥散加速)

Heisenberg 不确定性原理的戏剧性展示。
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm

print("=" * 60)
print("  Quantum Measurement Collapse & Spreading")
print("=" * 60)

# 网格
grid = WaveGrid(x_min=-40, x_max=40, N=1024)
x = grid.x
dx = grid.dx

# ============================================================
# Part 1: 测量前 — 宽波包自由演化
# ============================================================
sigma_wide = 3.0           # 宽波包，动量不确定度小
p0 = 2.0
t_measure = 4.0            # 测量时刻

psi0 = gaussian_wavepacket(grid, x0=-8.0, p0=p0, sigma=sigma_wide)
psi0 /= np.sqrt(np.trapezoid(np.abs(psi0)**2, x))

print(f"\nPre-measurement wavepacket:")
print(f"  σ = {sigma_wide}, p₀ = {p0}")
dx_init = np.sqrt(np.trapezoid(x**2 * np.abs(psi0)**2, x) -
                  np.trapezoid(x * np.abs(psi0)**2, x)**2)
print(f"  Δx(0) = {dx_init:.3f}")

# 演化到测量时刻
result_pre = evolve_ssfm(psi0, grid, dt=0.005, t_max=t_measure, snapshots=80)

psi_before = result_pre['psi'][-1]
x_mean_before = np.trapezoid(x * np.abs(psi_before)**2, x)
dx_before = np.sqrt(np.trapezoid(x**2 * np.abs(psi_before)**2, x) - x_mean_before**2)
print(f"\n  At t={t_measure} (before measurement):")
print(f"  ⟨x⟩ = {x_mean_before:.3f}, Δx = {dx_before:.3f}")

# ============================================================
# Part 2: 位置测量 — 坍缩
# ============================================================
# 模拟测量：以当前概率密度为权重随机选择测量结果
prob = np.abs(psi_before)**2
prob /= prob.sum()

# 为可复现性固定随机种子
np.random.seed(42)
measured_x = np.random.choice(x, p=prob)
print(f"\n  ⚡ Measurement result: x = {measured_x:.3f}")

# 坍缩：用窄高斯近似测量后的态
sigma_narrow = 0.3  # 位置测量精度
psi_after = np.exp(-(x - measured_x)**2 / (2 * sigma_narrow**2))
psi_after = psi_after.astype(complex)
psi_after /= np.sqrt(np.trapezoid(np.abs(psi_after)**2, x))

dx_after = np.sqrt(np.trapezoid(x**2 * np.abs(psi_after)**2, x) - measured_x**2)
print(f"  Collapsed wavefunction: Δx = {dx_after:.3f}")

# 坍缩前后的动量不确定度
def momentum_uncertainty(psi, grid):
    psi_k = np.fft.fft(psi)
    prob_k = np.abs(psi_k)**2
    k = grid.k
    p_mean = np.trapezoid(k * prob_k, k) / np.trapezoid(prob_k, k)
    p2_mean = np.trapezoid(k**2 * prob_k, k) / np.trapezoid(prob_k, k)
    return np.sqrt(p2_mean - p_mean**2), p_mean

dp_before, p_mean_before = momentum_uncertainty(psi_before, grid)
dp_after, p_mean_after = momentum_uncertainty(psi_after, grid)
print(f"  Δp before: {dp_before:.3f}, after: {dp_after:.3f}")
print(f"  Δx·Δp before: {dx_before * dp_before:.3f}  (≥ 0.5)")
print(f"  Δx·Δp after:  {dx_after * dp_after:.3f}  (≥ 0.5)")

# ============================================================
# Part 3: 坍缩后 — 快速弥散
# ============================================================
t_after = 5.0
# 特征弥散时间: τ = 2mσ²/ℏ
tau_before = 2 * sigma_wide**2
tau_after = 2 * sigma_narrow**2
print(f"\n  Spreading timescale: τ_before = {tau_before:.1f}, τ_after = {tau_after:.2f}")
print(f"  Ratio: τ_before/τ_after = {tau_before/tau_after:.0f}× faster after collapse!")

result_after = evolve_ssfm(psi_after, grid, dt=0.002, t_max=t_after, snapshots=100)

psi_final = result_after['psi'][-1]
x_final = np.trapezoid(x * np.abs(psi_final)**2, x)
dx_final = np.sqrt(np.trapezoid(x**2 * np.abs(psi_final)**2, x) - x_final**2)
print(f"\n  At t={t_measure + t_after} (final):")
print(f"  ⟨x⟩ = {x_final:.3f}, Δx = {dx_final:.3f}")

# ============================================================
# Part 4: 生成动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 合并时间线
all_times_pre = np.array(result_pre['times'])
all_probs_pre = np.array(result_pre['prob'])
all_times_after = np.array(result_after['times']) + t_measure
all_probs_after = np.array(result_after['prob'])

all_times = np.concatenate([all_times_pre, all_times_after])
all_probs = np.concatenate([all_probs_pre, all_probs_after])
measure_idx = len(all_times_pre) - 1

fig, (ax_top, ax_main) = plt.subplots(2, 1, figsize=(14, 8),
                                       gridspec_kw={'height_ratios': [1, 3]},
                                       facecolor='#0d1117')

# Top: Δx vs t
ax_top.set_facecolor('#0d1117')
times_dx_pre = all_times_pre
times_dx_after = all_times_after
dx_vals_pre = [np.sqrt(np.trapezoid(x**2 * p, x) - np.trapezoid(x * p, x)**2)
               for p in all_probs_pre]
dx_vals_after = [np.sqrt(np.trapezoid(x**2 * p, x) - np.trapezoid(x * p, x)**2)
                 for p in all_probs_after]

ax_top.plot(times_dx_pre, dx_vals_pre, color='#79c0ff', linewidth=2, label=r'$\Delta x$')
ax_top.plot(times_dx_after, dx_vals_after, color='#ff7b72', linewidth=2)
ax_top.axvline(t_measure, color='#f0883e', linewidth=2, linestyle='--', alpha=0.8)
ax_top.text(t_measure + 0.1, max(dx_vals_pre) * 0.9, 'MEASUREMENT',
            color='#f0883e', fontsize=10, fontweight='bold')
ax_top.set_ylabel(r'$\Delta x$', color='#e6edf3', fontsize=12)
ax_top.tick_params(colors='#e6edf3')
ax_top.grid(True, alpha=0.15, color='#30363d')
ax_top.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3')

# 时间指示线
t_line = ax_top.axvline(0, color='#e6edf3', linewidth=1, alpha=0.4)

# Main: |ψ|²
ax_main.set_facecolor('#0d1117')
fill = ax_main.fill_between(x, 0, all_probs[0], alpha=0.5, color='#79c0ff')
line, = ax_main.plot(x, all_probs[0], color='#79c0ff', linewidth=1.5)
ax_main.set_xlabel('x', color='#e6edf3', fontsize=12)
ax_main.set_ylabel(r'$|\psi(x,t)|^2$', color='#e6edf3', fontsize=12)
ax_main.tick_params(colors='#e6edf3')
ax_main.grid(True, alpha=0.15, color='#30363d')
title = ax_main.set_title('', color='#e6edf3', fontsize=14)

# 测量位置标记
meas_line = ax_main.axvline(measured_x, color='#f0883e', linewidth=1.5,
                              linestyle='--', alpha=0)

ax_main.set_xlim(x[0], x[-1])
all_max = max(p.max() for p in all_probs)
ax_main.set_ylim(0, all_max * 1.2)

for ax in [ax_top, ax_main]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

fill_ref = [fill]

def update(i):
    fill_ref[0].remove()
    is_after = i > measure_idx
    color = '#ff7b72' if is_after else '#79c0ff'
    fill_ref[0] = ax_main.fill_between(x, 0, all_probs[i], alpha=0.5, color=color)
    line.set_ydata(all_probs[i])
    line.set_color(color)
    t_line.set_xdata([all_times[i], all_times[i]])

    # 测量线
    if is_after:
        meas_line.set_alpha(0.6)
    else:
        meas_line.set_alpha(0)

    t = all_times[i]
    dx_val = dx_vals_pre[i] if i <= measure_idx else dx_vals_after[i - measure_idx - 1]
    phase = "FREE" if i <= measure_idx else "COLLAPSED"
    title.set_text(f't = {t:.2f}  |  Δx = {dx_val:.3f}  |  Phase: {phase}')

    return [line, t_line, title]

total_frames = len(all_times)

# 在测量帧暂停 (重复几帧)
extended_indices = []
for i in range(total_frames):
    extended_indices.append(i)
    if i == measure_idx:
        for _ in range(15):  # 暂停
            extended_indices.append(i)

ani = animation.FuncAnimation(fig, update, frames=extended_indices,
                               interval=40, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'measurement_collapse.mp4')

try:
    writer = animation.FFMpegWriter(fps=25, bitrate=2000)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    writer = animation.PillowWriter(fps=25)
    ani.save(save_path, writer=writer, dpi=150)

plt.close(fig)
print(f"\n{'='*60}")
print(f"  Animation saved: {save_path}")
print(f"{'='*60}")

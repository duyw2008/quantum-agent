#!/usr/bin/env python3
"""
自由粒子 — 动量测量坍缩

物理过程:
  1. 窄位置波包 (Δx 小, Δp 大) 自由演化
  2. 动量测量 → 坍缩为窄动量分布 (Δp 小)
  3. 位置空间剧烈展宽 (Δx 大 — 几乎变成平面波)

与位置测量的对比:
  位置测量:  Δx 骤降 → 快速弥散
  动量测量:  Δp 骤降 → 位置急剧展宽, 弥散反而变慢
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm

print("=" * 60)
print("  Momentum Measurement Collapse")
print("=" * 60)

# 大网格 (动量测量后位置展宽很大)
grid = WaveGrid(x_min=-80, x_max=80, N=2048)
x = grid.x
k = grid.k
dx = grid.dx

# ============================================================
# Part 1: 测量前 — 窄波包自由演化
# ============================================================
sigma_x = 0.5            # 窄位置 → 宽动量分布
p0 = 3.0
t_measure = 3.0

psi0 = gaussian_wavepacket(grid, x0=-5.0, p0=p0, sigma=sigma_x)
psi0 /= np.sqrt(np.trapezoid(np.abs(psi0)**2, x))

dx_init = np.sqrt(np.trapezoid(x**2 * np.abs(psi0)**2, x) -
                  np.trapezoid(x * np.abs(psi0)**2, x)**2)
print(f"\nPre-measurement wavepacket (position space):")
print(f"  σ_x = {sigma_x}, p₀ = {p0}")
print(f"  Δx(0) = {dx_init:.3f}")

result_pre = evolve_ssfm(psi0, grid, dt=0.003, t_max=t_measure, snapshots=80)

psi_before = result_pre['psi'][-1]
x_mean_before = np.trapezoid(x * np.abs(psi_before)**2, x)
dx_before = np.sqrt(np.trapezoid(x**2 * np.abs(psi_before)**2, x) - x_mean_before**2)

# 动量空间分析
psi_k_before = np.fft.fft(psi_before)
prob_k_before = np.abs(psi_k_before)**2
p_mean_before = np.trapezoid(k * prob_k_before, k) / np.trapezoid(prob_k_before, k)
dp_before = np.sqrt(np.trapezoid(k**2 * prob_k_before, k) /
                    np.trapezoid(prob_k_before, k) - p_mean_before**2)

print(f"\n  At t={t_measure} (before measurement):")
print(f"  Position:  ⟨x⟩={x_mean_before:.3f}, Δx={dx_before:.3f}")
print(f"  Momentum:  ⟨p⟩={p_mean_before:.3f}, Δp={dp_before:.3f}")
print(f"  Δx·Δp = {dx_before * dp_before:.3f}")

# ============================================================
# Part 2: 动量测量 — 坍缩在动量空间
# ============================================================
np.random.seed(42)
measured_p = np.random.choice(k, p=prob_k_before / prob_k_before.sum())
print(f"\n  ⚡ Momentum measurement: p = {measured_p:.3f}")

# 坍缩: 动量空间窄高斯
sigma_p = 0.3
psi_k_after = np.exp(-(k - measured_p)**2 / (2 * sigma_p**2))
psi_k_after = psi_k_after.astype(complex)
# 归一化
psi_k_after /= np.sqrt(np.trapezoid(np.abs(psi_k_after)**2, k))

# 回到坐标空间
psi_after = np.fft.ifft(psi_k_after)
psi_after /= np.sqrt(np.trapezoid(np.abs(psi_after)**2, x))

x_mean_after = np.trapezoid(x * np.abs(psi_after)**2, x)
dx_after = np.sqrt(np.trapezoid(x**2 * np.abs(psi_after)**2, x) - x_mean_after**2)
p_mean_after = np.trapezoid(k * np.abs(psi_k_after)**2, k) / np.trapezoid(np.abs(psi_k_after)**2, k)
dp_after = np.sqrt(np.trapezoid(k**2 * np.abs(psi_k_after)**2, k) /
                   np.trapezoid(np.abs(psi_k_after)**2, k) - p_mean_after**2)

print(f"\n  After momentum collapse:")
print(f"  Position:  ⟨x⟩={x_mean_after:.3f}, Δx={dx_after:.3f}")
print(f"  Momentum:  ⟨p⟩={p_mean_after:.3f}, Δp={dp_after:.3f}")
print(f"  Δx·Δp = {dx_after * dp_after:.3f}")

print(f"\n  Δp: {dp_before:.3f} → {dp_after:.3f}  (collapsed {dp_before/dp_after:.0f}×)")
print(f"  Δx: {dx_before:.3f} → {dx_after:.3f}  (expanded {dx_after/dx_before:.0f}×)")

# ============================================================
# Part 3: 坍缩后演化
# ============================================================
t_after = 4.0
tau = 2 * (1/(2*sigma_p))**2  # 特征弥散时间 (从 Δx ≈ 1/Δp 估算)
print(f"\n  Spreading timescale after collapse: τ ≈ {tau:.1f}")

result_after = evolve_ssfm(psi_after, grid, dt=0.005, t_max=t_after, snapshots=100)

psi_final = result_after['psi'][-1]
x_final = np.trapezoid(x * np.abs(psi_final)**2, x)
dx_final = np.sqrt(np.trapezoid(x**2 * np.abs(psi_final)**2, x) - x_final**2)
print(f"\n  At t={t_measure + t_after} (final):")
print(f"  ⟨x⟩={x_final:.3f}, Δx={dx_final:.3f}")

# ============================================================
# Part 4: 双面板动画 — 位置空间 + 动量空间
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 合并时间线
all_times_pre = np.array(result_pre['times'])
all_times_after = np.array(result_after['times']) + t_measure

# 计算所有帧的动量空间分布
def compute_momentum_probs(psi_list):
    probs = []
    for psi in psi_list:
        pk = np.abs(np.fft.fft(psi))**2
        probs.append(pk)
    return probs

prob_k_pre = compute_momentum_probs(result_pre['psi'])
prob_k_after = compute_momentum_probs(result_after['psi'])

all_times = np.concatenate([all_times_pre, all_times_after])
all_probs_x = np.concatenate([np.array(result_pre['prob']),
                               np.array(result_after['prob'])])
all_probs_k = np.concatenate([np.array(prob_k_pre), np.array(prob_k_after)])
measure_idx = len(all_times_pre) - 1

# 图形
fig, (ax_x, ax_k) = plt.subplots(2, 1, figsize=(14, 8),
                                   facecolor='#0d1117')

# === 上: 位置空间 ===
ax_x.set_facecolor('#0d1117')
fill_x = ax_x.fill_between(x, 0, all_probs_x[0], alpha=0.5, color='#58a6ff')
line_x, = ax_x.plot(x, all_probs_x[0], color='#58a6ff', linewidth=1.2)
ax_x.set_ylabel(r'$|\psi(x)|^2$', color='#e6edf3', fontsize=12)
ax_x.set_title('', color='#e6edf3', fontsize=14)
ax_x.tick_params(colors='#e6edf3')
ax_x.grid(True, alpha=0.1, color='#30363d')
ax_x.set_xlim(x[0], x[-1])

# === 下: 动量空间 ===
ax_k.set_facecolor('#0d1117')
fill_k = ax_k.fill_between(k, 0, all_probs_k[0], alpha=0.5, color='#d2a8ff')
line_k, = ax_k.plot(k, all_probs_k[0], color='#d2a8ff', linewidth=1.2)
ax_k.set_xlabel('p (momentum)', color='#e6edf3', fontsize=12)
ax_k.set_ylabel(r'$|\tilde{\psi}(p)|^2$', color='#e6edf3', fontsize=12)
ax_k.tick_params(colors='#e6edf3')
ax_k.grid(True, alpha=0.1, color='#30363d')

# 动量测量线
meas_line_k = ax_k.axvline(measured_p, color='#f0883e', linewidth=1.5,
                            linestyle='--', alpha=0)
meas_label_k = ax_k.text(measured_p + 0.5, 0, '', color='#f0883e',
                          fontsize=9, fontweight='bold', va='bottom')

# 统一 y 轴范围
all_x_max = max(p.max() for p in all_probs_x)
all_k_max = max(p.max() for p in all_probs_k)
ax_x.set_ylim(0, all_x_max * 1.2)
ax_k.set_ylim(0, all_k_max * 1.2)
ax_k.set_xlim(k[0], k[-1])

for ax in [ax_x, ax_k]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

fill_x_ref = [fill_x]
fill_k_ref = [fill_k]

plt.tight_layout()

def update(i):
    fill_x_ref[0].remove()
    fill_k_ref[0].remove()

    is_after = i > measure_idx
    color_x = '#ff7b72' if is_after else '#58a6ff'
    color_k = '#ff7b72' if is_after else '#d2a8ff'

    fill_x_ref[0] = ax_x.fill_between(x, 0, all_probs_x[i], alpha=0.5, color=color_x)
    line_x.set_ydata(all_probs_x[i])
    line_x.set_color(color_x)

    fill_k_ref[0] = ax_k.fill_between(k, 0, all_probs_k[i], alpha=0.5, color=color_k)
    line_k.set_ydata(all_probs_k[i])
    line_k.set_color(color_k)

    if is_after:
        meas_line_k.set_alpha(0.6)
        meas_label_k.set_text(f'p={measured_p:.2f}')
    else:
        meas_line_k.set_alpha(0)
        meas_label_k.set_text('')

    t = all_times[i]
    phase = "MOMENTUM COLLAPSED" if is_after else "FREE EVOLUTION"
    dx_val = np.sqrt(np.trapezoid(x**2 * all_probs_x[i], x) -
                     np.trapezoid(x * all_probs_x[i], x)**2)
    dp_val = np.sqrt(np.trapezoid(k**2 * all_probs_k[i], k) /
                     np.trapezoid(all_probs_k[i], k) -
                     (np.trapezoid(k * all_probs_k[i], k) /
                      np.trapezoid(all_probs_k[i], k))**2)
    ax_x.set_title(f't={t:.2f}  |  Δx={dx_val:.2f}  |  Δp={dp_val:.2f}  |  {phase}',
                   color='#e6edf3', fontsize=14)

    return [line_x, line_k]

total_frames = len(all_times)
extended = []
for i in range(total_frames):
    extended.append(i)
    if i == measure_idx:
        for _ in range(15):
            extended.append(i)

ani = animation.FuncAnimation(fig, update, frames=extended,
                               interval=40, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'momentum_collapse.mp4')

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

#!/usr/bin/env python3
"""
自由粒子 — 能量测量坍缩

能量本征态 ≠ 动量本征态！
  动量测量 → 坍缩为行波 e^{ikx} (有方向)
  能量测量 → 坍缩为驻波 cos(kx) ∝ |k⟩+|-k⟩ (无方向, E=k²/2m)

过程:
  1. 高斯波包 (行波, 有确定动量方向) 自由演化
  2. 能量测量 → 随机测得 E
  3. 坍缩为驻波 cos(kx), k=√(2mE)/ℏ
  4. ⟨p⟩=0! 波包不再移动, 但能量确定
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm

print("=" * 60)
print("  Energy Measurement Collapse")
print("=" * 60)

# ============================================================
# 参数
# ============================================================
grid = WaveGrid(x_min=-60, x_max=60, N=2048)
x, k_grid = grid.x, grid.k
hbar, mass = 1.0, 1.0

sigma_x, p0, t_measure = 1.0, 4.0, 4.0
psi0 = gaussian_wavepacket(grid, x0=-8.0, p0=p0, sigma=sigma_x)
psi0 /= np.sqrt(np.trapezoid(np.abs(psi0)**2, x))

print(f"Wavepacket: σ_x={sigma_x}, p₀={p0}")
print(f"E_kinetic ≈ p₀²/2m = {p0**2/2:.1f}")

# 演化到测量时刻
result_pre = evolve_ssfm(psi0, grid, dt=0.005, t_max=t_measure, snapshots=80)
psi_before = result_pre['psi'][-1]

# ============================================================
# 能量分布与测量
# ============================================================
psi_k_before = np.fft.fft(psi_before)
prob_k = np.abs(psi_k_before)**2

# 能量 E = k²/2 → 能量分布 P(E) = P(k)dk/dE + P(-k)dk/dE
# dk/dE = 1/√(2E) = 1/|k|
E_grid = k_grid**2 / 2
# 只对 k>0 采样 (E 简并)
k_pos = k_grid[k_grid > 0]; E_pos = k_pos**2 / 2
prob_E = np.zeros_like(E_pos)
for i, (kp, Ep) in enumerate(zip(k_pos, E_pos)):
    idx_p = np.argmin(np.abs(k_grid - kp))
    idx_m = np.argmin(np.abs(k_grid + kp))
    prob_E[i] = prob_k[idx_p] + prob_k[idx_m]

prob_E /= prob_E.sum()

np.random.seed(12345)
measured_E = np.random.choice(E_pos, p=prob_E)
measured_k = np.sqrt(2 * mass * measured_E) / hbar
print(f"\n⚡ Energy measurement: E = {measured_E:.3f}  →  k = ±{measured_k:.3f}")
print(f"   Expected E ≈ p₀²/2m = {p0**2/2:.1f}")

# ============================================================
# 坍缩为驻波
# ============================================================
# |E⟩ ∝ cos(kx) with k = √(2mE)/ℏ
psi_after = np.cos(measured_k * x).astype(complex)
# 加高斯包络 (有限空间局域化)
envelope = np.exp(-x**2 / (2 * (grid.x[-1]/4)**2))
psi_after *= envelope
psi_after /= np.sqrt(np.trapezoid(np.abs(psi_after)**2, x))

# 动量空间: 两个峰在 ±k
psi_k_after = np.fft.fft(psi_after)

dx_before = np.sqrt(np.trapezoid(x**2 * np.abs(psi_before)**2, x) -
                    np.trapezoid(x * np.abs(psi_before)**2, x)**2)
dx_after = np.sqrt(np.trapezoid(x**2 * np.abs(psi_after)**2, x) -
                   np.trapezoid(x * np.abs(psi_after)**2, x)**2)

# 动量/能量统计
p_mean_before = np.trapezoid(k_grid * prob_k, k_grid) / np.trapezoid(prob_k, k_grid)
p_mean_after = np.trapezoid(k_grid * np.abs(psi_k_after)**2, k_grid) / \
               np.trapezoid(np.abs(psi_k_after)**2, k_grid)

print(f"\nBefore: ⟨x⟩≈{np.trapezoid(x*np.abs(psi_before)**2, x):.1f}, "
      f"Δx={dx_before:.2f}, ⟨p⟩={p_mean_before:.2f}")
print(f"After:  ⟨x⟩≈{np.trapezoid(x*np.abs(psi_after)**2, x):.1f}, "
      f"Δx={dx_after:.1f}, ⟨p⟩={p_mean_after:.3f}")
print(f"⟨p⟩ after collapse ≈ 0 — standing wave doesn't propagate!")

# ============================================================
# 坍缩后演化 (驻波在自由空间不传播但会缓慢弥散)
# ============================================================
t_after = 6.0
result_after = evolve_ssfm(psi_after, grid, dt=0.005, t_max=t_after, snapshots=100)

# ============================================================
# 动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def to_momentum(psi):
    pk = np.abs(np.fft.fft(psi))**2
    return pk

prob_k_pre = [to_momentum(p) for p in result_pre['psi']]
prob_k_after = [to_momentum(p) for p in result_after['psi']]

all_times_pre = np.array(result_pre['times'])
all_times_after = np.array(result_after['times']) + t_measure
all_times = np.concatenate([all_times_pre, all_times_after])
all_probs_x = np.concatenate([np.array(result_pre['prob']),
                               np.array(result_after['prob'])])
all_probs_k = np.concatenate([np.array(prob_k_pre), np.array(prob_k_after)])
measure_idx = len(all_times_pre) - 1

fig = plt.figure(figsize=(16, 9), facecolor='#0d1117')
gs = fig.add_gridspec(2, 2, height_ratios=[2.5, 1.5], hspace=0.4, wspace=0.3)

ax_x = fig.add_subplot(gs[0, 0], facecolor='#0d1117')      # 位置空间
ax_k = fig.add_subplot(gs[0, 1], facecolor='#0d1117')      # 动量/能量空间
ax_t  = fig.add_subplot(gs[1, :], facecolor='#0d1117')      # ⟨x⟩, ⟨p⟩ 时间线

# --- 位置空间 ---
fill_x = ax_x.fill_between(x, 0, all_probs_x[0], alpha=0.5, color='#58a6ff')
line_x, = ax_x.plot(x, all_probs_x[0], color='#58a6ff', linewidth=1.5)
# 驻波节点标记
node_lines = [ax_x.axvline((n+0.5)*np.pi/measured_k, color='#f0883e',
                           linewidth=0.5, linestyle=':', alpha=0)
              for n in range(-8, 8)]
ax_x.set_xlabel('x', color='#e6edf3', fontsize=11)
ax_x.set_ylabel(r'$|\psi(x)|^2$', color='#e6edf3', fontsize=11)
ax_x.set_title('Position Space', color='#e6edf3', fontsize=13, fontweight='bold')
ax_x.tick_params(colors='#e6edf3'); ax_x.grid(True, alpha=0.1, color='#30363d')

# --- 动量/能量空间 ---
# 双 x 轴: 下=动量 k, 上=能量 E=k²/2
fill_k = ax_k.fill_between(k_grid, 0, all_probs_k[0], alpha=0.5, color='#d2a8ff')
line_k, = ax_k.plot(k_grid, all_probs_k[0], color='#d2a8ff', linewidth=1.5)
# ±k 标记线
k_line_p = ax_k.axvline(measured_k, color='#f0883e', linewidth=2, linestyle='--', alpha=0)
k_line_m = ax_k.axvline(-measured_k, color='#f0883e', linewidth=2, linestyle='--', alpha=0)
k_label = ax_k.text(measured_k+0.3, 0, '', color='#f0883e', fontsize=10, fontweight='bold', va='bottom')
ax_k.set_xlabel('k (momentum)', color='#e6edf3', fontsize=11)
ax_k.set_ylabel(r'$|\tilde{\psi}(k)|^2$', color='#e6edf3', fontsize=11)
ax_k.set_title('Momentum / Energy Space', color='#e6edf3', fontsize=13, fontweight='bold')
ax_k.tick_params(colors='#e6edf3'); ax_k.grid(True, alpha=0.1, color='#30363d')
ax_k.set_xlim(-10, 10)
# 能量轴 (顶部)
ax_e = ax_k.twiny()
ax_e.set_xlim(ax_k.get_xlim())
E_ticks = [0, 10, 20, 30, 40, 50]
ax_e.set_xticks([np.sqrt(2*e) for e in E_ticks])
ax_e.set_xticklabels([str(e) for e in E_ticks])
ax_e.set_xlabel('E = k²/2  (energy)', color='#d2a8ff', fontsize=9)
ax_e.tick_params(colors='#d2a8ff', labelsize=8)

# --- 时间演化 ---
# 计算每帧 ⟨x⟩, ⟨p⟩
x_exp_vals, p_exp_vals = [], []
for i in range(len(all_times)):
    px, pk = all_probs_x[i], all_probs_k[i]
    mx = np.trapezoid(x*px, x)
    x_exp_vals.append(mx)
    mp = np.trapezoid(k_grid*pk, k_grid) / np.trapezoid(pk, k_grid)
    p_exp_vals.append(mp)

x_exp_vals = np.array(x_exp_vals); p_exp_vals = np.array(p_exp_vals)

ax_t.plot(all_times_pre, x_exp_vals[:measure_idx+1], color='#58a6ff', linewidth=2, label=r'$\langle x \rangle$')
ax_t.plot(all_times_pre, p_exp_vals[:measure_idx+1], color='#d2a8ff', linewidth=2, linestyle='--', label=r'$\langle p \rangle$')
ax_t.plot(all_times_after, x_exp_vals[measure_idx+1:], color='#ff7b72', linewidth=2)
ax_t.plot(all_times_after, p_exp_vals[measure_idx+1:], color='#ff7b72', linewidth=2, linestyle='--')
ax_t.axvline(t_measure, color='#f0883e', linewidth=2.5, linestyle='--', alpha=0.8)
ax_t.text(t_measure+0.15, max(x_exp_vals)*0.85, 'ENERGY\nMEASUREMENT',
          color='#f0883e', fontsize=10, fontweight='bold')
# ⟨p⟩→0 标注
ax_t.annotate(r'$\langle p \rangle \rightarrow 0$', xy=(t_measure+0.5, 0.2),
              xytext=(t_measure+2, 1.5), arrowprops=dict(arrowstyle='->', color='#f0883e'),
              color='#f0883e', fontsize=11, fontweight='bold')
ax_t.annotate('STATIONARY\n(standing wave)', xy=(t_measure+2, 0),
              xytext=(t_measure+3, -1), arrowprops=dict(arrowstyle='->', color='#ff7b72'),
              color='#ff7b72', fontsize=10)
t_cur = ax_t.axvline(0, color='#e6edf3', linewidth=1, alpha=0.3)
ax_t.set_xlabel('Time', color='#e6edf3', fontsize=11)
ax_t.set_ylabel('Expectation', color='#e6edf3', fontsize=11)
ax_t.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3', loc='upper right')
ax_t.tick_params(colors='#e6edf3'); ax_t.grid(True, alpha=0.1, color='#30363d')

fig_title = fig.suptitle('', color='#e6edf3', fontsize=15, fontweight='bold', y=0.99)

for ax in [ax_x, ax_k, ax_t]:
    for spine in ax.spines.values(): spine.set_color('#30363d')

fill_x_r = [fill_x]; fill_k_r = [fill_k]

def update(i):
    fill_x_r[0].remove(); fill_k_r[0].remove()
    is_after = i > measure_idx
    cx = '#ff7b72' if is_after else '#58a6ff'
    ck = '#ff7b72' if is_after else '#d2a8ff'

    fill_x_r[0] = ax_x.fill_between(x, 0, all_probs_x[i], alpha=0.5, color=cx)
    line_x.set_ydata(all_probs_x[i]); line_x.set_color(cx)

    fill_k_r[0] = ax_k.fill_between(k_grid, 0, all_probs_k[i], alpha=0.5, color=ck)
    line_k.set_ydata(all_probs_k[i]); line_k.set_color(ck)

    # ±k 线 + 节点
    a = 0.7 if is_after else 0
    k_line_p.set_alpha(a); k_line_m.set_alpha(a)
    k_label.set_text(f'±k={measured_k:.2f}' if is_after else '')
    k_label.set_y(max(all_probs_k[i])*0.85 if is_after else 0)
    for nl in node_lines:
        nl.set_alpha(0.3 if is_after else 0)

    t_cur.set_xdata([all_times[i], all_times[i]])
    t = all_times[i]
    phase = '◉ ENERGY EIGENSTATE (standing wave)' if is_after else '○ TRAVELING WAVE'
    fig_title.set_text(
        f't = {t:.2f}    |    ⟨x⟩ = {x_exp_vals[i]:.2f}    |    '
        f'⟨p⟩ = {p_exp_vals[i]:.3f}    |    E = {measured_E:.3f}    |    {phase}')

    return [line_x, line_k, t_cur]

frames_list = []
for i in range(len(all_times)):
    frames_list.append(i)
    if i == measure_idx:
        for _ in range(25): frames_list.append(i)

ani = animation.FuncAnimation(fig, update, frames=frames_list, interval=50, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'energy_collapse.mp4')

try:
    writer = animation.FFMpegWriter(fps=20, bitrate=2800)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=20), dpi=100)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

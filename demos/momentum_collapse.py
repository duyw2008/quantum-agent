#!/usr/bin/env python3
"""
自由粒子 — 动量测量坍缩 (增强版)

3 面板动画:
  左上: 位置空间 |ψ(x)|²  — 坍缩后剧烈展宽
  右上: 动量空间 |ψ̃(p)|² — 坍缩后收窄为单一频率
  下方: Δx 和 Δp 随时间演化 — 测量瞬间的跳跃
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm

print("=" * 60)
print("  Momentum Measurement Collapse (Enhanced)")
print("=" * 60)

# ============================================================
# 参数与演化
# ============================================================
grid = WaveGrid(x_min=-80, x_max=80, N=2048)
x, k = grid.x, grid.k

sigma_x, p0, t_measure = 0.5, 3.0, 3.0
psi0 = gaussian_wavepacket(grid, x0=-5.0, p0=p0, sigma=sigma_x)
psi0 /= np.sqrt(np.trapezoid(np.abs(psi0)**2, x))

result_pre = evolve_ssfm(psi0, grid, dt=0.003, t_max=t_measure, snapshots=80)
psi_before = result_pre['psi'][-1]

# 动量测量 + 坍缩
psi_k_before = np.fft.fft(psi_before)
prob_k_before = np.abs(psi_k_before)**2
np.random.seed(42)
measured_p = np.random.choice(k, p=prob_k_before / prob_k_before.sum())

sigma_p = 0.3
psi_k_after = np.exp(-(k - measured_p)**2 / (2 * sigma_p**2)).astype(complex)
psi_k_after /= np.sqrt(np.trapezoid(np.abs(psi_k_after)**2, k))
psi_after = np.fft.ifft(psi_k_after)
psi_after /= np.sqrt(np.trapezoid(np.abs(psi_after)**2, x))

t_after = 4.0
result_after = evolve_ssfm(psi_after, grid, dt=0.005, t_max=t_after, snapshots=100)

# ============================================================
# 统计
# ============================================================
dx_before = np.sqrt(np.trapezoid(x**2 * np.abs(psi_before)**2, x) -
                    np.trapezoid(x * np.abs(psi_before)**2, x)**2)
dp_before = np.sqrt(np.trapezoid(k**2 * prob_k_before, k) /
                    np.trapezoid(prob_k_before, k) -
                    (np.trapezoid(k * prob_k_before, k) /
                     np.trapezoid(prob_k_before, k))**2)

dx_after = np.sqrt(np.trapezoid(x**2 * np.abs(psi_after)**2, x) -
                   np.trapezoid(x * np.abs(psi_after)**2, x)**2)
dp_after = sigma_p / np.sqrt(2)  # 高斯 σ_p 在概率密度中的 Δp

print(f"\nBefore measurement: Δx={dx_before:.2f}, Δp={dp_before:.2f}")
print(f"Measured momentum:  p={measured_p:.2f}")
print(f"After collapse:     Δx={dx_after:.1f}, Δp={dp_after:.2f}")
print(f"Δp collapsed {dp_before/dp_after:.0f}×,  Δx expanded {dx_after/dx_before:.0f}×")

# ============================================================
# 3 面板动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 准备数据
def compute_momentum_probs(psi_list):
    return [np.abs(np.fft.fft(psi))**2 for psi in psi_list]

prob_k_pre = compute_momentum_probs(result_pre['psi'])
prob_k_after = compute_momentum_probs(result_after['psi'])

all_times_pre = np.array(result_pre['times'])
all_times_after = np.array(result_after['times']) + t_measure
all_times = np.concatenate([all_times_pre, all_times_after])
all_probs_x = np.concatenate([np.array(result_pre['prob']),
                               np.array(result_after['prob'])])
all_probs_k = np.concatenate([np.array(prob_k_pre), np.array(prob_k_after)])
measure_idx = len(all_times_pre) - 1

# 计算每帧的 Δx, Δp
dx_vals, dp_vals = [], []
for i in range(len(all_times)):
    px, pk = all_probs_x[i], all_probs_k[i]
    mx = np.trapezoid(x * px, x)
    dx_vals.append(np.sqrt(max(0, np.trapezoid(x**2 * px, x) - mx**2)))
    mk = np.trapezoid(k * pk, k) / np.trapezoid(pk, k)
    dp_vals.append(np.sqrt(max(0, np.trapezoid(k**2 * pk, k) /
                                np.trapezoid(pk, k) - mk**2)))

# === 图形布局 ===
fig = plt.figure(figsize=(16, 9), facecolor='#0d1117')
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.5],
                       hspace=0.35, wspace=0.25)

ax_x = fig.add_subplot(gs[0, 0], facecolor='#0d1117')   # 位置空间
ax_k = fig.add_subplot(gs[0, 1], facecolor='#0d1117')   # 动量空间
ax_t = fig.add_subplot(gs[1, :], facecolor='#0d1117')    # 时间演化

# --- 位置空间 ---
fill_x = ax_x.fill_between(x, 0, all_probs_x[0], alpha=0.5, color='#58a6ff')
line_x, = ax_x.plot(x, all_probs_x[0], color='#58a6ff', linewidth=1.2)
ax_x.set_xlabel('x (position)', color='#e6edf3', fontsize=11)
ax_x.set_ylabel(r'$|\psi(x)|^2$', color='#e6edf3', fontsize=11)
ax_x.set_title('Position Space', color='#e6edf3', fontsize=13, fontweight='bold')
ax_x.tick_params(colors='#e6edf3')
ax_x.grid(True, alpha=0.1, color='#30363d')
ax_x.set_xlim(x[0], x[-1])
ax_x.set_ylim(0, max(p.max() for p in all_probs_x) * 1.15)

# --- 动量空间 ---
# 聚焦在测量动量附近的区域
k_zoom = 8  # 显示范围 [measured_p - k_zoom, measured_p + k_zoom]
fill_k = ax_k.fill_between(k, 0, all_probs_k[0], alpha=0.5, color='#d2a8ff')
line_k, = ax_k.plot(k, all_probs_k[0], color='#d2a8ff', linewidth=1.2)
ax_k.set_xlabel('p (momentum / frequency)', color='#e6edf3', fontsize=11)
ax_k.set_ylabel(r'$|\tilde{\psi}(p)|^2$', color='#e6edf3', fontsize=11)
ax_k.set_title('Momentum Space', color='#e6edf3', fontsize=13, fontweight='bold')
ax_k.tick_params(colors='#e6edf3')
ax_k.grid(True, alpha=0.1, color='#30363d')
ax_k.set_xlim(measured_p - k_zoom, measured_p + k_zoom)
ax_k.set_ylim(0, max(p.max() for p in all_probs_k) * 1.15)

# 测量动量标记
meas_line = ax_k.axvline(measured_p, color='#f0883e', linewidth=2,
                          linestyle='--', alpha=0)
meas_text = ax_k.text(measured_p + 0.3, 0, '', color='#f0883e',
                       fontsize=10, fontweight='bold', va='bottom')

# 动量分布宽度标注
width_arrow_before = ax_k.annotate('', xy=(measured_p - dp_before/2, 0.02),
                                    xytext=(measured_p + dp_before/2, 0.02),
                                    arrowprops=dict(arrowstyle='<->', color='#d2a8ff',
                                                    lw=2, alpha=0))
width_label_before = ax_k.text(measured_p, 0.025, '', color='#d2a8ff',
                                fontsize=9, ha='center', alpha=0)

width_arrow_after = ax_k.annotate('', xy=(measured_p - dp_after/2, 0.06),
                                   xytext=(measured_p + dp_after/2, 0.06),
                                   arrowprops=dict(arrowstyle='<->', color='#ff7b72',
                                                   lw=3, alpha=0))
width_label_after = ax_k.text(measured_p, 0.065, '', color='#ff7b72',
                               fontsize=10, ha='center', fontweight='bold', alpha=0)

# --- 时间演化面板 ---
ax_t.plot(all_times_pre, dx_vals[:measure_idx+1], color='#58a6ff', linewidth=2, label=r'$\Delta x$')
ax_t.plot(all_times_pre, dp_vals[:measure_idx+1], color='#d2a8ff', linewidth=2, label=r'$\Delta p$',
          linestyle='--')
ax_t.plot(all_times_after, dx_vals[measure_idx+1:], color='#ff7b72', linewidth=2)
ax_t.plot(all_times_after, dp_vals[measure_idx+1:], color='#ff7b72', linewidth=2, linestyle='--')

# 测量时刻垂直线
ax_t.axvline(t_measure, color='#f0883e', linewidth=2.5, linestyle='--', alpha=0.8)
ax_t.text(t_measure + 0.1, max(max(dx_vals), max(dp_vals)) * 0.85,
          'MEASUREMENT', color='#f0883e', fontsize=11, fontweight='bold')

# 标注 Δp 跳跃
mid_y = (dp_vals[measure_idx] + dp_vals[measure_idx+1]) / 2
ax_t.annotate(r'$\Delta p$ ↓' + f'{dp_before/dp_after:.0f}×',
              xy=(t_measure, dp_vals[measure_idx+1]),
              xytext=(t_measure + 0.8, dp_vals[measure_idx+1] * 1.5),
              arrowprops=dict(arrowstyle='->', color='#f0883e', lw=1.5),
              color='#f0883e', fontsize=10, fontweight='bold')
ax_t.annotate(r'$\Delta x$ ↑' + f'{dx_after/dx_before:.0f}×',
              xy=(t_measure, dx_vals[measure_idx+1]),
              xytext=(t_measure + 0.8, dx_vals[measure_idx+1] * 0.6),
              arrowprops=dict(arrowstyle='->', color='#f0883e', lw=1.5),
              color='#f0883e', fontsize=10, fontweight='bold')

t_cursor = ax_t.axvline(0, color='#e6edf3', linewidth=1, alpha=0.3)
ax_t.set_xlabel('Time', color='#e6edf3', fontsize=11)
ax_t.set_ylabel('Uncertainty', color='#e6edf3', fontsize=11)
ax_t.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3',
            loc='upper left')
ax_t.tick_params(colors='#e6edf3')
ax_t.grid(True, alpha=0.1, color='#30363d')

# 全局标题
fig_title = fig.suptitle('', color='#e6edf3', fontsize=15, fontweight='bold', y=0.98)

for ax in [ax_x, ax_k, ax_t]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

# --- 更新函数 ---
fill_x_ref, fill_k_ref = [fill_x], [fill_k]

def update(i):
    fill_x_ref[0].remove()
    fill_k_ref[0].remove()

    is_after = i > measure_idx
    cx = '#ff7b72' if is_after else '#58a6ff'
    ck = '#ff7b72' if is_after else '#d2a8ff'

    fill_x_ref[0] = ax_x.fill_between(x, 0, all_probs_x[i], alpha=0.5, color=cx)
    line_x.set_ydata(all_probs_x[i]); line_x.set_color(cx)

    fill_k_ref[0] = ax_k.fill_between(k, 0, all_probs_k[i], alpha=0.5, color=ck)
    line_k.set_ydata(all_probs_k[i]); line_k.set_color(ck)

    t_cursor.set_xdata([all_times[i], all_times[i]])

    if is_after:
        meas_line.set_alpha(0.7)
        meas_text.set_text(f'p₀={measured_p:.2f}')
        meas_text.set_y(all_probs_k[i].max() * 0.85)
        width_arrow_before.set_alpha(0.3)
        width_label_before.set_alpha(0.3)
        width_label_before.set_text(f'Δp={dp_before:.2f} (before)')
        width_arrow_after.set_alpha(1.0)
        width_label_after.set_alpha(1.0)
        width_label_after.set_text(f'Δp={dp_after:.2f} (now)')
    else:
        meas_line.set_alpha(0)
        meas_text.set_text('')
        width_arrow_before.set_alpha(0.7)
        width_label_before.set_alpha(0.7)
        width_label_before.set_text(f'Δp={dp_vals[i]:.2f}')
        width_arrow_after.set_alpha(0)
        width_label_after.set_alpha(0)

    t = all_times[i]
    phase = '◉ MOMENTUM COLLAPSED' if is_after else '○ FREE EVOLUTION'
    fig_title.set_text(
        f't = {t:.2f}    |    Δx = {dx_vals[i]:.2f}    |    '
        f'Δp = {dp_vals[i]:.3f}    |    {phase}')

    return [line_x, line_k, t_cursor]

total_frames = len(all_times)
extended = []
for i in range(total_frames):
    extended.append(i)
    if i == measure_idx:
        for _ in range(20):
            extended.append(i)

ani = animation.FuncAnimation(fig, update, frames=extended,
                               interval=40, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'momentum_collapse.mp4')

try:
    writer = animation.FFMpegWriter(fps=25, bitrate=2500)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=25), dpi=150)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

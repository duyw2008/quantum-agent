#!/usr/bin/env python3
"""
海森堡不确定性原理 — 动画演示

Δx·Δp ≥ ℏ/2 的连续时间演化：
  自由高斯波包: Δx(t) 增长, Δp 恒定 → Δx·Δp 持续增大
  最小不确定态仅在 t=0 达到

4 面板布局:
  左上: 位置空间 |ψ(x)|²  +  Δx 标注
  右上: 动量空间 |ψ̃(p)|²  +  Δp 标注
  左下: Δx(t), Δp(t) 时间演化
  右下: Δx·Δp 时间演化, 标注 ℏ/2 下界
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qm.wave import WaveGrid, gaussian_wavepacket, evolve_ssfm

print("=" * 60)
print("  Heisenberg Uncertainty Principle — Animation")
print("=" * 60)

# ============================================================
# 参数: 最小不确定高斯波包
# ============================================================
sigma = 1.0        # 初始 Δx
hbar = 1.0
mass = 1.0

# 特征弥散时间
tau = 2 * mass * sigma**2 / hbar
print(f"Spreading timescale τ = {tau:.1f}")

# 网格 — 需要足够大容纳弥散后的波包
grid = WaveGrid(x_min=-40, x_max=40, N=1024)
x, k = grid.x, grid.k

psi0 = gaussian_wavepacket(grid, x0=0.0, p0=0.0, sigma=sigma)
psi0 /= np.sqrt(np.trapezoid(np.abs(psi0)**2, x))

# 演化 5 个 τ
t_max = 5 * tau
result = evolve_ssfm(psi0, grid, dt=0.01, t_max=t_max, snapshots=200)

# ============================================================
# 计算每帧的不确定度
# ============================================================
times = np.array(result['times'])
n_frames = len(times)

dx_vals = np.zeros(n_frames)
dp_vals = np.zeros(n_frames)

for i in range(n_frames):
    psi = result['psi'][i]
    prob_x = np.abs(psi)**2
    mx = np.trapezoid(x * prob_x, x)
    dx_vals[i] = np.sqrt(np.trapezoid(x**2 * prob_x, x) - mx**2)

    psi_k = np.fft.fft(psi)
    prob_k = np.abs(psi_k)**2
    mk = np.trapezoid(k * prob_k, k) / np.trapezoid(prob_k, k)
    dp_vals[i] = np.sqrt(np.trapezoid(k**2 * prob_k, k) /
                          np.trapezoid(prob_k, k) - mk**2)

product = dx_vals * dp_vals

# 理论值
dx_theory = sigma * np.sqrt(1 + (times / tau)**2)
dp_theory = np.full_like(times, hbar / (2 * sigma))
product_theory = dx_theory * dp_theory

print(f"Initial:  Δx={dx_vals[0]:.3f}, Δp={dp_vals[0]:.3f}, Δx·Δp={product[0]:.4f}")
print(f"Final:    Δx={dx_vals[-1]:.3f}, Δp={dp_vals[-1]:.3f}, Δx·Δp={product[-1]:.4f}")
print(f"Minimum:  ℏ/2 = {hbar/2:.4f}")

# ============================================================
# 4 面板动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(16, 10), facecolor='#0d1117')
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3,
                       height_ratios=[2.5, 1.5])

ax_x = fig.add_subplot(gs[0, 0], facecolor='#0d1117')    # 位置空间
ax_k = fig.add_subplot(gs[0, 1], facecolor='#0d1117')    # 动量空间
ax_d = fig.add_subplot(gs[1, 0], facecolor='#0d1117')    # Δx, Δp vs t
ax_h = fig.add_subplot(gs[1, 1], facecolor='#0d1117')    # Δx·Δp vs t

# --- 位置空间 ---
prob_x = np.abs(result['psi'][0])**2
fill_x = ax_x.fill_between(x, 0, prob_x, alpha=0.5, color='#58a6ff')
line_x, = ax_x.plot(x, prob_x, color='#58a6ff', linewidth=1.5)
# Δx 标注线
dx_left = ax_x.axvline(-dx_vals[0], color='#f0883e', linewidth=1.5,
                         linestyle='--', alpha=0.7)
dx_right = ax_x.axvline(dx_vals[0], color='#f0883e', linewidth=1.5,
                          linestyle='--', alpha=0.7)
dx_label = ax_x.text(0.05, 0.92, '', transform=ax_x.transAxes,
                      color='#f0883e', fontsize=11, fontweight='bold')
ax_x.set_xlabel('x', color='#e6edf3', fontsize=12)
ax_x.set_ylabel(r'$|\psi(x)|^2$', color='#e6edf3', fontsize=12)
ax_x.set_title('Position Space', color='#e6edf3', fontsize=14, fontweight='bold')
ax_x.tick_params(colors='#e6edf3')
ax_x.grid(True, alpha=0.1, color='#30363d')
ax_x.set_xlim(x[0], x[-1])

# --- 动量空间 ---
k_zoom = 6
prob_k0 = np.abs(np.fft.fft(result['psi'][0]))**2
fill_k = ax_k.fill_between(k, 0, prob_k0, alpha=0.5, color='#d2a8ff')
line_k, = ax_k.plot(k, prob_k0, color='#d2a8ff', linewidth=1.5)
dp_left = ax_k.axvline(-dp_vals[0], color='#f0883e', linewidth=1.5,
                         linestyle='--', alpha=0.7)
dp_right = ax_k.axvline(dp_vals[0], color='#f0883e', linewidth=1.5,
                          linestyle='--', alpha=0.7)
dp_label = ax_k.text(0.05, 0.92, '', transform=ax_k.transAxes,
                      color='#f0883e', fontsize=11, fontweight='bold')
ax_k.set_xlabel('p', color='#e6edf3', fontsize=12)
ax_k.set_ylabel(r'$|\tilde{\psi}(p)|^2$', color='#e6edf3', fontsize=12)
ax_k.set_title('Momentum Space', color='#e6edf3', fontsize=14, fontweight='bold')
ax_k.tick_params(colors='#e6edf3')
ax_k.grid(True, alpha=0.1, color='#30363d')
ax_k.set_xlim(-k_zoom, k_zoom)

# --- 不确定度时间演化 ---
ax_d.plot(times, dx_vals, color='#58a6ff', linewidth=2.5, label=r'$\Delta x$')
ax_d.plot(times, dp_vals, color='#d2a8ff', linewidth=2.5, label=r'$\Delta p$')
# 理论线
ax_d.plot(times, dx_theory, color='#58a6ff', linewidth=1, linestyle=':', alpha=0.5)
ax_d.plot(times, dp_theory, color='#d2a8ff', linewidth=1, linestyle=':', alpha=0.5)
t_cur = ax_d.axvline(0, color='#e6edf3', linewidth=1, alpha=0.3)
ax_d.set_xlabel('Time', color='#e6edf3', fontsize=12)
ax_d.set_ylabel('Uncertainty', color='#e6edf3', fontsize=12)
ax_d.set_title(r'$\Delta x(t)$ and $\Delta p(t)$', color='#e6edf3', fontsize=13)
ax_d.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3')
ax_d.tick_params(colors='#e6edf3')
ax_d.grid(True, alpha=0.1, color='#30363d')

# --- Δx·Δp 乘积 ---
ax_h.axhline(hbar/2, color='#f0883e', linewidth=2, linestyle='-', alpha=0.8)
ax_h.text(times[-1] * 0.6, hbar/2 * 1.05, r'$\hbar/2$ (Heisenberg limit)',
          color='#f0883e', fontsize=11, ha='center')
ax_h.plot(times, product, color='#79c0ff', linewidth=2.5)
ax_h.plot(times, product_theory, color='#79c0ff', linewidth=1, linestyle=':', alpha=0.5)
ax_h.fill_between(times, hbar/2, product, alpha=0.15, color='#f0883e')
ax_h.fill_between(times, 0, hbar/2, alpha=0.3, color='#ff7b72')
ax_h.text(times[len(times)//2], hbar/4, 'FORBIDDEN', color='#ff7b72',
          fontsize=10, ha='center', alpha=0.5)
t_cur2 = ax_h.axvline(0, color='#e6edf3', linewidth=1, alpha=0.3)
ax_h.set_xlabel('Time', color='#e6edf3', fontsize=12)
ax_h.set_ylabel(r'$\Delta x \cdot \Delta p$', color='#e6edf3', fontsize=12)
ax_h.set_title(r'Uncertainty Product', color='#e6edf3', fontsize=13)
ax_h.tick_params(colors='#e6edf3')
ax_h.grid(True, alpha=0.1, color='#30363d')

# 全局标题
fig_title = fig.suptitle('', color='#e6edf3', fontsize=16, fontweight='bold', y=0.99)

# 样式
for ax in [ax_x, ax_k, ax_d, ax_h]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

# y 轴范围
all_x_max = max(np.max(np.abs(result['psi'][i]))**2 for i in range(0, n_frames, 5))
ax_x.set_ylim(0, all_x_max * 1.1)
all_k_max = max(np.max(np.abs(np.fft.fft(result['psi'][i])))**2 for i in range(0, n_frames, 5))
ax_k.set_ylim(0, all_k_max * 1.1)

fill_x_ref = [fill_x]; fill_k_ref = [fill_k]

def update(i):
    fill_x_ref[0].remove(); fill_k_ref[0].remove()

    psi = result['psi'][i]
    prob_x = np.abs(psi)**2
    prob_k = np.abs(np.fft.fft(psi))**2

    fill_x_ref[0] = ax_x.fill_between(x, 0, prob_x, alpha=0.5, color='#58a6ff')
    line_x.set_ydata(prob_x)
    dx_left.set_xdata([-dx_vals[i], -dx_vals[i]])
    dx_right.set_xdata([dx_vals[i], dx_vals[i]])
    dx_label.set_text(f'Δx = {dx_vals[i]:.3f}')

    fill_k_ref[0] = ax_k.fill_between(k, 0, prob_k, alpha=0.5, color='#d2a8ff')
    line_k.set_ydata(prob_k)
    dp_left.set_xdata([-dp_vals[i], -dp_vals[i]])
    dp_right.set_xdata([dp_vals[i], dp_vals[i]])
    dp_label.set_text(f'Δp = {dp_vals[i]:.3f}')

    t_cur.set_xdata([times[i], times[i]])
    t_cur2.set_xdata([times[i], times[i]])

    t = times[i]
    ratio = product[i] / (hbar/2)
    fig_title.set_text(
        f't = {t:.2f}    |    Δx = {dx_vals[i]:.3f}    |    Δp = {dp_vals[i]:.3f}    |    '
        f'Δx·Δp = {product[i]:.4f} = {ratio:.2f} × ℏ/2')

    return [line_x, line_k, t_cur, t_cur2]

ani = animation.FuncAnimation(fig, update, frames=n_frames,
                               interval=50, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'heisenberg_uncertainty.mp4')

try:
    writer = animation.FFMpegWriter(fps=20, bitrate=2500)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=20), dpi=150)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

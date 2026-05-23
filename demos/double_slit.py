#!/usr/bin/env python3
"""
双缝干涉实验 — 2D 含时薛定谔方程模拟

高斯波包从左侧入射 → 穿过双缝 → 右侧形成干涉条纹。

方法: 2D Split-Step Fourier Method
网格: 256×128, 势函数 = 遮挡屏 (x=0 处, y 有两缝开口)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np

print("=" * 60)
print("  Double-Slit Interference — 2D TDSE Simulation")
print("=" * 60)

# ============================================================
# 2D 网格与势函数
# ============================================================
Nx, Ny = 256, 128
x = np.linspace(-20, 30, Nx)
y = np.linspace(-8, 8, Ny)
dx, dy = x[1] - x[0], y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing='ij')

# 动量空间
kx = 2 * np.pi * np.fft.fftfreq(Nx, dx)
ky = 2 * np.pi * np.fft.fftfreq(Ny, dy)
KX, KY = np.meshgrid(kx, ky, indexing='ij')

# 势函数: 屏障在 x=0, 两条缝在 y = ±d
barrier_x = 0.0
slit_sep = 4.0
slit_width = 0.6
barrier_height = 50.0

V = np.zeros((Nx, Ny))
barrier_mask = np.abs(x - barrier_x) < 0.3
V[barrier_mask, :] = barrier_height
# 开缝
for slit_center in [-slit_sep/2, slit_sep/2]:
    slit_mask = np.abs(y - slit_center) < slit_width / 2
    V[np.ix_(barrier_mask, slit_mask)] = 0.0

print(f"Grid: {Nx}×{Ny}, dx={dx:.3f}, dy={dy:.3f}")
print(f"Slits: y = ±{slit_sep/2:.1f}, width = {slit_width}")
print(f"Barrier: x = {barrier_x}, height = {barrier_height}")

# ============================================================
# 初始波包
# ============================================================
x0, p0 = -8.0, 4.0
sigma_x, sigma_y = 1.5, 2.5
hbar, mass = 1.0, 1.0

psi = np.exp(-(X - x0)**2 / (2 * sigma_x**2) -
             Y**2 / (2 * sigma_y**2) +
             1j * p0 * (X - x0) / hbar)
psi /= np.sqrt(np.trapezoid(np.trapezoid(np.abs(psi)**2, y), x))

E0 = p0**2 / (2 * mass)
print(f"Wavepacket: x₀={x0}, p₀={p0}, E={E0:.1f}")
print(f"σ_x={sigma_x}, σ_y={sigma_y}")

# ============================================================
# 2D SSFM 演化
# ============================================================
dt = 0.003
t_max = 12.0
n_steps = int(t_max / dt)
snapshots = 200
interval = max(1, n_steps // snapshots)

# 预计算相位因子
pe_half = np.exp(-0.5j * V * dt / hbar)
ke_full = np.exp(-1.0j * hbar * (KX**2 + KY**2) * dt / (2 * mass))

times, probs, slices = [], [], []  # slices = |ψ|² at y=0 (center)

def save_state(t):
    prob = np.abs(psi)**2
    times.append(t)
    probs.append(prob)
    slices.append(prob[:, Ny//2])  # 中心切片

save_state(0)
print(f"Evolving {n_steps} steps...")

for step in range(1, n_steps + 1):
    psi = pe_half * psi
    psi_k = np.fft.fft2(psi)
    psi_k *= ke_full
    psi = np.fft.ifft2(psi_k)
    psi = pe_half * psi

    if step % interval == 0:
        save_state(step * dt)
    if step % (n_steps // 10) == 0:
        print(f"  {step/n_steps*100:.0f}%", end=' ', flush=True)

print(f"\nDone: {len(times)} snapshots")

# ============================================================
# 提取干涉图案 — 右侧屏幕上的概率密度
# ============================================================
screen_x = 8.0  # 屏幕位置
screen_idx = np.argmin(np.abs(x - screen_x))
screen_pattern = probs[-1][screen_idx, :]
print(f"Screen at x = {x[screen_idx]:.1f}")

# ============================================================
# 3 面板动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(16, 9), facecolor='#0d1117')
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.5],
                       hspace=0.4, wspace=0.3)

ax_2d = fig.add_subplot(gs[0, :], facecolor='#0d1117')      # 2D 概率密度
ax_scr = fig.add_subplot(gs[1, 0], facecolor='#0d1117')     # 屏幕图案
ax_slc = fig.add_subplot(gs[1, 1], facecolor='#0d1117')     # 中心切片

# --- 2D 热图 ---
vmax = max(p.max() for p in probs)
im = ax_2d.imshow(probs[0].T, extent=[x[0], x[-1], y[0], y[-1]],
                   origin='lower', cmap='plasma', aspect='auto',
                   vmin=0, vmax=vmax * 0.8)
# 画屏障位置
ax_2d.axvline(barrier_x, color='white', linewidth=1, alpha=0.3, linestyle='--')
# 画屏幕位置
ax_2d.axvline(screen_x, color='#f0883e', linewidth=1.5, alpha=0.6, linestyle='--')
ax_2d.text(screen_x + 0.2, y[-1] * 0.9, 'SCREEN', color='#f0883e',
           fontsize=9, fontweight='bold')
ax_2d.set_xlabel('x', color='#e6edf3', fontsize=12)
ax_2d.set_ylabel('y', color='#e6edf3', fontsize=12)
ax_2d.set_title('Double-Slit Interference', color='#e6edf3',
                fontsize=14, fontweight='bold')
ax_2d.tick_params(colors='#e6edf3')
cbar = plt.colorbar(im, ax=ax_2d, fraction=0.015)
cbar.set_label(r'$|\psi|^2$', color='#e6edf3')
cbar.ax.tick_params(colors='#e6edf3')

# --- 屏幕图案 ---
ax_scr.set_facecolor('#0d1117')
line_scr, = ax_scr.plot(y, screen_pattern, color='#79c0ff', linewidth=2)
ax_scr.fill_between(y, 0, screen_pattern, alpha=0.3, color='#79c0ff')
# 理论干涉极大位置 (双缝干涉: d sinθ = nλ)
lam = 2 * np.pi * hbar / p0  # de Broglie 波长
L = screen_x - barrier_x      # 缝到屏幕距离
for n in range(-3, 4):
    if n == 0:
        continue
    y_max = n * lam * L / slit_sep
    if abs(y_max) < y[-1]:
        ax_scr.axvline(y_max, color='#f0883e', linewidth=0.8,
                       linestyle=':', alpha=0.4)
ax_scr.set_xlabel('y (screen position)', color='#e6edf3', fontsize=11)
ax_scr.set_ylabel(r'$|\psi|^2$ on screen', color='#e6edf3', fontsize=11)
ax_scr.set_title(f'Interference Pattern at x={screen_x}', color='#e6edf3', fontsize=12)
ax_scr.tick_params(colors='#e6edf3')
ax_scr.grid(True, alpha=0.1, color='#30363d')

# --- 中心切片 ---
ax_slc.set_facecolor('#0d1117')
line_slc, = ax_slc.plot(x, slices[0], color='#58a6ff', linewidth=1.5)
ax_slc.fill_between(x, 0, slices[0], alpha=0.3, color='#58a6ff')
ax_slc.set_xlabel('x', color='#e6edf3', fontsize=11)
ax_slc.set_ylabel(r'$|\psi(x,0)|^2$', color='#e6edf3', fontsize=11)
ax_slc.set_title('Center Slice (y=0)', color='#e6edf3', fontsize=12)
ax_slc.tick_params(colors='#e6edf3')
ax_slc.grid(True, alpha=0.1, color='#30363d')
ax_slc.set_ylim(0, max(s.max() for s in slices) * 1.1)
t_cursor = ax_2d.axvline(x[0], color='white', linewidth=2, alpha=0.5)

# 全局标题
fig_title = fig.suptitle('', color='#e6edf3', fontsize=15, fontweight='bold', y=0.99)

for ax in [ax_2d, ax_scr, ax_slc]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

def update(i):
    im.set_data(probs[i].T)
    line_slc.set_ydata(slices[i])
    # 仅最后一帧更新屏幕图案
    if i == len(times) - 1:
        line_scr.set_ydata(probs[i][screen_idx, :])
    t_cursor.set_xdata([x[0] + p0 * times[i] / mass,
                         x[0] + p0 * times[i] / mass])
    t = times[i]
    fig_title.set_text(f't = {t:.2f}    |    de Broglie λ = {lam:.3f}')
    return [im, line_slc, t_cursor]

ani = animation.FuncAnimation(fig, update, frames=len(times),
                               interval=50, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'double_slit.mp4')

try:
    writer = animation.FFMpegWriter(fps=20, bitrate=3000)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=20), dpi=120)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

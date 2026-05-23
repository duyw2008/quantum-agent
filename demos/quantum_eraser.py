#!/usr/bin/env python3
"""
延迟选择量子擦除实验 — Wheeler's Delayed Choice

双缝后放置 BBO 晶体 → 产生纠缠光子对 → 分别检测:
  "which-path" (非相干): 可区分路径 → |ψ₁|² + |ψ₂|² → 无干涉
  "eraser"     (相干):   路径信息擦除 → |ψ₁ + ψ₂|² → 干涉条纹

分屏对比: 上半 = eraser (干涉), 下半 = which-path (无干涉)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np

print("=" * 60)
print("  Delayed-Choice Quantum Eraser")
print("=" * 60)

# ============================================================
# 网格与势函数 (同双缝实验)
# ============================================================
Nx, Ny = 256, 128
x = np.linspace(-20, 30, Nx)
y = np.linspace(-8, 8, Ny)
dx, dy = x[1] - x[0], y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing='ij')

kx = 2 * np.pi * np.fft.fftfreq(Nx, dx)
ky = 2 * np.pi * np.fft.fftfreq(Ny, dy)
KX, KY = np.meshgrid(kx, ky, indexing='ij')

hbar, mass = 1.0, 1.0

# 屏障
barrier_x = 0.0; slit_sep = 4.0; slit_width = 0.6
barrier_height = 50.0

def make_potential(open_upper=True, open_lower=True):
    V = np.zeros((Nx, Ny))
    barrier_mask = np.abs(x - barrier_x) < 0.3
    V[barrier_mask, :] = barrier_height
    if open_upper:
        V[np.ix_(barrier_mask, np.abs(y - slit_sep/2) < slit_width/2)] = 0.0
    if open_lower:
        V[np.ix_(barrier_mask, np.abs(y + slit_sep/2) < slit_width/2)] = 0.0
    return V

V_both  = make_potential(True, True)    # 双缝 (相干)
V_upper = make_potential(True, False)   # 仅上缝
V_lower = make_potential(False, True)   # 仅下缝

# 初始波包
x0, p0 = -10.0, 6.0
sigma_x = sigma_y = 2.0

def make_psi():
    psi = np.exp(-(X - x0)**2 / (2 * sigma_x**2) -
                 Y**2 / (2 * sigma_y**2) +
                 1j * p0 * (X - x0) / hbar)
    return psi / np.sqrt(np.trapezoid(np.trapezoid(np.abs(psi)**2, y), x))

psi0 = make_psi()

# SSFM 参数
dt = 0.003; t_max = 12.0
n_steps = int(t_max / dt)
snapshots = 150
interval = max(1, n_steps // snapshots)

def evolve(psi0, V):
    psi = psi0.copy()
    pe_half = np.exp(-0.5j * V * dt / hbar)
    ke_full = np.exp(-1.0j * hbar * (KX**2 + KY**2) * dt / (2 * mass))
    times, probs = [], []
    def save(t):
        times.append(t); probs.append(np.abs(psi)**2)
    save(0)
    for step in range(1, n_steps + 1):
        psi = pe_half * psi
        psi_k = np.fft.fft2(psi); psi_k *= ke_full
        psi = np.fft.ifft2(psi_k); psi = pe_half * psi
        if step % interval == 0:
            save(step * dt)
    return np.array(times), np.array(probs)

print(f"Running 3 simulations ({n_steps} steps each)...")

# 三个模拟: 相干(双缝) + 非相干分量(单缝×2)
print("  [1/3] Coherent (both slits)...")
t_coh, p_coh = evolve(psi0.copy(), V_both)
print("  [2/3] Upper slit only...")
t_upp, p_upp = evolve(psi0.copy(), V_upper)
print("  [3/3] Lower slit only...")
t_low, p_low = evolve(psi0.copy(), V_lower)

# 非相干 = |ψ_upper|² + |ψ_lower|² (概率相加, 非振幅)
p_incoh = p_upp + p_low

print(f"Done: {len(t_coh)} snapshots each")

# ============================================================
# 屏幕图案
# ============================================================
screen_x = 10.0
screen_idx = np.argmin(np.abs(x - screen_x))
coh_pattern = p_coh[-1, screen_idx, :]
incoh_pattern = p_incoh[-1, screen_idx, :]
lam = 2 * np.pi * hbar / p0
L = screen_x - barrier_x

print(f"Screen at x={x[screen_idx]:.1f}, λ={lam:.3f}, L={L:.1f}")
print(f"Fringe spacing: {lam*L/slit_sep:.2f}")

# ============================================================
# 动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(18, 10), facecolor='#0d1117')
gs = fig.add_gridspec(3, 2, height_ratios=[2.5, 2.5, 1.2],
                       hspace=0.45, wspace=0.3)

# 上: 相干 (eraser) — 有干涉
ax_coh = fig.add_subplot(gs[0, :], facecolor='#0d1117')
# 中: 非相干 (which-path) — 无干涉
ax_inc = fig.add_subplot(gs[1, :], facecolor='#0d1117')
# 下左: 屏幕图案
ax_scr = fig.add_subplot(gs[2, 0], facecolor='#0d1117')
# 下右: 对比
ax_cmp = fig.add_subplot(gs[2, 1], facecolor='#0d1117')

gamma = 0.45
vmax = max(p_coh.max(), p_incoh.max())

def to_display(prob):
    return (np.clip(prob / (vmax * 0.12), 0, 1)) ** gamma

# --- 相干面板 ---
im_coh = ax_coh.imshow(to_display(p_coh[0]).T,
                        extent=[x[0], x[-1], y[0], y[-1]],
                        origin='lower', cmap='inferno', aspect='equal')
ax_coh.axvline(barrier_x, color='white', linewidth=1, alpha=0.3, linestyle='--')
ax_coh.axvline(screen_x, color='#f0883e', linewidth=1.5, alpha=0.6, linestyle='--')
ax_coh.text(screen_x + 0.2, y[-1]*0.85, 'DETECTOR', color='#f0883e', fontsize=9, fontweight='bold')
# BBO 标记
ax_coh.fill_between([barrier_x+0.5, barrier_x+1.5], y[0], y[-1],
                     color='#d2a8ff', alpha=0.08)
ax_coh.text(barrier_x+0.7, y[0]+0.5, 'BBO', color='#d2a8ff', fontsize=8, alpha=0.6)
ax_coh.set_ylabel('y', color='#e6edf3', fontsize=11)
ax_coh.set_title('QUANTUM ERASER  (coherent: |ψ₁+ψ₂|²)  →  INTERFERENCE',
                 color='#58a6ff', fontsize=13, fontweight='bold')
ax_coh.tick_params(colors='#e6edf3')
cbar1 = plt.colorbar(im_coh, ax=ax_coh, fraction=0.012)
cbar1.ax.tick_params(colors='#e6edf3')

# --- 非相干面板 ---
im_inc = ax_inc.imshow(to_display(p_incoh[0]).T,
                        extent=[x[0], x[-1], y[0], y[-1]],
                        origin='lower', cmap='inferno', aspect='equal')
ax_inc.axvline(barrier_x, color='white', linewidth=1, alpha=0.3, linestyle='--')
ax_inc.axvline(screen_x, color='#f0883e', linewidth=1.5, alpha=0.6, linestyle='--')
ax_inc.text(screen_x + 0.2, y[-1]*0.85, 'DETECTOR', color='#f0883e', fontsize=9, fontweight='bold')
ax_inc.fill_between([barrier_x+0.5, barrier_x+1.5], y[0], y[-1],
                     color='#d2a8ff', alpha=0.08)
ax_inc.text(barrier_x+0.7, y[0]+0.5, 'BBO', color='#d2a8ff', fontsize=8, alpha=0.6)
ax_inc.set_xlabel('x', color='#e6edf3', fontsize=11)
ax_inc.set_ylabel('y', color='#e6edf3', fontsize=11)
ax_inc.set_title('WHICH-PATH  (incoherent: |ψ₁|²+|ψ₂|²)  →  NO INTERFERENCE',
                 color='#ff7b72', fontsize=13, fontweight='bold')
ax_inc.tick_params(colors='#e6edf3')
cbar2 = plt.colorbar(im_inc, ax=ax_inc, fraction=0.012)
cbar2.ax.tick_params(colors='#e6edf3')

# --- 屏幕图案 ---
ax_scr.set_facecolor('#0d1117')
line_coh, = ax_scr.plot(y, coh_pattern, color='#58a6ff', linewidth=2, label='Eraser')
line_inc, = ax_scr.plot(y, incoh_pattern, color='#ff7b72', linewidth=2, label='Which-path')
ax_scr.fill_between(y, 0, coh_pattern, alpha=0.2, color='#58a6ff')
ax_scr.set_xlabel('y (screen)', color='#e6edf3', fontsize=10)
ax_scr.set_ylabel(r'$|\psi|^2$', color='#e6edf3', fontsize=10)
ax_scr.set_title(f'Screen Pattern at x={screen_x}', color='#e6edf3', fontsize=11)
ax_scr.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
ax_scr.tick_params(colors='#e6edf3')
ax_scr.grid(True, alpha=0.1, color='#30363d')
# 理论条纹位置
for n in range(-4, 5):
    ym = n * lam * L / slit_sep
    if abs(ym) < y[-1]:
        ax_scr.axvline(ym, color='#f0883e', linewidth=0.6, linestyle=':', alpha=0.25)

# --- 下方右: 可见度 vs 时间 ---
ax_cmp.set_facecolor('#0d1117')
# 计算每帧的可见度
visibility = []
for i in range(len(t_coh)):
    pat_coh = p_coh[i, screen_idx, :]
    pat_inc = p_incoh[i, screen_idx, :]
    vis = max(0, (pat_coh.max() - pat_coh.min()) / (pat_coh.max() + pat_coh.min() + 1e-15))
    visibility.append(vis)

ax_cmp.plot(t_coh, visibility, color='#58a6ff', linewidth=2)
ax_cmp.axhline(0, color='#30363d', linewidth=0.5)
ax_cmp.set_xlabel('Time', color='#e6edf3', fontsize=10)
ax_cmp.set_ylabel('Fringe Visibility', color='#e6edf3', fontsize=10)
ax_cmp.set_title('Interference Visibility vs Time', color='#e6edf3', fontsize=11)
ax_cmp.tick_params(colors='#e6edf3')
ax_cmp.grid(True, alpha=0.1, color='#30363d')
t_cur = ax_cmp.axvline(0, color='#e6edf3', linewidth=1, alpha=0.3)

# 全局标题
fig_title = fig.suptitle('', color='#e6edf3', fontsize=16, fontweight='bold', y=0.99)

for ax in [ax_coh, ax_inc, ax_scr, ax_cmp]:
    for spine in ax.spines.values():
        spine.set_color('#30363d')

def update(i):
    im_coh.set_data(to_display(p_coh[i]).T)
    im_inc.set_data(to_display(p_incoh[i]).T)

    # 更新屏幕图案
    line_coh.set_ydata(p_coh[i, screen_idx, :])
    line_inc.set_ydata(p_incoh[i, screen_idx, :])

    t_cur.set_xdata([t_coh[i], t_coh[i]])
    t = t_coh[i]
    vis = visibility[i]
    fig_title.set_text(
        f't = {t:.2f}    |    Fringe Visibility = {vis:.3f}    |    '
        f'λ = {lam:.3f}    |    Δy = {lam*L/slit_sep:.2f}')

    return [im_coh, im_inc, line_coh, line_inc, t_cur]

ani = animation.FuncAnimation(fig, update, frames=len(t_coh),
                               interval=60, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'quantum_eraser.mp4')

try:
    writer = animation.FFMpegWriter(fps=16, bitrate=3500)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=16), dpi=100)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

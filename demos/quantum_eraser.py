#!/usr/bin/env python3
"""
量子擦除实验 — 增强版 (直观对比)

改进:
  1. 单一场景, 中途切换: t < t_erase → which-path (无条纹)
                         t > t_erase → eraser (条纹浮现!)
  2. 右下角: 干涉项 2Re(ψ₁*ψ₂) — 直观展示条纹来源
  3. 屏幕图案实时累积, 模拟探测器逐帧记录
  4. 切换时刻视觉提示
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np

print("=" * 60)
print("  Quantum Eraser — Enhanced")
print("=" * 60)

# ============================================================
# 网格与模拟 (同前)
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

barrier_x = 0.0; slit_sep = 4.0; slit_width = 0.6
barrier_height = 50.0

def make_potential(open_upper, open_lower):
    V = np.zeros((Nx, Ny))
    bm = np.abs(x - barrier_x) < 0.3
    V[bm, :] = barrier_height
    if open_upper: V[np.ix_(bm, np.abs(y - slit_sep/2) < slit_width/2)] = 0.0
    if open_lower: V[np.ix_(bm, np.abs(y + slit_sep/2) < slit_width/2)] = 0.0
    return V

V_both  = make_potential(True, True)
V_upper = make_potential(True, False)
V_lower = make_potential(False, True)

x0, p0 = -10.0, 6.0; sigma_x = sigma_y = 2.0
def make_psi():
    psi = np.exp(-(X-x0)**2/(2*sigma_x**2) - Y**2/(2*sigma_y**2) + 1j*p0*(X-x0)/hbar)
    return psi / np.sqrt(np.trapezoid(np.trapezoid(np.abs(psi)**2, y), x))

psi0 = make_psi()
dt, t_max = 0.003, 14.0; n_steps = int(t_max/dt)
snapshots = 200; interval = max(1, n_steps // snapshots)

def evolve(psi0, V):
    psi = psi0.copy()
    pe_half = np.exp(-0.5j*V*dt/hbar)
    ke_full = np.exp(-1.0j*hbar*(KX**2+KY**2)*dt/(2*mass))
    times, probs = [], []
    def save(t): times.append(t); probs.append(np.abs(psi)**2)
    save(0)
    for step in range(1, n_steps+1):
        psi = pe_half*psi; psi_k = np.fft.fft2(psi); psi_k *= ke_full
        psi = np.fft.ifft2(psi_k); psi = pe_half*psi
        if step % interval == 0: save(step*dt)
    return np.array(times), np.array(probs)

print(f"Running 3 simulations ({n_steps} steps each)...")
print("  [1/3] Coherent...");   t_coh, p_coh = evolve(psi0.copy(), V_both)
print("  [2/3] Upper slit..."); t_upp, p_upp = evolve(psi0.copy(), V_upper)
print("  [3/3] Lower slit..."); t_low, p_low = evolve(psi0.copy(), V_lower)
p_incoh = p_upp + p_low

# 干涉项 (纯量子效应!)
interference_term = p_coh - p_incoh  # = 2 Re(ψ₁* ψ₂)
print(f"Done: {len(t_coh)} snapshots")

# 擦除时间点
t_erase = 8.0
erase_idx = np.argmin(np.abs(t_coh - t_erase))
print(f"Erase at t={t_coh[erase_idx]:.1f} (frame {erase_idx}/{len(t_coh)})")

# 屏幕
screen_x = 10.0; screen_idx = np.argmin(np.abs(x - screen_x))
lam = 2*np.pi*hbar/p0; L = screen_x - barrier_x

# ============================================================
# 动画
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(18, 10), facecolor='#0d1117')
gs = fig.add_gridspec(2, 2, height_ratios=[2.5, 1.5], hspace=0.4, wspace=0.3)

ax_main = fig.add_subplot(gs[0, :], facecolor='#0d1117')    # 主 2D 图
ax_scr  = fig.add_subplot(gs[1, 0], facecolor='#0d1117')    # 屏幕图案
ax_diff = fig.add_subplot(gs[1, 1], facecolor='#0d1117')    # 干涉项

gamma = 0.4; vmax = max(p_coh.max(), p_incoh.max())
def to_disp(p): return (np.clip(p/(vmax*0.1), 0, 1))**gamma

# --- 主图 ---
im = ax_main.imshow(to_disp(p_incoh[0]).T, extent=[x[0],x[-1],y[0],y[-1]],
                     origin='lower', cmap='inferno', aspect='equal')
ax_main.axvline(barrier_x, color='white', linewidth=1, alpha=0.25, linestyle='--')
ax_main.axvline(screen_x, color='#f0883e', linewidth=2, alpha=0.5, linestyle='--')
ax_main.text(screen_x+0.2, y[-1]*0.88, 'SCREEN', color='#f0883e', fontsize=9, fontweight='bold')
bbo = ax_main.fill_between([barrier_x+0.3, barrier_x+1.8], y[0], y[-1],
                            color='#d2a8ff', alpha=0.0)
ax_main.text(barrier_x+0.5, y[0]+0.4, 'BBO', color='#d2a8ff', fontsize=8, alpha=0)
mode_text = ax_main.text(0.02, 0.95, '', transform=ax_main.transAxes,
                          color='#e6edf3', fontsize=14, fontweight='bold', va='top')
ax_main.set_ylabel('y', color='#e6edf3', fontsize=11)
ax_main.tick_params(colors='#e6edf3')
cbar = plt.colorbar(im, ax=ax_main, fraction=0.01)
cbar.ax.tick_params(colors='#e6edf3')

# --- 屏幕图案 ---
ax_scr.set_facecolor('#0d1117')
line_s, = ax_scr.plot(y, np.zeros_like(y), color='#79c0ff', linewidth=2.5)
fill_s = ax_scr.fill_between(y, 0, np.zeros_like(y), alpha=0.3, color='#79c0ff')
for n in range(-4, 5):
    ym = n*lam*L/slit_sep
    if abs(ym) < y[-1]:
        ax_scr.axvline(ym, color='#f0883e', linewidth=0.6, linestyle=':', alpha=0.3)
ax_scr.set_xlabel('y (screen)', color='#e6edf3', fontsize=10)
ax_scr.set_ylabel('Intensity', color='#e6edf3', fontsize=10)
ax_scr.set_title('Detector Screen (accumulating)', color='#e6edf3', fontsize=12)
ax_scr.tick_params(colors='#e6edf3'); ax_scr.grid(True, alpha=0.1, color='#30363d')

# --- 干涉项 (纯量子) ---
ax_diff.set_facecolor('#0d1117')
int_max = max(abs(interference_term.max()), abs(interference_term.min()))
im_diff = ax_diff.imshow(interference_term[0].T,
                          extent=[x[0], x[-1], y[0], y[-1]],
                          origin='lower', cmap='RdBu_r', aspect='equal',
                          vmin=-int_max*0.8, vmax=int_max*0.8)
ax_diff.axvline(barrier_x, color='white', linewidth=1, alpha=0.25, linestyle='--')
ax_diff.axvline(screen_x, color='#f0883e', linewidth=1.5, alpha=0.4, linestyle='--')
ax_diff.set_xlabel('x', color='#e6edf3', fontsize=10)
ax_diff.set_ylabel('y', color='#e6edf3', fontsize=10)
ax_diff.set_title(r'Interference Term  $2\,\mathrm{Re}(\psi_1^*\psi_2)$',
                  color='#e6edf3', fontsize=12)
ax_diff.tick_params(colors='#e6edf3')
cbar_d = plt.colorbar(im_diff, ax=ax_diff, fraction=0.03)
cbar_d.ax.tick_params(colors='#e6edf3')

fig_title = fig.suptitle('', color='#e6edf3', fontsize=16, fontweight='bold', y=0.99)

for ax in [ax_main, ax_scr, ax_diff]:
    for spine in ax.spines.values(): spine.set_color('#30363d')

# 累积屏幕图案
accumulated = np.zeros(Ny)

def update(i):
    is_erased = i >= erase_idx

    # 根据阶段选择显示数据
    if is_erased:
        data = p_coh[i]
        mode_text.set_text('◉ QUANTUM ERASER  (interference restored!)')
        mode_text.set_color('#58a6ff')
        bbo.set_alpha(0.08)
    else:
        data = p_incoh[i]
        mode_text.set_text('○ WHICH-PATH  (no interference)')
        mode_text.set_color('#ff7b72')
        bbo.set_alpha(0.0)

    # 主图
    im.set_data(to_disp(data).T)

    # 干涉项
    im_diff.set_data(interference_term[i].T)

    # 屏幕累积 (仅在波包到达屏幕后开始)
    if t_coh[i] > 3.0:
        alpha = 0.3
        accumulated[:] = (1-alpha)*accumulated + alpha*p_coh[i, screen_idx, :]
        line_s.set_ydata(accumulated)
        # 更新 fill
        for coll in ax_scr.collections[:]:
            coll.remove()
        ax_scr.fill_between(y, 0, accumulated, alpha=0.3, color='#79c0ff')

    t = t_coh[i]
    fig_title.set_text(
        f't = {t:.2f}    |    {"ERASER MODE" if is_erased else "WHICH-PATH MODE"}    |    '
        f'λ = {lam:.3f}    |    Δy = {lam*L/slit_sep:.2f}')

    return [im, im_diff, line_s]

# 在切换帧暂停
frames = []
for i in range(len(t_coh)):
    frames.append(i)
    if i == erase_idx:
        for _ in range(30): frames.append(i)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=55, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'quantum_eraser.mp4')

try:
    writer = animation.FFMpegWriter(fps=18, bitrate=3500)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=18), dpi=100)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

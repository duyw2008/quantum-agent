#!/usr/bin/env python3
"""
量子场论 Demo — 自由标量场

展示:
  1. 场对易子 [φ̂(x), φ̂(y)] — 类空区域非零, 微观因果性的量子体现
  2. 真空涨落 ⟨0|φ̂(x)²|0⟩ — 平移不变, 紫外截断依赖
  3. Feynman 传播子 D_F(x-y, t) — 粒子传播的量子描述

3 面板动画: 对易子 vs 距离 | 传播子实部/虚部 vs 距离 | 传播子 vs t
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qft import ScalarField

print("=" * 60)
print("  QFT: Free Scalar Field")
print("=" * 60)

# 创建场
mass = 1.0; L = 20.0; N_modes = 40
sf = ScalarField(mass, L, N_modes)
print(sf.summary())

# 空间点
x_vals = np.linspace(-L/2, L/2, 300)

# ================================================================
# 计算
# ================================================================
comm = sf.commutator_profile(0, x_vals)
prop_t0 = sf.propagator_profile(0, x_vals)  # 等时传播子
vac = sf.vacuum_fluctuation_profile(x_vals)

# 传播子 vs 时间 (在固定 x)
t_vals = np.linspace(0, 5, 100)
prop_vs_t = np.array([sf.feynman_propagator(2.0, 0, t) for t in t_vals])

print(f"\nVacuum fluctuation ⟨φ²⟩ = {sf.vacuum_fluctuation(0):.4f}")
print(f"Vacuum energy density = {sf.vacuum_energy_density():.4f}")
print(f"Commutator at x=1: [φ̂(0),φ̂(1)] = {sf.commutator(0,1):.4f}")
print(f"Propagator at x=1, t=0: D_F = {sf.feynman_propagator(1,0,0):.4f}")

# ================================================================
# 动画
# ================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, (ax_comm, ax_prop, ax_time) = plt.subplots(1, 3, figsize=(18, 5.5),
                                                  facecolor='#0d1117')

# --- 对易子 ---
ax_comm.set_facecolor('#0d1117')
ax_comm.plot(x_vals, comm, color='#58a6ff', linewidth=2)
ax_comm.axhline(0, color='#30363d', linewidth=0.5)
ax_comm.fill_between(x_vals, comm, 0, alpha=0.15, color='#58a6ff')
cursor_c = ax_comm.axvline(0, color='#f0883e', linewidth=2, alpha=0)
dot_c, = ax_comm.plot([0], [sf.commutator(0,0)], 'o', color='#f0883e', markersize=8, alpha=0)
ax_comm.set_xlabel('x - y', color='#e6edf3', fontsize=11)
ax_comm.set_ylabel(r'$\langle 0|[\hat{\phi}(0),\hat{\phi}(x)]|0\rangle$', color='#e6edf3', fontsize=11)
ax_comm.set_title('Field Commutator', color='#e6edf3', fontsize=13, fontweight='bold')
ax_comm.tick_params(colors='#e6edf3'); ax_comm.grid(True, alpha=0.1, color='#30363d')

# --- 传播子 ---
ax_prop.set_facecolor('#0d1117')
line_p_real, = ax_prop.plot(x_vals, np.real(prop_t0), color='#79c0ff', linewidth=2, label='Re D_F')
line_p_imag, = ax_prop.plot(x_vals, np.imag(prop_t0), color='#ff7b72', linewidth=1.5, linestyle='--', label='Im D_F')
ax_prop.axhline(0, color='#30363d', linewidth=0.5)
cursor_p = ax_prop.axvline(0, color='#f0883e', linewidth=2, alpha=0)
ax_prop.set_xlabel('x - y', color='#e6edf3', fontsize=11)
ax_prop.set_ylabel(r'$D_F(x-y)$', color='#e6edf3', fontsize=11)
ax_prop.set_title('Feynman Propagator (t=0)', color='#e6edf3', fontsize=13, fontweight='bold')
ax_prop.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3')
ax_prop.tick_params(colors='#e6edf3'); ax_prop.grid(True, alpha=0.1, color='#30363d')

# --- 传播子 vs t ---
ax_time.set_facecolor('#0d1117')
line_t_real, = ax_time.plot([], [], color='#79c0ff', linewidth=2)
line_t_imag, = ax_time.plot([], [], color='#ff7b72', linewidth=1.5, linestyle='--')
ax_time.set_xlim(0, t_vals[-1])
ax_time.set_ylim(-0.2, 0.3)
ax_time.axhline(0, color='#30363d', linewidth=0.5)
cursor_t = ax_time.axvline(0, color='#f0883e', linewidth=2, alpha=0)
ax_time.set_xlabel('t', color='#e6edf3', fontsize=11)
ax_time.set_ylabel(r'$D_F(\Delta x=2, t)$', color='#e6edf3', fontsize=11)
ax_time.set_title('Propagator vs Time', color='#e6edf3', fontsize=13, fontweight='bold')
ax_time.tick_params(colors='#e6edf3'); ax_time.grid(True, alpha=0.1, color='#30363d')

fig_title = fig.suptitle('', color='#e6edf3', fontsize=15, fontweight='bold', y=0.99)

for ax in [ax_comm, ax_prop, ax_time]:
    for spine in ax.spines.values(): spine.set_color('#30363d')

frames = len(t_vals)
extended = []
for i in range(frames):
    extended.append(i)
    if i == 0:
        for _ in range(20): extended.append(0)

def update(i):
    t = t_vals[i]
    # 传播子 vs t
    prop_now = np.array([sf.feynman_propagator(2.0, 0, tt) for tt in t_vals[:i+1]])
    line_t_real.set_data(t_vals[:i+1], np.real(prop_now))
    line_t_imag.set_data(t_vals[:i+1], np.imag(prop_now))

    # 在 x=2 处标记
    cursor_c.set_xdata([2, 2])
    cursor_c.set_alpha(0.5)
    dot_c.set_data([2], [sf.commutator(0, 2)])
    dot_c.set_alpha(1)
    cursor_p.set_xdata([2, 2])
    cursor_p.set_alpha(0.5)
    cursor_t.set_xdata([t, t])
    cursor_t.set_alpha(0.5)

    fig_title.set_text(
        f'Free Scalar Field  |  m={mass}  |  t = {t:.2f}  |  '
        f'[φ̂(0),φ̂(2)] = {sf.commutator(0,2):.4f}')

    return [line_t_real, line_t_imag, cursor_c, cursor_t]

ani = animation.FuncAnimation(fig, update, frames=extended,
                               interval=60, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'qft_scalar_field.mp4')

try:
    writer = animation.FFMpegWriter(fps=18, bitrate=2800)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=18), dpi=100)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

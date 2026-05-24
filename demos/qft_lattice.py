#!/usr/bin/env python3
"""
格点 φ⁴ 理论 — 相互作用量子场论 Demo

展示:
  1. 基态能量 E₀ vs 耦合常数 λ
  2. 关联函数 ⟨φ₀ φ_d⟩ vs 距离 d
  3. 能隙 Δ = E₁ - E₀ vs λ
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qft import LatticePhi4

print("=" * 60)
print("  Lattice φ⁴ Theory")
print("=" * 60)

# 创建格点
N_sites, mass, N_fock = 5, 0.5, 3
lpt = LatticePhi4(N_sites, mass, 0.0, N_fock)
print(lpt.summary())
print(f"  Hilbert dim = {N_fock}^{N_sites} = {lpt.dim}")

# 扫描耦合常数
couplings = np.linspace(0, 3.0, 40)
print(f"\nScanning λ ∈ [0, {couplings[-1]}] ...")
result = lpt.scan_coupling(couplings)
print(f"  E₀(0) = {result['E0'][0]:.4f},  E₀({couplings[-1]}) = {result['E0'][-1]:.4f}")

# 粒子数分布
for lam in [0.0, 1.0, 2.0]:
    n_j = lpt.particle_number_distribution(lam)
    print(f"  λ={lam}: ⟨Nⱼ⟩ = {n_j.round(3)}")

# 关联函数
print(f"\nCorrelation ⟨φ₀ φ_d⟩:")
for lam in [0.0, 1.0, 3.0]:
    corr = lpt.correlation_function(lam)
    print(f"  λ={lam}: {corr.round(3)}")

# ================================================================
# 动画
# ================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, (ax_e, ax_c, ax_g) = plt.subplots(1, 3, figsize=(18, 5.5),
                                        facecolor='#0d1117')

# --- 基态能量 ---
ax_e.set_facecolor('#0d1117')
ax_e.set_xlim(0, couplings[-1])
E_all = result['E0']
e_margin = (E_all.max() - E_all.min()) * 0.1
ax_e.set_ylim(E_all.min() - e_margin, E_all.max() + e_margin)
line_e, = ax_e.plot([], [], color='#58a6ff', linewidth=2.5)
dot_e, = ax_e.plot([], [], 'o', color='#f0883e', markersize=8)
ax_e.set_xlabel('λ (coupling)', color='#e6edf3', fontsize=11)
ax_e.set_ylabel('E₀', color='#e6edf3', fontsize=11)
ax_e.set_title('Ground State Energy', color='#e6edf3', fontsize=13, fontweight='bold')
ax_e.tick_params(colors='#e6edf3'); ax_e.grid(True, alpha=0.1, color='#30363d')

# --- 关联函数 ---
ax_c.set_facecolor('#0d1117')
d_vals = np.arange(N_sites)
colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(couplings)))
line_corr = [ax_c.plot([], [], 'o-', color=colors[i], markersize=5, alpha=0)[0]
             for i in range(len(couplings))]
ax_c.set_xlabel('distance d', color='#e6edf3', fontsize=11)
ax_c.set_ylabel(r'$\langle\phi_0 \phi_d\rangle$', color='#e6edf3', fontsize=11)
ax_c.set_title('Correlation Function', color='#e6edf3', fontsize=13, fontweight='bold')
ax_c.tick_params(colors='#e6edf3'); ax_c.grid(True, alpha=0.1, color='#30363d')
all_corr = result['correlation']
c_max = max(abs(all_corr.min()), abs(all_corr.max())) * 1.1
ax_c.set_ylim(-c_max, c_max)
ax_c.set_xlim(-0.5, N_sites - 0.5)

# --- 能隙 ---
ax_g.set_facecolor('#0d1117')
ax_g.set_xlim(0, couplings[-1])
gap_all = result['gap']
g_margin = (gap_all.max() - gap_all.min()) * 0.1
ax_g.set_ylim(max(0, gap_all.min() - g_margin), gap_all.max() + g_margin)
line_g, = ax_g.plot([], [], color='#d2a8ff', linewidth=2.5)
dot_g, = ax_g.plot([], [], 'o', color='#f0883e', markersize=8)
ax_g.set_xlabel('λ (coupling)', color='#e6edf3', fontsize=11)
ax_g.set_ylabel('Δ = E₁ - E₀', color='#e6edf3', fontsize=11)
ax_g.set_title('Energy Gap', color='#e6edf3', fontsize=13, fontweight='bold')
ax_g.tick_params(colors='#e6edf3'); ax_g.grid(True, alpha=0.1, color='#30363d')

fig_title = fig.suptitle('', color='#e6edf3', fontsize=15, fontweight='bold', y=0.99)

for ax in [ax_e, ax_c, ax_g]:
    for spine in ax.spines.values(): spine.set_color('#30363d')

frames = len(couplings)
extended = list(range(frames)) + [frames-1]*15

def update(i):
    lam = couplings[i]
    line_e.set_data(couplings[:i+1], E_all[:i+1])
    dot_e.set_data([lam], [E_all[i]])

    # 显示当前 λ 的关联函数
    for j, lc in enumerate(line_corr):
        lc.set_alpha(0.15 if j != i else 1.0)
        lc.set_data(d_vals, all_corr[j])
        lc.set_markersize(3 if j != i else 7)

    line_g.set_data(couplings[:i+1], gap_all[:i+1])
    dot_g.set_data([lam], [gap_all[i]])

    fig_title.set_text(
        f'Lattice φ⁴  |  N={N_sites}, m={mass}  |  '
        f'λ={lam:.2f}  |  E₀={E_all[i]:.4f}  |  Δ={gap_all[i]:.4f}')

    return [line_e, dot_e, line_g, dot_g] + line_corr

ani = animation.FuncAnimation(fig, update, frames=extended,
                               interval=80, blit=False)

save_dir = os.path.join(os.path.dirname(__file__), 'output', 'animations')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'qft_lattice_phi4.mp4')

try:
    writer = animation.FFMpegWriter(fps=15, bitrate=2500)
    ani.save(save_path, writer=writer, dpi=150)
except (FileNotFoundError, RuntimeError):
    save_path = save_path.replace('.mp4', '.gif')
    ani.save(save_path, writer=animation.PillowWriter(fps=15), dpi=100)

plt.close(fig)
print(f"\nAnimation saved: {save_path}")

#!/usr/bin/env python3
"""氢原子 1D 库仑势 — 电子波函数图像生成"""

import sys, os
sys.path.insert(0, '/home/duyw/quantum_agent')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

from src.core import Grid, create_potential, PotentialType
from src.viz.animate import QM_DARK_THEME

theme = QM_DARK_THEME

# ============================================================
# 参数
# ============================================================
Z = 1.0
softening = 0.1
hbar = 1.0
mass = 1.0
x_range = (-30, 30)
n_points = 2048

grid = Grid(x_range[0], x_range[1], n_points)
V = create_potential(PotentialType.COULOMB_1D, Z=Z, softening=softening)
x, dx = grid.x, grid.dx

print(f"Grid: {grid}")
print(f"Potential: {V.label}")
print(f"Soften parameter a = {softening}")

# ============================================================
# 哈密顿量对角化
# ============================================================
N = grid.n_points
ke_diag = np.ones(N) * hbar**2 / (mass * dx**2)
ke_off = np.ones(N - 1) * (-hbar**2 / (2 * mass * dx**2))
Vx = V(x)
diag = ke_diag + Vx

n_states = 12
eigvals, eigvecs = eigh_tridiagonal(diag, ke_off,
                                     select='i', select_range=(0, n_states - 1))

print(f"\nEigenvalues:")
for n, E in enumerate(eigvals):
    n_qn = n + 1
    E_3d = -Z**2 / (2 * n_qn**2)
    print(f"  n={n_qn:2d}: E_1D={E:+.6f}  (3D: {E_3d:+.6f})")

bound_states = sum(1 for E in eigvals if E < 0)
print(f"\nBound states: {bound_states}")

# ============================================================
# 图 1: 全景 (势函数 + 波函数 + 能级)
# ============================================================
fig, (ax_pot, ax_wf, ax_en) = plt.subplots(1, 3, figsize=(20, 8),
                                             facecolor=theme['bg'],
                                             gridspec_kw={'width_ratios': [1.2, 2.5, 0.8]})

# Left: potential + levels
ax_pot.set_facecolor(theme['bg'])
Vx_clipped = np.clip(Vx, -15, 5)
ax_pot.fill_between(x, Vx_clipped, Vx_clipped.min(), alpha=0.3, color=theme['potential'])
ax_pot.plot(x, Vx_clipped, color=theme['potential'], linewidth=2)

for n, E in enumerate(eigvals[:bound_states]):
    ax_pot.axhline(y=E, color=theme['energy'], linewidth=1.5, linestyle='--', alpha=0.8)
    ax_pot.text(x_range[1] * 0.85, E, f'n={n+1}  {E:+.3f}',
                color=theme['fg'], va='center', fontsize=9, fontfamily='monospace')

ax_pot.set_xlabel('x (a.u.)', color=theme['fg'], fontsize=12)
ax_pot.set_ylabel('V(x) / Energy', color=theme['fg'], fontsize=12)
ax_pot.set_title(f'Hydrogen 1D Coulomb\nZ={Z}, a={softening}', color=theme['fg'], fontsize=12)
ax_pot.set_ylim(-15, 2)
ax_pot.tick_params(colors=theme['fg'])
ax_pot.grid(True, alpha=0.15, color=theme['grid'])

# Middle: wavefunctions
ax_wf.set_facecolor(theme['bg'])
show_states = min(bound_states, 8)
colors_wf = plt.cm.plasma(np.linspace(0.15, 0.9, show_states))

for n in range(show_states):
    E = eigvals[n]
    psi_n = eigvecs[:, n]
    prob_n = psi_n**2
    scale = 1.5 if n == 0 else 1.0 / (n + 0.5)
    prob_scaled = prob_n / prob_n.max() * scale * 2
    ax_wf.fill_between(x, E, E + prob_scaled, alpha=0.5, color=colors_wf[n], label=f'n={n+1}')
    ax_wf.plot(x, E + psi_n * scale * np.sqrt(prob_n.max()),
               color=colors_wf[n], linewidth=0.8, alpha=0.7)

ax_wf.set_xlabel('x (a.u.)', color=theme['fg'], fontsize=12)
ax_wf.set_ylabel('Energy / |ψ|²', color=theme['fg'], fontsize=12)
ax_wf.set_title('Eigenstate Wavefunctions', color=theme['fg'], fontsize=14)
ax_wf.tick_params(colors=theme['fg'])
ax_wf.grid(True, alpha=0.15, color=theme['grid'])
ax_wf.legend(loc='upper right', facecolor=theme['bg'],
             edgecolor=theme['grid'], labelcolor=theme['fg'], fontsize=9, ncol=2)

# Right: energy levels
ax_en.set_facecolor(theme['bg'])
for n in range(show_states):
    E = eigvals[n]
    ax_en.hlines(y=E, xmin=0.1, xmax=0.9, color=colors_wf[n], linewidth=3)
    ax_en.text(0.95, E, f'  {n+1}', color=theme['fg'],
               va='center', fontsize=11, fontfamily='monospace')

ax_en.axhline(y=0, color=theme['grid'], linewidth=1, linestyle='-')
ax_en.text(0.5, 0.15, 'continuum', color=theme['grid'],
           ha='center', fontsize=10, fontstyle='italic')
ax_en.set_xlim(0, 1.5)
ax_en.set_ylim(ax_pot.get_ylim())
ax_en.set_title('Levels', color=theme['fg'], fontsize=12)
ax_en.set_xticks([])
ax_en.tick_params(colors=theme['fg'])
ax_en.grid(True, axis='y', alpha=0.15, color=theme['grid'])

for ax in [ax_pot, ax_wf, ax_en]:
    for spine in ax.spines.values():
        spine.set_color(theme['grid'])

plt.tight_layout()
save_dir = '/home/duyw/quantum_agent/output/figures'
os.makedirs(save_dir, exist_ok=True)
fig.savefig(os.path.join(save_dir, 'hydrogen_coulomb_1d.png'),
            dpi=200, bbox_inches='tight', facecolor=theme['bg'])
plt.close(fig)
print(f"✓ Saved: {save_dir}/hydrogen_coulomb_1d.png")

# ============================================================
# 图 2: 概率密度对比
# ============================================================
fig2, ax2 = plt.subplots(figsize=(14, 6), facecolor=theme['bg'])
ax2.set_facecolor(theme['bg'])

for n in range(min(bound_states, 6)):
    psi_n = eigvecs[:, n]
    prob_n = psi_n**2
    E = eigvals[n]
    offset = n * 0.8
    prob_norm = prob_n / prob_n.max()
    ax2.fill_between(x, offset, offset + prob_norm, alpha=0.5, color=colors_wf[n])
    ax2.plot(x, offset + prob_norm, color=colors_wf[n], linewidth=1.2)
    nodes = n
    ax2.text(x_range[1] * 0.92, offset + 0.5,
             f'|{n+1}⟩  E={E:+.3f}  nodes={nodes}',
             color=theme['fg'], va='center', fontsize=10, fontfamily='monospace')

ax2.set_xlabel('x (a.u.)', color=theme['fg'], fontsize=12)
ax2.set_ylabel('Probability Density |ψ(x)|²', color=theme['fg'], fontsize=12)
ax2.set_title(f'Hydrogen 1D Coulomb — Electron Probability Densities (Z={Z})',
              color=theme['fg'], fontsize=14)
ax2.set_yticks([])
ax2.tick_params(colors=theme['fg'])
ax2.grid(True, alpha=0.15, color=theme['grid'])
for spine in ax2.spines.values():
    spine.set_color(theme['grid'])

plt.tight_layout()
fig2.savefig(os.path.join(save_dir, 'hydrogen_density.png'),
             dpi=200, bbox_inches='tight', facecolor=theme['bg'])
plt.close(fig2)
print(f"✓ Saved: {save_dir}/hydrogen_density.png")

# ============================================================
# 图 3: 径向概率分布风格 (|ψ|² 在正半轴)
# ============================================================
fig3, ax3 = plt.subplots(figsize=(12, 6), facecolor=theme['bg'])
ax3.set_facecolor(theme['bg'])

# 只看 x ≥ 0
pos_mask = x >= 0
x_pos = x[pos_mask]
for n in range(min(bound_states, 6)):
    psi_n = eigvecs[:, n]
    prob_n = psi_n**2
    prob_pos = prob_n[pos_mask]
    r2_factor = x_pos**2  # 模拟 r² 因子 (3D 径向分布: P(r) = r²|R(r)|²)
    radial_prob = prob_pos * r2_factor
    radial_prob /= radial_prob.max()
    ax3.plot(x_pos, radial_prob, color=colors_wf[n], linewidth=2, label=f'n={n+1}')
    ax3.fill_between(x_pos, 0, radial_prob, alpha=0.25, color=colors_wf[n])

ax3.set_xlabel('r (a.u.)  [x ≥ 0 half-space]', color=theme['fg'], fontsize=12)
ax3.set_ylabel(r'$x^2 |\psi(x)|^2$  (radial distribution analog)', color=theme['fg'], fontsize=12)
ax3.set_title(f'Hydrogen 1D — Radial Distribution (× x² weight)',
              color=theme['fg'], fontsize=14)
ax3.legend(loc='upper right', facecolor=theme['bg'],
           edgecolor=theme['grid'], labelcolor=theme['fg'])
ax3.tick_params(colors=theme['fg'])
ax3.grid(True, alpha=0.15, color=theme['grid'])
ax3.set_xlim(0, 15)
for spine in ax3.spines.values():
    spine.set_color(theme['grid'])

plt.tight_layout()
fig3.savefig(os.path.join(save_dir, 'hydrogen_radial.png'),
             dpi=200, bbox_inches='tight', facecolor=theme['bg'])
plt.close(fig3)
print(f"✓ Saved: {save_dir}/hydrogen_radial.png")

print(f"\n{'='*60}")
print(f"  Hydrogen 1D Coulomb — Image Generation Complete")
print(f"{'='*60}")
print(f"  3 images saved to {save_dir}/")
print(f"    - hydrogen_coulomb_1d.png   (full overview)")
print(f"    - hydrogen_density.png      (probability densities)")
print(f"    - hydrogen_radial.png       (radial distribution analog)")
print(f"{'='*60}")

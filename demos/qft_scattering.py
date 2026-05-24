#!/usr/bin/env python3
"""
散射与 Feynman 图 Demo

展示:
  1. Wick 定理 — 场算符的所有收缩方式
  2. Feynman 图 — φ⁴ 2→2 树图 + 截面
  3. 跃迁概率 — Dyson 级数逐阶计算
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qft import (
    wick_expand, feynman_amplitude_phi4_2to2,
    differential_cross_section, draw_feynman_phi4_2to2,
    LatticePhi4,
)

print("=" * 60)
print("  QFT Scattering & Feynman Diagrams")
print("=" * 60)

# ================================================================
# 1. Wick 定理
# ================================================================
print("\n--- 1. Wick's Theorem ---")

phi_list = ['φ₁', 'φ₂', 'φ₃', 'φ₄']
contractions = wick_expand(phi_list)
print(f"Field operators: {' '.join(phi_list)}")
print(f"All contractions ({len(contractions)} ways):")
for i, contr in enumerate(contractions):
    pairs = ' + '.join(f'{phi_list[a]}{phi_list[b]}' for a, b in contr)
    print(f"  {i+1}. {pairs}")

# ================================================================
# 2. Feynman 图
# ================================================================
print("\n--- 2. Feynman Diagram ---")
print(draw_feynman_phi4_2to2())

lam = 0.5
amp = feynman_amplitude_phi4_2to2(lam)
print(f"  λ = {lam}")
print(f"  iM = {amp:.3f}  (= -iλ)")
print(f"  |M|² = {abs(amp)**2:.3f}")

# 截面 vs 质心能量
print(f"\n  Cross section vs √s:")
for s in [5, 10, 20, 50, 100]:
    dsigma = differential_cross_section(lam, s)
    print(f"    √s={np.sqrt(s):5.1f}  →  dσ/dΩ = {dsigma:.6f}")

# ================================================================
# 3. 跃迁概率 (用格点 φ⁴)
# ================================================================
print("\n--- 3. Transition Probability ---")

from src.qft.scattering import transition_probability

# 用格点模型计算 |⟨f|U(t)|i⟩|²
N_sites = 3
lpt = LatticePhi4(N_sites, mass=1.0, coupling=lam, N_fock=3)

# |i⟩ = |1,0,0⟩ (一个粒子在格点 0)
# |f⟩ = |0,1,0⟩ (粒子跃迁到格点 1)
initial = np.zeros(lpt.dim); initial[0] = 1.0  # 简化: 用基矢
# 构建简单的初态和末态
# 在 3 格点系统中: |n₀,n₁,n₂⟩ 索引 = n₀*9 + n₁*3 + n₂
i_idx = 1 * 9 + 0 * 3 + 0  # |1,0,0⟩
f_idx = 0 * 9 + 1 * 3 + 0  # |0,1,0⟩
initial = np.zeros(lpt.dim); initial[i_idx] = 1.0
final = np.zeros(lpt.dim); final[f_idx] = 1.0

H0 = lpt.hamiltonian(0.0)    # 自由理论
V  = lpt.hamiltonian(lam) - H0  # 相互作用

print(f"  N_sites={N_sites}, λ={lam}")
print(f"  |i⟩ = |1,0,0⟩, |f⟩ = |0,1,0⟩")
print(f"  Computing transition probability...")
probs = transition_probability(H0, V, initial, final, t=1.0, max_order=2)

for order, prob in sorted(probs.items()):
    bar = '█' * int(prob * 50) if prob > 0 else ''
    print(f"  Order {order}: P = {prob:.6f}  {bar}")

# ================================================================
# 生成图像
# ================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, (ax_w, ax_s) = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='#0d1117')

# --- Wick 收缩示意图 ---
ax_w.set_facecolor('#0d1117')
ax_w.set_xlim(0, 3); ax_w.set_ylim(0, 3)
ax_w.axis('off')
ax_w.set_title("Wick Contractions: ⟨φ₁φ₂φ₃φ₄⟩", color='#e6edf3', fontsize=14, fontweight='bold')

# 画 4 个点 (场算符)
positions = [(0.5, 1.5), (0.5, 0.5), (2.5, 1.5), (2.5, 0.5)]
labels = ['φ₁', 'φ₂', 'φ₃', 'φ₄']
colors_pt = ['#58a6ff', '#58a6ff', '#ff7b72', '#ff7b72']

for (px, py), label, c in zip(positions, labels, colors_pt):
    ax_w.plot(px, py, 'o', color=c, markersize=25, markeredgecolor='white', markeredgewidth=1.5)
    ax_w.text(px, py, label, color='white', fontsize=11, ha='center', va='center', fontweight='bold')

# 画收缩线
contr = contractions[0]  # 第一种收缩方式
line_styles = ['-', '--']
colors_line = ['#79c0ff', '#d2a8ff']
for idx, (a, b) in enumerate(contr):
    x1, y1 = positions[a]
    x2, y2 = positions[b]
    # 曲线连接
    mid_x = (x1 + x2) / 2
    mid_y = max(y1, y2) + 0.3 + idx * 0.2
    curve = plt.matplotlib.patches.FancyArrowPatch(
        (x1, y1), (x2, y2),
        connectionstyle=f"arc3,rad={0.3 + idx*0.2}",
        color=colors_line[idx], linewidth=2.5, alpha=0.8
    )
    ax_w.add_patch(curve)

# --- 截面图 ---
ax_s.set_facecolor('#0d1117')
s_vals = np.linspace(5, 200, 100)
sigma_vals = [differential_cross_section(lam, s) for s in s_vals]
ax_s.plot(np.sqrt(s_vals), sigma_vals, color='#58a6ff', linewidth=2.5)
ax_s.fill_between(np.sqrt(s_vals), 0, sigma_vals, alpha=0.2, color='#58a6ff')
ax_s.set_xlabel(r'$\sqrt{s}$ (CM energy)', color='#e6edf3', fontsize=11)
ax_s.set_ylabel(r'$d\sigma/d\Omega$', color='#e6edf3', fontsize=11)
ax_s.set_title(f'φ⁴ 2→2 Cross Section (λ={lam})', color='#e6edf3', fontsize=13, fontweight='bold')
ax_s.tick_params(colors='#e6edf3'); ax_s.grid(True, alpha=0.1, color='#30363d')

for ax in [ax_s]:
    for spine in ax.spines.values(): spine.set_color('#30363d')

plt.tight_layout()
save_dir = os.path.join(os.path.dirname(__file__), 'output', 'figures')
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'qft_scattering.png')
fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
plt.close(fig)

print(f"\nFigure saved: {save_path}")
print(f"\n{'='*60}")
print(f"  Scattering Demo Complete")
print(f"{'='*60}")

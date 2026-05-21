#!/usr/bin/env python3
"""
Demo 3: 势垒隧穿与散射

展示一个高斯波包遇到势垒时的量子隧穿效应。
比较经典禁区 (E < V₀) 和经典允许区 (E > V₀) 的行为。
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import animate_evolution

print("=" * 60)
print("  Demo 3: Quantum Tunneling through a Barrier")
print("=" * 60)

hbar = 1.0
mass = 1.0

# 势垒参数
barrier_height = 10.0
barrier_width = 0.5
V = create_potential(PotentialType.POTENTIAL_BARRIER,
                     height=barrier_height, width=barrier_width)
print(f"\nBarrier: V₀ = {barrier_height}, width = {barrier_width}")

# 网格 (宽区域以适应散射)
grid = Grid(-8.0, 8.0, 2048)
print(f"Grid: {grid}")

# ============================================================
# 初始波包: 从左侧入射
# ============================================================
x0 = -3.0  # 在势垒左侧
p0 = 4.0   # 动量 (动能 E = p²/2m = 8.0 < V₀ = 10.0 → 经典禁区!)
sigma = 0.6

wf = WaveFunction(grid)
wf.set_gaussian(x0=x0, p0=p0, sigma=sigma)

E_kinetic = p0**2 / (2 * mass)
print(f"\nIncident wavepacket:")
print(f"  Kinetic energy E = {E_kinetic:.4f}")
print(f"  Barrier V₀ = {barrier_height}")
print(f"  E < V₀: {'YES — quantum tunneling regime!' if E_kinetic < barrier_height else 'NO — classical regime'}")
print(f"  ⟨x⟩ = {wf.expectation_x():.4f}, ⟨p⟩ = {wf.expectation_p():.4f}")
print(f"  Δx·Δp = {wf.uncertainty_product():.4f}")

# 隧穿概率的 WKB 近似
kappa = np.sqrt(2 * mass * (barrier_height - E_kinetic)) / hbar
T_wkb = np.exp(-2 * kappa * barrier_width)
print(f"  WKB tunneling probability ≈ exp(-2κa) = {T_wkb:.4e}")
print(f"  (κ = {kappa:.4f})")

# ============================================================
# 演化
# ============================================================
t_max = 5.0
dt = 0.002
print(f"\nEvolving for t_max = {t_max}, dt = {dt}")

solver = create_solver('ssfm', grid, V, hbar, mass)
result = solver.evolve(wf, t_max, dt, snapshot_interval=10)

print(f"\n✓ Evolution complete: {len(result.times)} snapshots")

# 分析透射和反射
# 在势垒右侧 (x > barrier_width/2) 的概率
prob = result.prob_density_snapshots
x = grid.x
barrier_edge = barrier_width / 2 + 0.5  # 势垒右侧一点
right_mask = x > barrier_edge
left_mask = x < -barrier_edge

transmitted = np.array([np.trapz(p[right_mask], x[right_mask]) for p in prob])
reflected = np.array([np.trapz(p[left_mask], x[left_mask]) for p in prob])

print(f"  Final transmitted probability: {transmitted[-1]:.4f}")
print(f"  Final reflected probability:  {reflected[-1]:.4f}")
print(f"  Total: {transmitted[-1] + reflected[-1]:.4f}")

# 检查 x=0 附近 (透射+反射 ≠ 1 因为波包还在势垒区域)
barrier_prob = np.array([np.trapz(p[(x > -barrier_edge) & (x < barrier_edge)],
                                   x[(x > -barrier_edge) & (x < barrier_edge)])
                          for p in prob])
print(f"  In barrier region: {barrier_prob[-1]:.4f}")

# ============================================================
# 动画
# ============================================================
output_dir = os.environ.get('output_dir',
    os.path.join(os.path.dirname(__file__), '..', 'output'))

anim_path = os.path.join(output_dir, 'animations', 'quantum_tunneling.mp4')
print(f"\nGenerating animation → {anim_path}")
animate_evolution(result, V, save_path=anim_path, fps=30, dark=True)
print(f"  ✓ Animation saved!")

print(f"\n{'='*60}")
print(f"  Demo complete!")
print(f"  Quantum tunneling demonstrated: some probability leaks through")
print(f"  the barrier even though the classical particle would be reflected.")
print(f"{'='*60}")

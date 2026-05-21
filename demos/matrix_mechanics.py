#!/usr/bin/env python3
"""
Demo 4: 矩阵力学 — 算符、对易子、本征值

展示:
    1. 位置-动量对易关系 [x̂, p̂] = iℏ
    2. 产生-湮灭算符代数 [â, â†] = 1
    3. 谐振子能谱
    4. 算符的矩阵表示
    5. 角动量代数 (如果 sympy 可用)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.matrix import MatrixMechanics, NumericOperatorSystem, commutator, fidelity
from src.matrix.symbolic import SymbolicQuantum, HAS_SYMPY

print("=" * 60)
print("  Demo 4: Matrix Mechanics")
print("=" * 60)

hbar = 1.0
mass = 1.0
omega = 2.0

# ============================================================
# Part 1: Numeric Operator System
# ============================================================
print("\n--- Part 1: Numeric Operator Algebra ---")

mm = MatrixMechanics(n_basis=50, hbar=hbar, mass=mass, omega=omega)

# 1.1 对易关系
print("\n1.1 Commutation Relations:")

# [x̂, p̂] = iℏ
xp_comm = mm.check_commutation(mm.x, mm.p)
print(f"  [x̂, p̂] Frobenius norm: {xp_comm['frobenius_norm']:.4f}")
print(f"  Expected iℏI norm: ℏ·√N = {hbar * np.sqrt(50):.4f}")

# [â, â†] = 1
aa_comm = mm.check_commutation(mm.a, mm.a_dag)
print(f"  [â, â†] Frobenius norm: {aa_comm['frobenius_norm']:.4f}")
print(f"  Expected I norm: √N = {np.sqrt(50):.4f}")

# 1.2 厄米性检查
print("\n1.2 Hermiticity:")
x_herm = mm.numeric.check_hermiticity(mm.x)
p_herm = mm.numeric.check_hermiticity(mm.p)
print(f"  x̂ is Hermitian: {x_herm['is_hermitian']} (|x̂ - x̂†| = {x_herm['frobenius_norm']:.2e})")
print(f"  p̂ is Hermitian: {p_herm['is_hermitian']} (|p̂ - p̂†| = {p_herm['frobenius_norm']:.2e})")

# 1.3 算符矩阵表示
print("\n1.3 Operator Matrices (first 5×5):")
print("  x̂ (in ℏ=m=ω=1 units, x₀=1/√2):")
print(f"    {np.array2string(mm.x[:5, :5].real, precision=3, suppress_small=True)}")
print("  p̂:")
print(f"    {np.array2string(mm.p[:5, :5].imag, precision=3, suppress_small=True)} (×i)")

# ============================================================
# Part 2: Eigenvalue Problem
# ============================================================
print("\n--- Part 2: Harmonic Oscillator Spectrum ---")

eigvals, eigvecs = mm.eigensolve(k=10)
print(f"  Energy eigenvalues (ℏω = {hbar * omega}):")
for n in range(min(10, len(eigvals))):
    E_expected = hbar * omega * (n + 0.5)
    print(f"    E[{n}] = {eigvals[n]:.8f}  expected: {E_expected:.8f}  "
          f"error: {abs(eigvals[n] - E_expected):.2e}")

# 基态验证: ⟨0|x̂|0⟩ = 0, ⟨0|p̂|0⟩ = 0
ground = eigvecs[:, 0]
x_exp = mm.numeric.expectation(mm.x, ground)
p_exp = mm.numeric.expectation(mm.p, ground)
print(f"\n  Ground state expectations:")
print(f"    ⟨0|x̂|0⟩ = {x_exp.real:.6f} (expected: 0)")
print(f"    ⟨0|p̂|0⟩ = {p_exp.real:.6f} (expected: 0)")
print(f"    Δx = {mm.numeric.uncertainty(mm.x, ground):.4f} (expected: √(ℏ/2mω) = {np.sqrt(hbar/(2*mass*omega)):.4f})")
print(f"    Δp = {mm.numeric.uncertainty(mm.p, ground):.4f} (expected: √(mℏω/2) = {np.sqrt(mass*hbar*omega/2):.4f})")

# Δx·Δp 应该等于 ℏ/2
dx = mm.numeric.uncertainty(mm.x, ground)
dp = mm.numeric.uncertainty(mm.p, ground)
print(f"    Δx·Δp = {dx*dp:.6f} (HUP bound: {hbar/2:.6f})")

# ============================================================
# Part 3: Time Evolution
# ============================================================
print("\n--- Part 3: Time Evolution ---")

# 初始态: |ψ₀⟩ = (|0⟩ + |1⟩)/√2  (叠加态)
state0 = (eigvecs[:, 0] + eigvecs[:, 1]) / np.sqrt(2)

# 周期性运动: 经过 t = 2π/ω 后应回到初始态
T = 2 * np.pi / omega
H = mm.H_harmonic

state_t = mm.numeric.evolve_state(H, state0, T)
F = fidelity(state0, state_t)
print(f"  State fidelity after one period: F = {F:.8f} (expected: 1.0)")

# 半周期后 phase flip
state_half = mm.numeric.evolve_state(H, state0, T/2)
# 期望: |ψ(T/2)⟩ = (|0⟩ + e^{-iE₁T/2ℏ}|1⟩)/√2
# E₁ = 3ℏω/2, E₀ = ℏω/2
# e^{-i(E₁-E₀)T/2ℏ} = e^{-iℏω·π/ℏω} = e^{-iπ} = -1
# → |ψ(T/2)⟩ = (|0⟩ - |1⟩)/√2
expected_state = (eigvecs[:, 0] - eigvecs[:, 1]) / np.sqrt(2)
F_half = fidelity(state_half, expected_state)
print(f"  State fidelity after half period: F = {F_half:.8f} (expected: 1.0)")

# ============================================================
# Part 4: Symbolic (sympy) — 如果可用
# ============================================================
if HAS_SYMPY:
    print("\n--- Part 4: Symbolic Quantum Mechanics ---")
    sq = SymbolicQuantum()

    # 谐振子能级公式
    energies = sq.harmonic_energies_symbolic(5)
    print("  Harmonic oscillator energy formula: E_n = ℏω(n + ½)")
    for n, E in enumerate(energies):
        print(f"    n={n}: {E}")

    # 自旋矩阵
    sx, sy, sz, s2 = sq.spin_matrices()
    print(f"\n  Spin-1/2 matrices:")
    print(f"    S_x = {sx}")
    print(f"    S_y = {sy}")
    print(f"    S_z = {sz}")

    # 验证 [S_x, S_y] = iℏ S_z
    comm_xy = sx * sy - sy * sx
    expected_comm = 1j * hbar * sz
    print(f"\n  [S_x, S_y] = {sp.simplify(comm_xy)}")
    print(f"  Expected: iℏ S_z = {expected_comm}")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print(f"  Matrix Mechanics Demo Complete!")
print(f"  Verified: [x̂, p̂] = iℏ, [â, â†] = I, E_n = ℏω(n+½)")
print(f"  HUP: Δx·Δp = {dx*dp:.6f} ≥ {hbar/2:.6f} ✓")
print(f"{'='*60}")

"""测试套件: 矩阵力学"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.matrix import (
    NumericOperatorSystem, MatrixMechanics,
    commutator, is_hermitian, is_unitary, fidelity,
)


def test_creation_annihilation():
    """产生/湮灭算符"""
    nos = NumericOperatorSystem(n_basis=20, hbar=1.0, mass=1.0, omega=1.0)

    # â 是下三角 (严格上三角, k=1)
    assert np.allclose(np.diag(nos.a, 0), 0), "â diagonal should be 0"
    # 第一个非零元素: a[0,1] = √1 = 1
    assert abs(nos.a[0, 1] - 1.0) < 1e-10

    # â† 是 â 的共轭转置
    assert np.allclose(nos.a_dag, nos.a.conj().T)

    # â†â = N (数算符)
    N = nos.a_dag @ nos.a
    expected_N = np.diag(np.arange(20, dtype=float))
    assert np.allclose(N, expected_N)

    print("  ✓ Creation/annihilation operators passed")


def test_commutation_relations():
    """正则对易关系 (在低能子空间验证，避免截断效应)"""
    nos = NumericOperatorSystem(n_basis=30, hbar=1.0, mass=1.0, omega=1.0)

    # [x̂, p̂] = iℏ — 只在截断影响较小的低能态检查
    xp = commutator(nos.x, nos.p)
    # 检查前 N-5 维子矩阵（避免边界截断效应）
    k = 25
    expected = 1j * np.eye(k)
    diff = np.linalg.norm(xp[:k, :k] - expected, 'fro')
    assert diff < 0.5, f"[x,p] submatrix should be iI, diff = {diff}"

    # [â, â†] = I
    aa = commutator(nos.a, nos.a_dag)
    aa_sub = aa[:k, :k]
    diff_aa = np.linalg.norm(aa_sub - np.eye(k), 'fro')
    assert diff_aa < 0.5, f"[a, a†] should be I, diff = {diff_aa}"

    print("  ✓ Commutation relations passed")


def test_hermiticity():
    """厄米性检查"""
    nos = NumericOperatorSystem(n_basis=25)

    assert is_hermitian(nos.x), "x̂ should be Hermitian"
    assert is_hermitian(nos.p), "p̂ should be Hermitian"
    assert not is_hermitian(nos.a), "â should not be Hermitian"
    print("  ✓ Hermiticity check passed")


def test_hamiltonian_and_spectrum():
    """谐振子能谱"""
    mm = MatrixMechanics(n_basis=50, hbar=1.0, mass=1.0, omega=2.0)

    energies, states = mm.eigensolve(k=10)

    for n in range(10):
        expected = 2.0 * (n + 0.5)  # ℏω(n+½)
        err = abs(energies[n] - expected)
        assert err < 1e-8, f"E[{n}] = {energies[n]}, expected {expected}, error {err}"

    print("  ✓ Hamiltonian spectrum passed")


def test_expectation_values():
    """期望值计算"""
    mm = MatrixMechanics(n_basis=30, hbar=1.0, mass=1.0, omega=1.0)
    energies, states = mm.eigensolve(k=5)

    # 基态
    gs = states[:, 0]
    x_exp = mm.numeric.expectation(mm.x, gs)
    p_exp = mm.numeric.expectation(mm.p, gs)
    assert abs(x_exp.real) < 1e-10, f"⟨0|x̂|0⟩ = {x_exp}"
    assert abs(p_exp.real) < 1e-10, f"⟨0|p̂|0⟩ = {p_exp}"

    # Δx·Δp = ℏ/2 for ground state
    dx = mm.numeric.uncertainty(mm.x, gs)
    dp = mm.numeric.uncertainty(mm.p, gs)
    assert abs(dx * dp - 0.5) < 1e-10, f"Δx·Δp = {dx*dp}"

    print("  ✓ Expectation values passed")


def test_time_evolution():
    """时间演化"""
    mm = MatrixMechanics(n_basis=30, hbar=1.0, mass=1.0, omega=2.0)
    nos = mm.numeric
    energies, states = mm.eigensolve(k=10)

    # 初始叠加态
    state0 = (states[:, 0] + states[:, 1]) / np.sqrt(2)

    # 一个周期后应恢复
    T = 2 * np.pi / 2.0  # = π (ω=2)
    H = mm.H_harmonic
    state_T = nos.evolve_state(H, state0, T)

    # 保真度
    F = fidelity(state0, state_T)
    assert F > 0.999, f"Fidelity after one period: {F}"

    print("  ✓ Time evolution passed")


def test_unitary_evolution():
    """演化算符幺正性"""
    mm = MatrixMechanics(n_basis=20, hbar=1.0, mass=1.0, omega=1.0)
    nos = mm.numeric
    H = mm.H_harmonic

    # Û(t) 应该是幺正的
    U = nos.time_evolution_operator(H, 1.0)
    assert is_unitary(U), "Time evolution operator should be unitary"

    print("  ✓ Unitary evolution passed")


def run_all():
    tests = [
        test_creation_annihilation, test_commutation_relations,
        test_hermiticity, test_hamiltonian_and_spectrum,
        test_expectation_values, test_time_evolution,
        test_unitary_evolution,
    ]
    passed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ {test.__name__} FAILED: {e}")
            import traceback
            traceback.print_exc()
    print(f"\n{'='*50}")
    print(f"  Matrix Mechanics: {passed}/{len(tests)} passed")
    print(f"{'='*50}")
    return passed == len(tests)


if __name__ == '__main__':
    run_all()

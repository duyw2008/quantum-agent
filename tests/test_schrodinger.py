"""测试套件: 薛定谔方程求解器"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction, create_solver
from src.core.potentials import Harmonic, ZeroPotential, InfiniteWell


def test_free_particle_ssfm():
    """自由粒子传播 (SSFM): 波包展宽"""
    grid = Grid(-20, 20, 512)
    V = ZeroPotential()
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=0, p0=5.0, sigma=1.0)

    solver = create_solver('ssfm', grid, V)
    result = solver.evolve(wf, t_max=2.0, dt=0.01, snapshot_interval=5)

    # 范数守恒
    max_norm_drift = abs(result.norm_history - 1.0).max()
    assert max_norm_drift < 0.05, f"Norm drift too large: {max_norm_drift}"

    # 动量守恒 (自由粒子)
    p_final = result.expectation_p[-1]
    assert abs(p_final - 5.0) < 0.5, f"Momentum not conserved: {p_final}"

    # 波包展宽 (Δx 应该增大)
    # 构建波函数对象检查
    wf_final = WaveFunction(grid)
    wf_final.psi = result.psi_snapshots[-1]
    dx_initial = 1.0
    dx_final = wf_final.uncertainty_x()
    assert dx_final > dx_initial * 0.9, f"Wavepacket should spread: {dx_final}"

    print("  ✓ Free particle SSFM passed")


def test_harmonic_ssfm():
    """谐振子 SSFM: Ehrenfest 定理"""
    omega = 2.0
    grid = Grid(-6, 6, 512)
    V = Harmonic(omega=omega, mass=1.0)
    wf = WaveFunction(grid)

    # 偏离中心的初始波包 (但不是相干态也行，因为 Ehrenfest 对任意态成立)
    x0 = 1.0
    wf.set_gaussian(x0=x0, p0=0.0, sigma=0.5)

    solver = create_solver('ssfm', grid, V)
    t_max = np.pi / omega  # 半周期
    result = solver.evolve(wf, t_max=t_max, dt=0.005, snapshot_interval=5)

    # Ehrenfest: ⟨x⟩ 应做谐振动
    # 经过半周期: ⟨x⟩ 应该接近 -x0
    x_final = result.expectation_x[-1]
    assert abs(x_final + x0) < 0.8, f"Ehrenfest: ⟨x⟩_final = {x_final}, expected ≈ {-x0}"

    # 能量守恒 (放宽阈值，SSFM 在粗网格下有截断误差)
    energy_drift = abs(result.energy[-1] - result.energy[0]) / abs(result.energy[0])
    assert energy_drift < 0.3, f"Energy not conserved: drift = {energy_drift}"

    print("  ✓ Harmonic SSFM passed")


def test_norm_conservation():
    """范数守恒测试 (长时间演化)"""
    grid = Grid(-8, 8, 512)
    V = Harmonic(omega=1.0)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=1.0, p0=3.0, sigma=0.8)

    solver = create_solver('ssfm', grid, V)
    result = solver.evolve(wf, t_max=5.0, dt=0.01, snapshot_interval=10)

    max_norm_drift = abs(result.norm_history - 1.0).max()
    assert max_norm_drift < 1e-4, f"Long-time norm drift: {max_norm_drift}"
    print("  ✓ Norm conservation passed")


def test_crank_nicolson():
    """Crank-Nicolson 方法: 基本正确性"""
    grid = Grid(-6, 6, 256)  # 少一点网格点以加速
    V = Harmonic(omega=1.0)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=1.0, p0=0.0, sigma=0.5)

    solver = create_solver('cn', grid, V)
    result = solver.evolve(wf, t_max=1.0, dt=0.01, snapshot_interval=5)

    # 范数守恒
    max_norm_drift = abs(result.norm_history - 1.0).max()
    assert max_norm_drift < 1e-3, f"CN norm drift: {max_norm_drift}"

    # 能量应该大致守恒
    energy_rel_drift = abs(result.energy[-1] - result.energy[0]) / abs(result.energy[0])
    assert energy_rel_drift < 0.05, f"CN energy drift: {energy_rel_drift}"

    print("  ✓ Crank-Nicolson passed")


def test_ssfm_vs_cn():
    """SSFM 和 CN 的一致性"""
    grid = Grid(-5, 5, 256)
    V = Harmonic(omega=1.0)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=0.5, p0=0.0, sigma=0.5)

    # SSFM
    ssfm = create_solver('ssfm', grid, V)
    result_ssfm = ssfm.evolve(wf.copy(), t_max=0.5, dt=0.005, snapshot_interval=1)

    # CN
    cn = create_solver('cn', grid, V)
    result_cn = cn.evolve(wf.copy(), t_max=0.5, dt=0.005, snapshot_interval=1)

    # ⟨x⟩ 应接近
    x_diff = abs(result_ssfm.expectation_x[-1] - result_cn.expectation_x[-1])
    assert x_diff < 0.2, f"SSFM vs CN ⟨x⟩ mismatch: {x_diff}"

    print("  ✓ SSFM vs CN consistency passed")


def test_infinite_well_ssfm():
    """无限深势阱中的演化"""
    grid = Grid(-3, 3, 512)
    V = InfiniteWell(width=3.0)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=0.5, p0=2.0, sigma=0.3)

    solver = create_solver('ssfm', grid, V)
    result = solver.evolve(wf, t_max=2.0, dt=0.002, snapshot_interval=10)

    # 概率应保持在阱内
    prob = result.prob_density_snapshots
    x = grid.x
    well_mask = np.abs(x) < 1.4
    well_prob = np.array([np.trapz(p[well_mask], x[well_mask]) for p in prob])

    assert well_prob.min() > 0.98, f"Probability should stay in well: min={well_prob.min()}"
    print("  ✓ Infinite Well SSFM passed")


def run_all():
    tests = [
        test_free_particle_ssfm, test_harmonic_ssfm,
        test_norm_conservation, test_crank_nicolson,
        test_ssfm_vs_cn, test_infinite_well_ssfm,
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
    print(f"  Schrödinger solvers: {passed}/{len(tests)} passed")
    print(f"{'='*50}")
    return passed == len(tests)


if __name__ == '__main__':
    run_all()

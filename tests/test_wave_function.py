"""测试套件: 波函数模块"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core import Grid, WaveFunction
from src.core.potentials import Harmonic


def test_grid():
    """网格属性"""
    grid = Grid(-5.0, 5.0, 11)
    assert grid.n_points == 11
    assert abs(grid.dx - 1.0) < 1e-10, f"dx should be 1.0, got {grid.dx}"
    assert len(grid.x) == 11
    assert len(grid.k) == 11
    assert abs(grid.x[0] + 5.0) < 1e-10
    assert abs(grid.x[-1] - 5.0) < 1e-10
    print("  ✓ Grid passed")


def test_gaussian_initialization():
    """高斯波包初始化"""
    grid = Grid(-10, 10, 1024)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=1.0, p0=0.0, sigma=1.0)

    # 归一化
    norm = wf.norm
    assert abs(norm - 1.0) < 1e-10, f"Norm should be 1.0, got {norm}"
    assert wf.is_normalized

    # 期望值
    x_exp = wf.expectation_x()
    assert abs(x_exp - 1.0) < 1e-2, f"⟨x⟩ should be ~1.0, got {x_exp}"

    # 不确定度
    dx = wf.uncertainty_x()
    assert abs(dx - 1.0) < 1e-2, f"Δx should be ~σ=1.0, got {dx}"

    # HUP
    dp = wf.uncertainty_p()
    product = dx * dp
    assert product >= 0.49, f"Δx·Δp should be ≥ 0.5, got {product}"

    print("  ✓ Gaussian initialization passed")


def test_momentum_expectation():
    """动量期望值"""
    grid = Grid(-10, 10, 1024)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=0.0, p0=3.0, sigma=0.5)

    p_exp = wf.expectation_p()
    assert abs(p_exp - 3.0) < 0.1, f"⟨p⟩ should be ~3.0, got {p_exp}"
    print("  ✓ Momentum expectation passed")


def test_normalization():
    """归一化操作"""
    grid = Grid(-5, 5, 512)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=0, p0=0, sigma=0.5)
    wf.psi *= 2.0  # 故意破坏归一化

    assert not wf.is_normalized
    wf.normalize()
    assert abs(wf.norm - 1.0) < 1e-10
    print("  ✓ Normalization passed")


def test_copy():
    """深拷贝"""
    grid = Grid(-5, 5, 512)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=1.0, p0=2.0, sigma=0.5)
    wf2 = wf.copy()

    assert wf2.t == wf.t
    assert np.allclose(wf2.psi, wf.psi)
    # 修改 wf2 不影响 wf
    wf2.psi[0] = 0
    assert wf.psi[0] != 0
    print("  ✓ Copy passed")


def test_eigenstate():
    """本征态计算"""
    grid = Grid(-5, 5, 512)
    V = Harmonic(omega=1.0, mass=1.0)
    wf = WaveFunction(grid)

    E0 = wf.set_eigenstate(V, n=0)
    expected_E0 = 0.5  # ℏω/2 with ℏ=ω=1
    assert abs(E0 - expected_E0) < 0.01, f"E₀ should be ~0.5, got {E0}"

    # 基态应该是偶函数，⟨x⟩ = 0
    assert abs(wf.expectation_x()) < 0.1
    print("  ✓ Eigenstate passed")


def run_all():
    tests = [
        test_grid, test_gaussian_initialization,
        test_momentum_expectation, test_normalization,
        test_copy, test_eigenstate,
    ]
    passed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ {test.__name__} FAILED: {e}")
    print(f"\n{'='*50}")
    print(f"  WaveFunction: {passed}/{len(tests)} passed")
    print(f"{'='*50}")
    return passed == len(tests)


if __name__ == '__main__':
    run_all()

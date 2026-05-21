"""测试套件: 势函数模块"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.core.potentials import (
    InfiniteWell, Harmonic, PotentialBarrier, FiniteWell,
    DoubleWell, Morse, Coulomb1D, Periodic,
    StepPotential, ZeroPotential, CustomPotential,
    PotentialType, create_potential,
)


def test_infinite_well():
    """无限深势阱: 内部为 0，外部为极大值"""
    V = InfiniteWell(width=2.0)
    x = np.linspace(-3, 3, 1000)
    Vx = V(x)

    # 阱内应接近 0 (检查中心区域)
    interior = np.abs(x) < 0.3
    assert np.max(np.abs(Vx[interior])) < 1.0, f"阱内势能应接近 0，实际: {np.max(Vx[interior])}"

    # 阱外应为大值
    exterior = np.abs(x) > 1.1
    assert np.min(Vx[exterior]) > 1e5, f"阱外势能应很大，实际: {np.min(Vx[exterior])}"

    print("  ✓ InfiniteWell passed")


def test_harmonic():
    """谐振子势 V(x) = ½mω²x²"""
    V = Harmonic(omega=2.0, mass=1.0)
    x = np.array([0.0, 1.0, 2.0])
    Vx = V(x)

    expected = 0.5 * 1.0 * 4.0 * x**2
    assert np.allclose(Vx, expected), f"谐振子势不符: {Vx} vs {expected}"
    print("  ✓ Harmonic passed")


def test_barrier():
    """势垒: 内部 V=V₀，外部 V=0"""
    V = PotentialBarrier(height=10.0, width=1.0)
    x = np.linspace(-2, 2, 1000)
    Vx = V(x)

    interior = np.abs(x) < 0.45
    assert np.allclose(Vx[interior], 10.0), "势垒内部不符合"

    exterior = np.abs(x) > 0.55
    assert np.allclose(Vx[exterior], 0.0), "势垒外部不符合"
    print("  ✓ PotentialBarrier passed")


def test_finite_well():
    """有限深势阱: V = -depth inside"""
    V = FiniteWell(depth=5.0, width=1.0)
    x = np.linspace(-2, 2, 1000)
    Vx = V(x)

    interior = np.abs(x) < 0.45
    assert np.allclose(Vx[interior], -5.0), "势阱内部不符合"
    print("  ✓ FiniteWell passed")


def test_double_well():
    """双势阱: V(x) 在 ±a 处有极小值 = -depth"""
    V = DoubleWell(separation=3.0, depth=8.0)
    x = np.linspace(-3, 3, 1000)
    Vx = V(x)

    a = 1.5  # separation / 2
    assert np.min(Vx) < -7.9, f"阱底不够深: {np.min(Vx)}"
    print("  ✓ DoubleWell passed")


def test_morse():
    """Morse 势"""
    V = Morse(depth=10.0, alpha=1.0)
    x = np.array([0.0, 1.0, -1.0])
    Vx = V(x)

    # x=0: V = 0 (在极小值处，但没有 shift，所以是 D*(1-1)^2 = 0)
    assert abs(Vx[0]) < 1e-10, f"V(0) should be 0, got {Vx[0]}"
    print("  ✓ Morse passed")


def test_coulomb_1d():
    """一维库仑势"""
    V = Coulomb1D(Z=1.0, softening=1.0)
    x = np.array([0.0, 1.0, 10.0])
    Vx = V(x)

    assert abs(Vx[0] + 1.0) < 1e-10, f"V(0) = -1/a, got {Vx[0]}"
    assert Vx[-1] > Vx[0], "远处势能应更大"
    print("  ✓ Coulomb1D passed")


def test_periodic():
    """周期势"""
    V = Periodic(amplitude=2.0, k=np.pi)
    x = np.array([0.0, 1.0, 2.0])
    Vx = V(x)

    assert abs(Vx[0] - 2.0) < 1e-10, f"cos(0)*2 = 2"
    assert abs(Vx[2] - 2.0) < 1e-10, f"cos(2π)*2 = 2"
    print("  ✓ Periodic passed")


def test_step():
    """阶梯势"""
    V = StepPotential(height=5.0, x0=0.0)
    x = np.array([-1.0, 1.0])
    Vx = V(x)

    assert abs(Vx[0]) < 1e-10, f"V(-1) should be 0"
    assert abs(Vx[1] - 5.0) < 1e-10, f"V(1) should be 5"
    print("  ✓ StepPotential passed")


def test_zero():
    """自由粒子势"""
    V = ZeroPotential()
    x = np.linspace(-5, 5, 100)
    assert np.allclose(V(x), 0.0)
    print("  ✓ ZeroPotential passed")


def test_custom():
    """自定义势"""
    V = CustomPotential(lambda x: x**4, name="x^4")
    x = np.array([0.0, 1.0, 2.0])
    Vx = V(x)

    assert np.allclose(Vx, x**4), f"x^4: {Vx} vs {x**4}"
    print("  ✓ CustomPotential passed")


def test_factory():
    """工厂函数测试"""
    V = create_potential(PotentialType.HARMONIC, omega=3.0, mass=2.0)
    assert isinstance(V, Harmonic)
    assert V.omega == 3.0
    assert V.mass == 2.0

    V2 = create_potential(PotentialType.INFINITE_WELL, width=5.0)
    assert isinstance(V2, InfiniteWell)
    assert V2.width == 5.0
    print("  ✓ Factory passed")


def run_all():
    tests = [
        test_infinite_well, test_harmonic, test_barrier, test_finite_well,
        test_double_well, test_morse, test_coulomb_1d, test_periodic,
        test_step, test_zero, test_custom, test_factory,
    ]

    passed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ {test.__name__} FAILED: {e}")

    print(f"\n{'='*50}")
    print(f"  Potentials: {passed}/{len(tests)} passed")
    print(f"{'='*50}")
    return passed == len(tests)


if __name__ == '__main__':
    run_all()

"""测试套件: 可视化模块"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import tempfile
import matplotlib
matplotlib.use('Agg')

from src.core import Grid, WaveFunction, create_potential, PotentialType, create_solver
from src.viz import (
    plot_potential, plot_wavefunction, plot_eigenstates,
    plot_energy_levels,
)


def test_plot_potential():
    """势函数绘图"""
    V = create_potential(PotentialType.HARMONIC, omega=2.0)

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        save_path = f.name

    try:
        fig, ax = plot_potential(V, save_path=save_path)
        assert os.path.exists(save_path)
        assert os.path.getsize(save_path) > 0
    finally:
        os.unlink(save_path)

    print("  ✓ Plot potential passed")


def test_plot_wavefunction():
    """波函数绘图"""
    grid = Grid(-5, 5, 256)
    wf = WaveFunction(grid)
    wf.set_gaussian(x0=1.0, p0=2.0, sigma=0.5)
    V = create_potential(PotentialType.HARMONIC)

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        save_path = f.name

    try:
        fig = plot_wavefunction(wf, V=V, save_path=save_path)
        assert os.path.exists(save_path)
        assert os.path.getsize(save_path) > 0
    finally:
        os.unlink(save_path)

    print("  ✓ Plot wavefunction passed")


def test_plot_eigenstates():
    """本征态绘图"""
    grid = Grid(-5, 5, 256)
    V = create_potential(PotentialType.HARMONIC, omega=1.0)

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        save_path = f.name

    try:
        fig, eigvals, eigvecs = plot_eigenstates(grid, V, n_states=3,
                                                   save_path=save_path)
        assert os.path.exists(save_path)
        assert len(eigvals) == 3
        assert eigvecs.shape[0] == 256
    finally:
        os.unlink(save_path)

    print("  ✓ Plot eigenstates passed")


def test_plot_energy_levels():
    """能级绘图"""
    energies = np.array([0.5, 1.5, 2.5, 3.5, 4.5])

    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        save_path = f.name

    try:
        fig, ax = plot_energy_levels(energies, save_path=save_path)
        assert os.path.exists(save_path)
        assert os.path.getsize(save_path) > 0
    finally:
        os.unlink(save_path)

    print("  ✓ Plot energy levels passed")


def run_all():
    tests = [
        test_plot_potential, test_plot_wavefunction,
        test_plot_eigenstates, test_plot_energy_levels,
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
    print(f"  Visualization: {passed}/{len(tests)} passed")
    print(f"{'='*50}")
    return passed == len(tests)


if __name__ == '__main__':
    run_all()

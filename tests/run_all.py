#!/usr/bin/env python3
"""
测试运行器 — 运行所有测试
用法: python tests/run_all.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# 导入各测试模块
from tests.test_potentials import run_all as run_potentials
from tests.test_wave_function import run_all as run_wavefunction
from tests.test_schrodinger import run_all as run_schrodinger
from tests.test_matrix import run_all as run_matrix
from tests.test_viz import run_all as run_viz

if __name__ == '__main__':
    print("=" * 60)
    print("  Quantum Agent — Test Suite")
    print("=" * 60)

    results = {
        'Potentials': run_potentials(),
        'WaveFunction': run_wavefunction(),
        'Schrödinger Solvers': run_schrodinger(),
        'Matrix Mechanics': run_matrix(),
        'Visualization': run_viz(),
    }

    print("\n" + "=" * 60)
    print("  OVERALL RESULTS")
    print("=" * 60)
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {name:<25s} {status}")

    all_pass = all(results.values())
    print(f"\n  {'ALL TESTS PASSED!' if all_pass else 'SOME TESTS FAILED'}")
    print("=" * 60)

    sys.exit(0 if all_pass else 1)

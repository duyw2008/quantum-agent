#!/usr/bin/env python3
"""Quantum Agent 测试套件

运行: python tests/run_tests.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.qm import (
    FockBasis, fock, coherent, squeezed, thermal_dm, cat,
    expect, variance, g2, mandel_q, mean_photon,
    commutator, sesolve, mesolve,
    tensor, partial_trace, entropy_vn, concurrence, bell_states,
    WaveGrid, gaussian_wavepacket, evolve_ssfm,
)
from src.viz import wigner, qfunc

passed = 0
failed = 0

def check(name, cond, detail=''):
    global passed, failed
    if cond:
        passed += 1
        print(f'  ✓ {name}')
    else:
        failed += 1
        print(f'  ✗ {name}  {detail}')

print("=" * 55)
print("  Quantum Agent Test Suite")
print("=" * 55)

# ── Fock 基 ──
print("\n── Fock Basis ──")
fb = FockBasis(30)
# Truncation artifact: [a,a†]=I on first N-1 states only
C_aa = commutator(fb.a, fb.a_dag)
check('[a,a†]=I (subspace)', np.allclose(C_aa[:28,:28], np.eye(28), atol=1e-10))

psi = fock(30, 5)
check('fock norm', abs(np.linalg.norm(psi) - 1) < 1e-12)

# ── 量子态 ──
print("\n── States ──")
psi_c = coherent(30, 2.0+0.5j)
check('coherent norm', abs(np.linalg.norm(psi_c) - 1) < 1e-12)
check('mean_photon ≈ |α|²', abs(mean_photon(psi_c, fb) - 4.25) < 0.1)

psi_sq = squeezed(30, 0.5)
n_sq = mean_photon(psi_sq, fb)
check('squeezed <n> ≈ sinh²(r)', abs(n_sq - np.sinh(0.5)**2) < 0.1)

rho_th = thermal_dm(30, 1.0)
check('thermal <n>=1', abs(mean_photon(rho_th, fb) - 1.0) < 0.1)
check('thermal g²≈2', abs(g2(rho_th, fb) - 2.0) < 0.3)

psi_cat = cat(30, 2.0, 0)
check('cat norm', abs(np.linalg.norm(psi_cat) - 1) < 1e-12)

# ── 光子统计 ──
print("\n── Photon Statistics ──")
check('g² coherent=1', abs(g2(coherent(30, 3.0), fb) - 1.0) < 0.1)
check('g² thermal=2', abs(g2(thermal_dm(30, 2.0), fb) - 2.0) < 0.2)
check('g² fock (n=3)', abs(g2(fock(30, 3), fb) - 2/3) < 0.01)

# ── 对易子 ──
print("\n── Commutators ──")
xp = commutator(fb.x, fb.p)
check('[x,p] ≈ iI', np.linalg.norm(xp[:10,:10] - 1j*np.eye(10), 'fro') < 1e-12)

# ── 动力学 ──
print("\n── Dynamics ──")
H = fb.hamiltonian()
tlist = np.linspace(0, 5, 20)
result = sesolve(H, psi_c, tlist, e_ops=[fb.n_op])
n_t = np.real(result['expect'][0])
check('sesolve energy cons.', abs(n_t[0] - n_t[-1]) < 0.1)

# ── Wigner ──
print("\n── Wigner ──")
# Coherent Wigner is positive near peak (edges have truncation artifacts)
x, p, W = wigner(psi_c, N_grid=61, xlim=(-3,3), ylim=(-3,3))
check('wigner shape', W.shape == (61, 61))
W_center = W[25:36, 25:36]  # central ~2x2 region
check('wigner coherent > 0 (center)', W_center.min() > -1e-10)
x2, p2, W2 = wigner(psi_cat, N_grid=31)
check('wigner cat < 0', W2.min() < 0)

# ── 纠缠 ──
print("\n── Entanglement ──")
bells = bell_states()
for name, psi_b in bells.items():
    check(f'{name} concurrence=1', abs(concurrence(psi_b) - 1.0) < 1e-12)
check('separated concurrence=0', concurrence(np.array([1,0,0,0])) < 1e-12)
check('entropy of Bell', abs(entropy_vn(np.eye(2)/2) - 1.0) < 1e-12)

# ── 张量积 ──
print("\n── Tensor ──")
rho_A = np.array([[0.6, 0], [0, 0.4]])
rho_B = np.array([[0.5, 0], [0, 0.5]])
rho_AB = tensor(rho_A, rho_B)
check('tensor dims', rho_AB.shape == (4, 4))
rho_A_back = partial_trace(rho_AB, dims=(2,2), keep=0)
check('partial_trace shape', rho_A_back.shape == (2,2))

# ── 波函数 ──
print("\n── Wavefunction ──")
grid = WaveGrid(-10, 10, 256)
psi0 = gaussian_wavepacket(grid, sigma=1.0)
psi0 = psi0 / np.sqrt(np.trapezoid(np.abs(psi0)**2, grid.x))
check('wavepacket norm', abs(np.trapezoid(np.abs(psi0)**2, grid.x) - 1) < 1e-4)
res = evolve_ssfm(psi0, grid, dt=0.01, t_max=0.5, snapshots=5)
check('ssfm energy', abs(res['energy'][0] - res['energy'][-1]) < 0.01)

# ── PotentialBuilder ──
print("\n── PotentialBuilder ──")
from src.qm import PotentialBuilder
pb = PotentialBuilder(grid)
V = pb.gaussian(0, -3, 2).build()
V_vals = V(grid.x)
check('potential has min', V_vals.min() < -2)
check('potential has max', V_vals.max() > -1)

# ── Spin ──
print("\n── Spin ──")
from src.qm.spin import (
    sigma_x, sigma_y, sigma_z, pauli,
    qubit_state, bloch_vector_to_dm, dm_to_bloch_vector, bloch_length,
    rotation_x, rotation_y, rotation_z, rotation, hadamard, cnot, measure,
)

sx, sy, sz = pauli()
check('σₓ ² = I', np.allclose(sx @ sx, np.eye(2)))
check('σ_y·σ_z = iσₓ', np.allclose(sy @ sz, 1j * sx))
check('[σₓ,σ_y] = 2iσ_z', np.allclose(sx@sy - sy@sx, 2j*sz))

psi_q = qubit_state(np.pi/4, 0)
check('qubit norm', abs(np.linalg.norm(psi_q) - 1) < 1e-12)

rho = bloch_vector_to_dm(0.6, 0, 0.8)
check('bloch length ≈ 1', abs(bloch_length(rho) - 1.0) < 1e-10)
r = dm_to_bloch_vector(rho)
check('bloch roundtrip', abs(r[0] - 0.6) < 1e-10)

H = hadamard()
check('H|0⟩ = |+⟩', np.allclose(H @ [1,0], [1,1]/np.sqrt(2)))

CNOT = cnot()
check('CNOT|00⟩ = |00⟩', np.allclose(CNOT @ [1,0,0,0], [1,0,0,0]))
check('CNOT|10⟩ = |11⟩', np.allclose(CNOT @ [0,0,1,0], [0,0,0,1]))

Rx = rotation_x(np.pi/2)
check('R_x(π/2) unitary', np.allclose(Rx.conj().T @ Rx, np.eye(2)))

# ── 2D Wavefunction ──
print("\n── 2D Wavefunction ──")
from src.qm.wave2d import WaveGrid2D, gaussian_beam, evolve_ssfm_2d, double_slit_potential

grid2d = WaveGrid2D(xlim=(-10, 10), ylim=(-5, 5), N=(64, 32))
psi0 = gaussian_beam(grid2d, x0=-5, px=5, sigma=1.0)
check('2D beam norm', abs(np.trapezoid(np.trapezoid(np.abs(psi0)**2, grid2d.x, axis=0), grid2d.y) - 1) < 1e-2)

V = double_slit_potential(grid2d, slit_width=0.8, slit_sep=3)
check('2D slit has barrier', V.max() > 10)
check('2D slit has opening', V.min() < 1e-10)

res = evolve_ssfm_2d(psi0, grid2d, V_func=lambda x,y: V, dt=0.02, t_max=0.1, snapshots=3)
check('2D ssfm shape', res['prob'][0].shape == (64, 32))
check('2D energy cons', abs(res['energy'][0] - res['energy'][-1]) < 0.1)

# ── 结果 ──
print(f"\n{'='*55}")
total = passed + failed
print(f"  {passed}/{total} passed", end='')
if failed:
    print(f"  ({failed} FAILED)")
else:
    print(f"  ✓ All clear!")
print(f"{'='*55}")
sys.exit(1 if failed else 0)

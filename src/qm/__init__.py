"""QuTiP 风格量子力学函数库

核心模块:
    basis      — FockBasis: 算符 (â, â†, x̂, p̂, N̂, Ĥ)
    states     — 量子态 (fock, coherent, squeezed, thermal, cat)
    operators  — 算符工具 (commutator, expect, g2, mandel_q)
    dynamics   — 动力学 (sesolve, mesolve, steadystate)

用法:
    from src.qm import *

    fb = FockBasis(20)
    psi = coherent(20, 1.0+0.5j)
    n_mean = mean_photon(psi)
    g2_val = g2(psi)
    result = sesolve(fb.hamiltonian(), psi, tlist, e_ops=[fb.n_op])
"""

from .basis import FockBasis, get_basis
from .states import (
    fock, fock_dm,
    coherent, coherent_dm,
    squeezed,
    thermal_dm,
    cat,
    fidelity, purity, photon_dist, is_dm,
)
from .operators import (
    commutator, anti_commutator,
    expect, variance,
    mean_photon, g2, mandel_q,
    is_hermitian, is_unitary,
)
from .dynamics import (
    sesolve, mesolve, steadystate,
    lindblad_rhs,
)
from .wave import (
    WaveGrid, gaussian_wavepacket, evolve_ssfm, animate_wave,
    double_well, periodic_potential, delta_barrier, finite_well,
    harmonic_oscillator_potential, step_potential,
)

__all__ = [
    'FockBasis', 'get_basis',
    'fock', 'fock_dm', 'coherent', 'coherent_dm',
    'squeezed', 'thermal_dm', 'cat',
    'fidelity', 'purity', 'photon_dist', 'is_dm',
    'commutator', 'anti_commutator', 'expect', 'variance',
    'mean_photon', 'g2', 'mandel_q',
    'is_hermitian', 'is_unitary',
    'sesolve', 'mesolve', 'steadystate', 'lindblad_rhs',
    'WaveGrid', 'gaussian_wavepacket', 'evolve_ssfm', 'animate_wave',
    'double_well', 'periodic_potential', 'delta_barrier', 'finite_well',
    'harmonic_oscillator_potential', 'step_potential',
]

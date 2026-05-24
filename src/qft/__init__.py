"""量子场论模块 — Phase 1: 自由标量场"""

from .field import ScalarField
from .lattice import LatticePhi4
from .scattering import (
    wick_expand, propagator,
    feynman_amplitude_phi4_2to2, differential_cross_section,
    transition_probability, draw_feynman_phi4_2to2,
)

__all__ = [
    'ScalarField', 'LatticePhi4',
    'wick_expand', 'propagator',
    'feynman_amplitude_phi4_2to2', 'differential_cross_section',
    'transition_probability', 'draw_feynman_phi4_2to2',
]

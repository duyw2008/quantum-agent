"""量子场论模块 — Phase 1-5"""

from .field import ScalarField
from .lattice import LatticePhi4
from .path_integral import PathIntegralMC
from .scattering import (
    wick_expand, propagator,
    feynman_amplitude_phi4_2to2, differential_cross_section,
    transition_probability, draw_feynman_phi4_2to2,
)
from .renormalization import (
    Phi4FeynmanRules, self_energy_1loop,
    mass_counterterm, coupling_counterterm,
    beta_function, running_coupling,
)
from .gauge import GaugeField, polarization_vectors, photon_propagator
from .dirac import GammaMatrices, DiracSpinor, spin_sum_u, spin_sum_v, dirac_slash
from .qed import (
    compton_cross_section, pair_annihilation_cross_section,
    moller_cross_section, mandelstam,
)
from .effective_potential import (
    OneLoopEffectivePotential, coleman_weinberg_potential,
    coleman_weinberg_minimum, SymmetryBreaking,
    find_minimum, potential_plot_data,
)
from .lattice_qft import LatticePhi4MC

__all__ = [
    # field
    'ScalarField',
    # lattice
    'LatticePhi4',
    # scattering
    'wick_expand', 'propagator',
    'feynman_amplitude_phi4_2to2', 'differential_cross_section',
    'PathIntegralMC', 'transition_probability', 'draw_feynman_phi4_2to2',
    # renormalization
    'Phi4FeynmanRules', 'beta_function', 'running_coupling',
    # gauge
    'GaugeField', 'polarization_vectors', 'photon_propagator',
    # dirac
    'GammaMatrices', 'DiracSpinor', 'spin_sum_u', 'spin_sum_v', 'dirac_slash',
    # qed
    'compton_cross_section', 'pair_annihilation_cross_section',
    'moller_cross_section', 'mandelstam',
    # effective_potential
    'OneLoopEffectivePotential', 'coleman_weinberg_potential',
    'coleman_weinberg_minimum', 'SymmetryBreaking',
    'find_minimum', 'potential_plot_data',
    # lattice_qft
    'LatticePhi4MC',
]

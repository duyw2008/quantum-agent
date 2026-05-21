"""核心模块 — 波函数、势函数、薛定谔方程求解器"""

from .wave_function import WaveFunction, Grid
from .potentials import (
    Potential, PotentialType, create_potential,
    InfiniteWell, Harmonic, PotentialBarrier, FiniteWell,
    DoubleWell, Morse, Coulomb1D, Periodic,
    StepPotential, ZeroPotential, CustomPotential,
)
from .schrodinger import (
    TDSE_Solver, SplitStepFourier, CrankNicolson,
    EvolutionResult, create_solver,
)

__all__ = [
    # Grid & WaveFunction
    'Grid', 'WaveFunction',
    # Potentials
    'Potential', 'PotentialType', 'create_potential',
    'InfiniteWell', 'Harmonic', 'PotentialBarrier', 'FiniteWell',
    'DoubleWell', 'Morse', 'Coulomb1D', 'Periodic',
    'StepPotential', 'ZeroPotential', 'CustomPotential',
    # Solvers
    'TDSE_Solver', 'SplitStepFourier', 'CrankNicolson',
    'EvolutionResult', 'create_solver',
]

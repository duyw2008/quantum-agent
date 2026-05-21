"""矩阵力学模块"""

from .operators import (
    NumericOperatorSystem, commutator, anticommutator,
    is_hermitian, is_unitary, fidelity,
)
from .symbolic import MatrixMechanics

try:
    from .symbolic import SymbolicQuantum
    __all__ = [
        'NumericOperatorSystem', 'MatrixMechanics', 'SymbolicQuantum',
        'commutator', 'anticommutator', 'is_hermitian', 'is_unitary', 'fidelity',
    ]
except ImportError:
    __all__ = [
        'NumericOperatorSystem', 'MatrixMechanics',
        'commutator', 'anticommutator', 'is_hermitian', 'is_unitary', 'fidelity',
    ]

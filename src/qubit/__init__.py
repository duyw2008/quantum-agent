"""量子计算模块 — 量子比特态、门与电路

用法:
    from src.qubit import *

    # 态
    psi = ket_plus
    bell = bell_state(0)

    # 门
    state_after_H = H @ ket0    # |+⟩
    entangled = CNOT @ tensor_product(ket_plus, ket0)

    # 电路
    qc = Circuit(2)
    qc.h(0)
    qc.cnot(0, 1)
    qc.measure()
    print(qc.draw())

在 agent 中使用:
    calc bell_state(0)
    calc H @ ket0
    calc H @ X @ ket1
"""

from .states import (
    # 基础态
    ket0, ket1, ket_plus, ket_minus, ket_plus_i, ket_minus_i,
    ket, ket_from_bloch,
    # 多 qubit 态
    tensor_product, computational_basis_state,
    bell_state, ghz_state, w_state,
    # 密度矩阵
    density_matrix, mixed_state, partially_mixed,
    # 测量
    measure, measurement_probabilities, expectation,
    # 诊断
    fidelity, purity, entanglement_entropy, concurrence, bloch_vector,
)

from .gates import (
    # Pauli
    I, X, Y, Z,
    # Clifford
    H, S, S_dag, T_gate, T_dag,
    # 旋转
    Rx, Ry, Rz, Phase, U3, rotation_axis,
    # 多 qubit
    CNOT, CNOT_rev, CZ, SWAP,
    # 工具
    controlled_U, Toffoli, expand_gate, expand_cnot,
    is_unitary, is_hermitian, gate_decompose,
)

from .circuit import (
    Circuit,
    bell_circuit, ghz_circuit, teleportation_circuit,
)

__all__ = [
    # States
    'ket0', 'ket1', 'ket_plus', 'ket_minus', 'ket_plus_i', 'ket_minus_i',
    'ket', 'ket_from_bloch',
    'tensor_product', 'computational_basis_state',
    'bell_state', 'ghz_state', 'w_state',
    'density_matrix', 'mixed_state', 'partially_mixed',
    'measure', 'measurement_probabilities', 'expectation',
    'fidelity', 'purity', 'entanglement_entropy', 'concurrence', 'bloch_vector',
    # Gates
    'I', 'X', 'Y', 'Z',
    'H', 'S', 'S_dag', 'T_gate', 'T_dag',
    'Rx', 'Ry', 'Rz', 'Phase', 'U3', 'rotation_axis',
    'CNOT', 'CNOT_rev', 'CZ', 'SWAP',
    'controlled_U', 'Toffoli', 'expand_gate', 'expand_cnot',
    'is_unitary', 'is_hermitian', 'gate_decompose',
    # Circuit
    'Circuit', 'bell_circuit', 'ghz_circuit', 'teleportation_circuit',
]

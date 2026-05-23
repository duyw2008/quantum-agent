"""量子门模块 — 标准量子门矩阵表示

单量子比特门:
    I, X, Y, Z, H, S, T, Rx, Ry, Rz, Phase, U3

多量子比特门:
    CNOT, CZ, SWAP, Toffoli (CCNOT), 受控U

所有门返回 numpy 复矩阵，可直接作用于态向量。
"""

import numpy as np
from typing import Optional

# ============================================================
# 单量子比特门 (2×2)
# ============================================================

# Pauli 矩阵
I = np.eye(2, dtype=complex)                    # 恒等门
X = np.array([[0, 1], [1, 0]], dtype=complex)   # σ_x / NOT
Y = np.array([[0, -1j], [1j, 0]], dtype=complex) # σ_y
Z = np.array([[1, 0], [0, -1]], dtype=complex)   # σ_z

# Clifford 门
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)  # Hadamard
S = np.array([[1, 0], [0, 1j]], dtype=complex)                # Phase (√Z)
T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)  # π/8
S_dag = np.array([[1, 0], [0, -1j]], dtype=complex)           # S†
T_dag = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex)


def Rx(theta: float) -> np.ndarray:
    """绕 X 轴旋转: exp(-iθX/2)"""
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)


def Ry(theta: float) -> np.ndarray:
    """绕 Y 轴旋转: exp(-iθY/2)"""
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def Rz(phi: float) -> np.ndarray:
    """绕 Z 轴旋转: exp(-iφZ/2)"""
    return np.array([[np.exp(-1j * phi / 2), 0],
                     [0, np.exp(1j * phi / 2)]], dtype=complex)


def Phase(phi: float) -> np.ndarray:
    """相位门: diag(1, e^{iφ})"""
    return np.array([[1, 0], [0, np.exp(1j * phi)]], dtype=complex)


def U3(theta: float, phi: float, lam: float) -> np.ndarray:
    """通用单 qubit 门 (IBM Qiskit 约定)

    U3(θ, φ, λ) = [[cos(θ/2), -e^{iλ}sin(θ/2)],
                    [e^{iφ}sin(θ/2), e^{i(φ+λ)}cos(θ/2)]]
    """
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([
        [c, -np.exp(1j * lam) * s],
        [np.exp(1j * phi) * s, np.exp(1j * (phi + lam)) * c]
    ], dtype=complex)


def rotation_axis(n: np.ndarray, theta: float) -> np.ndarray:
    """绕任意轴旋转: exp(-iθ n·σ/2)

    参数:
        n:     单位方向向量 (3,)
        theta: 旋转角

    返回:
        2×2 酉矩阵
    """
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    nx, ny, nz = n
    return np.cos(theta / 2) * I - 1j * np.sin(theta / 2) * (nx * X + ny * Y + nz * Z)


# ============================================================
# 两量子比特门 (4×4)
# ============================================================

CNOT = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)  # 控制=qubit0, 目标=qubit1

CNOT_rev = np.array([
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0]
], dtype=complex)  # 控制=qubit1, 目标=qubit0

CZ = np.diag([1, 1, 1, -1]).astype(complex)  # 控制 Z

SWAP = np.array([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]
], dtype=complex)


def controlled_U(U: np.ndarray) -> np.ndarray:
    """受控 U 门: |0⟩⟨0|⊗I + |1⟩⟨1|⊗U

    控制=qubit0, 目标=qubit1

    参数:
        U: 2×2 酉矩阵

    返回:
        4×4 矩阵
    """
    result = np.zeros((4, 4), dtype=complex)
    result[:2, :2] = np.eye(2)
    result[2:, 2:] = U
    return result


# ============================================================
# 三量子比特门 (8×8)
# ============================================================

def Toffoli() -> np.ndarray:
    """Toffoli 门 (CCNOT): 两个控制 qubit 翻转目标 qubit

    |a,b,c⟩ → |a,b,c⊕(a·b)⟩
    """
    gate = np.eye(8, dtype=complex)
    gate[6, 6] = 0
    gate[6, 7] = 1
    gate[7, 7] = 0
    gate[7, 6] = 1
    return gate


# ============================================================
# 通用门构建
# ============================================================

def expand_gate(gate: np.ndarray, target: int, n_qubits: int) -> np.ndarray:
    """将单 qubit 门扩展到 n-qubit 系统

    G → I ⊗ ... ⊗ I ⊗ G ⊗ I ⊗ ... ⊗ I
                     ↑ target

    参数:
        gate:      单 qubit 门 (2×2)
        target:    目标 qubit 索引 (0-based)
        n_qubits:  总 qubit 数

    返回:
        (2^n, 2^n) 矩阵
    """
    if gate.shape != (2, 2):
        raise ValueError(f"Gate must be 2×2, got {gate.shape}")
    result = np.array([[1]], dtype=complex)
    for i in range(n_qubits):
        if i == target:
            result = np.kron(result, gate)
        else:
            result = np.kron(result, I)
    return result


def expand_cnot(control: int, target: int, n_qubits: int) -> np.ndarray:
    """将 CNOT 扩展到 n-qubit 系统 (任意控制/目标)

    量子比特 0 对应张量积的最左因子 (二进制索引的 MSB)。

    参数:
        control: 控制 qubit 索引 (0 = 最左)
        target:  目标 qubit 索引
        n_qubits: 总 qubit 数

    返回:
        (2^n, 2^n) 矩阵
    """
    dim = 2**n_qubits
    gate = np.zeros((dim, dim), dtype=complex)

    # 将 qubit 索引映射到二进制位 (0=MSB, n-1=LSB)
    c_bit = n_qubits - 1 - control
    t_bit = n_qubits - 1 - target

    for i in range(dim):
        if (i >> c_bit) & 1:
            j = i ^ (1 << t_bit)
            gate[j, i] = 1.0
        else:
            gate[i, i] = 1.0
    return gate


# ============================================================
# 门诊断
# ============================================================

def is_unitary(gate: np.ndarray, tol: float = 1e-10) -> bool:
    """检查门是否幺正: U†U = I"""
    return np.allclose(gate @ gate.conj().T, np.eye(gate.shape[0]), atol=tol)


def is_hermitian(gate: np.ndarray, tol: float = 1e-10) -> bool:
    """检查门是否厄米"""
    return np.allclose(gate, gate.conj().T, atol=tol)


def gate_decompose(U: np.ndarray) -> dict:
    """将任意单 qubit 酉矩阵分解为 ZYZ 旋转

    U = e^{iα} Rz(φ) Ry(θ) Rz(λ)

    返回:
        {'alpha': α, 'theta': θ, 'phi': φ, 'lambda': λ}
    """
    det = np.linalg.det(U)
    alpha = np.angle(det) / 2

    # 提取 Rz(φ) Ry(θ) Rz(λ)
    V = U * np.exp(-1j * alpha)
    theta = 2 * np.arctan2(abs(V[1, 0]), abs(V[0, 0]))
    if abs(np.cos(theta / 2)) > 1e-10:
        phi = np.angle(V[1, 0] / np.sin(theta / 2))
    else:
        phi = 0
    lam = np.angle(V[1, 1]) - phi

    return {'alpha': alpha, 'theta': theta, 'phi': phi, 'lambda': lam}

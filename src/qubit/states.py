"""量子态模块 — 量子比特态的创建与操作

纯态:
    |0⟩, |1⟩, |+⟩, |-⟩, |+i⟩, |-i⟩
    任意叠加态 α|0⟩ + β|1⟩

多量子比特:
    |00⟩, GHZ 态, W 态, Bell 态

混合态:
    密度矩阵 ρ = Σ p_i |ψ_i⟩⟨ψ_i|

所有态均为 numpy 复数组 (列向量) 或矩阵。
"""

import numpy as np
from typing import Tuple, List, Optional, Union

# ============================================================
# 单量子比特态
# ============================================================

# 计算基态
ket0 = np.array([1.0 + 0j, 0.0 + 0j])       # |0⟩
ket1 = np.array([0.0 + 0j, 1.0 + 0j])       # |1⟩

# X 本征态
ket_plus  = np.array([1.0, 1.0]) / np.sqrt(2)   # |+⟩ = (|0⟩+|1⟩)/√2
ket_minus = np.array([1.0, -1.0]) / np.sqrt(2)  # |-⟩ = (|0⟩-|1⟩)/√2

# Y 本征态
ket_plus_i  = np.array([1.0, 1j]) / np.sqrt(2)   # |+i⟩ = (|0⟩+i|1⟩)/√2
ket_minus_i = np.array([1.0, -1j]) / np.sqrt(2)  # |-i⟩ = (|0⟩-i|1⟩)/√2

# 布洛赫球上的任意纯态
def ket(theta: float, phi: float = 0.0) -> np.ndarray:
    """布洛赫球参数化纯态

    |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩

    参数:
        theta: 极角 [0, π]
        phi:   方位角 [0, 2π]

    返回:
        shape (2,) 复向量
    """
    return np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])


def ket_from_bloch(x: float, y: float, z: float) -> np.ndarray:
    """从布洛赫球坐标构建纯态

    r = (x, y, z) 是单位球面上的点

    返回:
        shape (2,) 复向量
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    x, y, z = x/r, y/r, z/r
    theta = np.arccos(z)
    phi = np.arctan2(y, x)
    return ket(theta, phi)


# ============================================================
# 多量子比特态
# ============================================================

def tensor_product(*states: np.ndarray) -> np.ndarray:
    """多量子比特张量积

    |ψ₁⟩ ⊗ |ψ₂⟩ ⊗ ... ⊗ |ψₙ⟩

    用法:
        tensor_product(ket0, ket1)  → |01⟩
        tensor_product(ket_plus, ket0)  → |+0⟩
    """
    result = states[0]
    for s in states[1:]:
        result = np.kron(result, s)
    return result


def computational_basis_state(n: int, index: int) -> np.ndarray:
    """n 量子比特计算基态 |index⟩

    |0⟩ → index=0, |1⟩ → index=1, |2⟩ → index=2 (|10⟩)

    返回:
        shape (2^n,) 向量
    """
    dim = 2**n
    ket_vec = np.zeros(dim, dtype=complex)
    ket_vec[index] = 1.0
    return ket_vec


def bell_state(i: int = 0) -> np.ndarray:
    """Bell 态 (最大纠缠态)

    i=0: |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
    i=1: |Φ⁻⟩ = (|00⟩ - |11⟩)/√2
    i=2: |Ψ⁺⟩ = (|01⟩ + |10⟩)/√2
    i=3: |Ψ⁻⟩ = (|01⟩ - |10⟩)/√2

    返回:
        shape (4,) 向量
    """
    if i == 0:
        return np.array([1, 0, 0, 1]) / np.sqrt(2)
    elif i == 1:
        return np.array([1, 0, 0, -1]) / np.sqrt(2)
    elif i == 2:
        return np.array([0, 1, 1, 0]) / np.sqrt(2)
    elif i == 3:
        return np.array([0, 1, -1, 0]) / np.sqrt(2)
    raise ValueError(f"i must be 0-3, got {i}")


def ghz_state(n: int) -> np.ndarray:
    """n 量子比特 GHZ 态

    |GHZ⟩ = (|0...0⟩ + |1...1⟩)/√2

    返回:
        shape (2^n,) 向量
    """
    dim = 2**n
    ket_vec = np.zeros(dim, dtype=complex)
    ket_vec[0] = 1.0
    ket_vec[-1] = 1.0
    return ket_vec / np.sqrt(2)


def w_state(n: int) -> np.ndarray:
    """n 量子比特 W 态

    |W⟩ = (|10...0⟩ + |010...0⟩ + ... + |0...01⟩)/√n

    返回:
        shape (2^n,) 向量
    """
    dim = 2**n
    ket_vec = np.zeros(dim, dtype=complex)
    for i in range(n):
        ket_vec[2**i] = 1.0
    return ket_vec / np.sqrt(n)


# ============================================================
# 密度矩阵
# ============================================================

def density_matrix(ket_vec: np.ndarray) -> np.ndarray:
    """纯态密度矩阵: ρ = |ψ⟩⟨ψ|

    返回:
        shape (d, d) 矩阵
    """
    ket_vec = np.asarray(ket_vec).reshape(-1, 1)
    return ket_vec @ ket_vec.conj().T


def mixed_state(probabilities: List[float],
                states: List[np.ndarray]) -> np.ndarray:
    """混合态密度矩阵: ρ = Σ p_i |ψ_i⟩⟨ψ_i|

    返回:
        shape (d, d) 矩阵
    """
    dim = len(states[0])
    rho = np.zeros((dim, dim), dtype=complex)
    for p, ket_vec in zip(probabilities, states):
        rho += p * density_matrix(ket_vec)
    return rho


def partially_mixed(purity: float = 0.5) -> np.ndarray:
    """部分混合态: ρ = p|0⟩⟨0| + (1-p)|1⟩⟨1|

    purity ∈ [0.5, 1] (对角混合度)
    p = (1 + √(2*purity - 1)) / 2

    返回:
        shape (2, 2) 密度矩阵
    """
    p = (1 + np.sqrt(2 * purity - 1)) / 2
    return mixed_state([p, 1 - p], [ket0, ket1])


# ============================================================
# 测量与可观测量
# ============================================================

def measure(state: np.ndarray) -> Tuple[int, np.ndarray]:
    """模拟计算基测量 (随机坍缩)

    对于纯态 |ψ⟩ = Σ c_i |i⟩:
        P(i) = |c_i|²

    返回:
        (测量结果整数, 坍缩后的态向量)
    """
    if state.ndim == 2:
        # 密度矩阵 → 对角元作为概率
        probs = np.real(np.diag(state))
    else:
        probs = np.abs(state)**2
    probs = probs / probs.sum()
    outcome = np.random.choice(len(probs), p=probs)
    # 坍缩
    collapsed = np.zeros_like(state)
    collapsed[outcome] = 1.0
    if state.ndim == 2:
        collapsed = np.outer(collapsed, collapsed.conj())
    return outcome, collapsed


def measurement_probabilities(state: np.ndarray) -> np.ndarray:
    """计算基测量概率分布

    返回:
        P(i) = |⟨i|ψ⟩|² 数组
    """
    if state.ndim == 2:
        return np.real(np.diag(state))
    return np.abs(state)**2


def expectation(state: np.ndarray, operator: np.ndarray) -> complex:
    """期望值: ⟨ψ|O|ψ⟩ 或 Tr[ρO]

    参数:
        state:    纯态向量或密度矩阵
        operator: 厄米算符矩阵
    """
    if state.ndim == 1:
        return np.conj(state) @ operator @ state
    else:
        return np.trace(state @ operator)


# ============================================================
# 量子态诊断
# ============================================================

def fidelity(state1: np.ndarray, state2: np.ndarray) -> float:
    """保真度 F = |⟨ψ₁|ψ₂⟩| (纯态) 或 F = Tr[√√ρ σ √ρ]²

    纯态简化: F = |⟨ψ₁|ψ₂⟩|²
    """
    if state1.ndim == 1 and state2.ndim == 1:
        return float(abs(np.conj(state1) @ state2)**2)
    # 密度矩阵保真度 (简化版，仅对角化)
    eigvals, _ = np.linalg.eigh(state1)
    sqrt_rho = np.sqrt(np.maximum(eigvals, 0))
    # 近似的 fidelity
    return float(np.real(np.trace(np.diag(sqrt_rho) @ state2 @ np.diag(sqrt_rho))))


def purity(rho: np.ndarray) -> float:
    """纯度: γ = Tr[ρ²] ∈ [1/d, 1]

    γ=1 → 纯态, γ<1 → 混合态
    """
    return float(np.real(np.trace(rho @ rho)))


def entanglement_entropy(state: np.ndarray, partition: int = 1) -> float:
    """纠缠熵 (von Neumann entropy of reduced density matrix)

    对 bipartite 系统，约化掉 partition 个 qubit，计算 S = -Tr[ρ_A log ρ_A]

    参数:
        state:    纯态向量 shape (2^n,)
        partition: 约化掉的 qubit 数 (约化第一个 partition 个 qubit)

    返回:
        von Neumann 熵 (0 if separable, log(2) per ebit)
    """
    dim = len(state)
    n_qubits = int(np.log2(dim))
    if 2**n_qubits != dim:
        raise ValueError(f"State dimension {dim} must be power of 2")

    # 构建约化密度矩阵
    dim_a = 2**partition
    dim_b = dim // dim_a
    psi = state.reshape(dim_a, dim_b)
    rho_a = psi @ psi.conj().T

    # von Neumann 熵
    eigvals = np.linalg.eigvalsh(rho_a)
    eigvals = eigvals[eigvals > 1e-15]  # 去掉零
    return float(-np.sum(eigvals * np.log2(eigvals)))


def concurrence(state: np.ndarray) -> float:
    """两体纠缠度 Concurrence (仅对 2-qubit)

    对于纯态 |ψ⟩ = a|00⟩ + b|01⟩ + c|10⟩ + d|11⟩:
        C = 2|ad - bc|

    返回:
        [0, 1], 0=separable, 1=maximally entangled
    """
    if len(state) != 4:
        raise ValueError("Concurrence only defined for 2-qubit states")
    a, b, c, d = state
    return float(2 * abs(a * d - b * c))


def bloch_vector(state: np.ndarray) -> np.ndarray:
    """单 qubit 布洛赫矢量

    r = (⟨X⟩, ⟨Y⟩, ⟨Z⟩) = (2Re[α*β], 2Im[α*β], |α|²-|β|²)

    返回:
        shape (3,) 实向量, 纯态时 |r|=1
    """
    if state.ndim == 2:
        # 密度矩阵
        r = np.array([
            2 * np.real(state[0, 1]),
            2 * np.imag(state[1, 0]),
            np.real(state[0, 0] - state[1, 1])
        ])
    else:
        alpha, beta = state[0], state[1]
        r = np.array([
            2 * np.real(np.conj(alpha) * beta),
            2 * np.imag(np.conj(alpha) * beta),
            abs(alpha)**2 - abs(beta)**2
        ])
    return r

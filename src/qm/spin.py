"""自旋 / 量子比特模块

Pauli 矩阵、Bloch 球、旋转算符、密度矩阵构造。
"""

import numpy as np


# ═══════════════════════════════════════════════════════════
# Pauli 矩阵
# ═══════════════════════════════════════════════════════════

def sigma_x():
    """Pauli-X σₓ"""
    return np.array([[0, 1], [1, 0]], dtype=complex)


def sigma_y():
    """Pauli-Y σ_y"""
    return np.array([[0, -1j], [1j, 0]], dtype=complex)


def sigma_z():
    """Pauli-Z σ_z"""
    return np.array([[1, 0], [0, -1]], dtype=complex)


def pauli():
    """返回 (σₓ, σ_y, σ_z) 三矩阵"""
    return sigma_x(), sigma_y(), sigma_z()


# ═══════════════════════════════════════════════════════════
# 量子态
# ═══════════════════════════════════════════════════════════

def qubit_state(theta=0.0, phi=0.0):
    """Bloch 球上的纯态

    |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩

    Parameters
    ----------
    theta : float
        极角 [0, π]
    phi : float
        方位角 [0, 2π]

    Returns
    -------
    np.ndarray (2,)
    """
    return np.array([np.cos(theta/2), np.exp(1j*phi) * np.sin(theta/2)])


def bloch_vector_to_dm(rx, ry, rz):
    """从 Bloch 向量构造密度矩阵

    ρ = ½(I + rx σ_x + ry σ_y + rz σ_z)

    Parameters
    ----------
    rx, ry, rz : float
        Bloch 向量分量。|r| ≤ 1: 合法量子态, |r| = 1: 纯态

    Returns
    -------
    np.ndarray (2, 2)
    """
    I = np.eye(2, dtype=complex)
    return 0.5 * (I + rx * sigma_x() + ry * sigma_y() + rz * sigma_z())


def dm_to_bloch_vector(rho):
    """密度矩阵 → Bloch 向量

    r_i = Tr[ρ σ_i]
    """
    return (np.real(np.trace(rho @ sigma_x())),
            np.real(np.trace(rho @ sigma_y())),
            np.real(np.trace(rho @ sigma_z())))


def bloch_length(rho):
    """Bloch 向量长度 |r| = √(2 Tr[ρ²] - 1)

    纯态: |r| = 1, 混合态: |r| < 1, 完全混合: |r| = 0
    """
    r = dm_to_bloch_vector(rho)
    return np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)


# ═══════════════════════════════════════════════════════════
# 旋转算符
# ═══════════════════════════════════════════════════════════

def rotation_x(angle):
    """绕 x 轴旋转 R_x(θ) = exp(-iθσ_x/2)"""
    c = np.cos(angle / 2)
    s = np.sin(angle / 2)
    return np.array([[c, -1j*s], [-1j*s, c]], dtype=complex)


def rotation_y(angle):
    """绕 y 轴旋转 R_y(θ) = exp(-iθσ_y/2)"""
    c = np.cos(angle / 2)
    s = np.sin(angle / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rotation_z(angle):
    """绕 z 轴旋转 R_z(θ) = exp(-iθσ_z/2)"""
    return np.array([[np.exp(-1j*angle/2), 0],
                     [0, np.exp(1j*angle/2)]], dtype=complex)


def rotation(n, angle):
    """绕任意轴 n (单位向量) 旋转 R_n(θ) = exp(-i θ/2 n·σ)"""
    nx, ny, nz = n / np.linalg.norm(n)
    c = np.cos(angle / 2)
    s = np.sin(angle / 2)
    return np.array([
        [c - 1j*s*nz, -1j*s*(nx - 1j*ny)],
        [-1j*s*(nx + 1j*ny), c + 1j*s*nz]
    ], dtype=complex)


# ═══════════════════════════════════════════════════════════
# 常用操作
# ═══════════════════════════════════════════════════════════

def hadamard():
    """Hadamard 门 H = (|0⟩⟨0| + |0⟩⟨1| + |1⟩⟨0| - |1⟩⟨1|)/√2"""
    return np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def cnot():
    """CNOT 门 (4×4)"""
    return np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)


# ═══════════════════════════════════════════════════════════
# 测量 + 统计
# ═══════════════════════════════════════════════════════════

def measure(psi, basis='z'):
    """单量子比特测量

    Parameters
    ----------
    psi : np.ndarray (2,)
        纯态
    basis : str
        'x', 'y', 'z' — 测量基

    Returns
    -------
    (result, prob_up, prob_down)
    """
    if basis == 'z':
        prob0 = np.abs(psi[0])**2
    elif basis == 'x':
        proj = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
        prob0 = np.abs(np.dot(proj.conj(), psi))**2
    elif basis == 'y':
        proj = np.array([1/np.sqrt(2), 1j/np.sqrt(2)])
        prob0 = np.abs(np.dot(proj.conj(), psi))**2
    else:
        raise ValueError(f"Unknown basis: {basis}")

    result = 0 if np.random.random() < prob0 else 1
    return result, prob0, 1 - prob0

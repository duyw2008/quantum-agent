"""量子态 — QuTiP 风格量子态函数库"""

import numpy as np
import math
from .basis import FockBasis, get_basis


def fock(N: int, n: int = 0) -> np.ndarray:
    """Fock 态 |n⟩"""
    if n >= N:
        raise ValueError(f"n={n} >= N={N}")
    ket = np.zeros(N, dtype=complex)
    ket[n] = 1.0
    return ket


def fock_dm(N: int, n: int = 0) -> np.ndarray:
    """Fock 态密度矩阵 ρ = |n⟩⟨n|"""
    k = fock(N, n).reshape(-1, 1)
    return k @ k.conj().T


def coherent(N: int, alpha: complex) -> np.ndarray:
    """相干态 |α⟩"""
    ket = np.zeros(N, dtype=complex)
    if abs(alpha) < 1e-15:
        ket[0] = 1.0
        return ket
    norm = np.exp(-0.5 * abs(alpha)**2)
    fact = 1.0
    ap = 1.0 + 0j
    for n in range(N):
        if n > 0:
            fact *= n
            ap *= alpha
        ket[n] = norm * ap / np.sqrt(fact)
    return ket


def coherent_dm(N: int, alpha: complex) -> np.ndarray:
    """相干态密度矩阵"""
    k = coherent(N, alpha).reshape(-1, 1)
    return k @ k.conj().T


def squeezed(N: int, zeta: complex) -> np.ndarray:
    """压缩真空 |ζ⟩ = Ŝ(ζ)|0⟩"""
    r = abs(zeta)
    theta = np.angle(zeta) if r > 1e-15 else 0.0
    ket = np.zeros(N, dtype=complex)
    norm = 1.0 / np.sqrt(np.cosh(r))
    for m in range(N // 2):
        n = 2 * m
        factor = ((-np.exp(1j * theta) * np.tanh(r))**m *
                  math.sqrt(math.factorial(2 * m)) /
                  (2**m * math.factorial(m)))
        ket[n] = norm * factor
    return ket


def thermal_dm(N: int, n_th: float) -> np.ndarray:
    """热态密度矩阵 ρ_th"""
    rho = np.zeros((N, N), dtype=complex)
    for n in range(N):
        rho[n, n] = n_th**n / (n_th + 1)**(n + 1)
    return rho


def cat(N: int, alpha: complex, phi: float = 0.0) -> np.ndarray:
    """薛定谔猫态 |ψ⟩ ∝ |α⟩ + e^{iφ}|-α⟩"""
    psi_p = coherent(N, alpha)
    psi_m = coherent(N, -alpha)
    psi = psi_p + np.exp(1j * phi) * psi_m
    return psi / np.linalg.norm(psi)


def fidelity(psi1: np.ndarray, psi2: np.ndarray) -> float:
    """保真度 F = |⟨ψ₁|ψ₂⟩|²"""
    if psi1.ndim == 1 and psi2.ndim == 1:
        return float(abs(np.conj(psi1) @ psi2)**2)
    return float(abs(np.trace(psi1 @ psi2)))


def purity(rho: np.ndarray) -> float:
    """纯度 Tr[ρ²]"""
    return float(np.real(np.trace(rho @ rho)))


def photon_dist(state: np.ndarray) -> np.ndarray:
    """光子数分布 P(n)"""
    if state.ndim == 1:
        return np.abs(state)**2
    return np.real(np.diag(state))


def is_dm(state: np.ndarray) -> bool:
    """是否为密度矩阵"""
    return state.ndim == 2

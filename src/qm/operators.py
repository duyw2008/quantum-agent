"""算符工具 — 对易子、期望值、关联函数

核心函数:
    commutator(A, B)          — 对易子 [A, B]
    expect(oper, state)       — 期望值 ⟨O⟩
    variance(oper, state)     — 方差 ΔO²
    g2(state)                 — 二阶关联 g²(0)
    mandel_q(state)           — Mandel Q 参数
"""

import numpy as np
from .basis import FockBasis, get_basis


# ============================================================
# 对易子
# ============================================================

def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """对易子 [A, B] = AB - BA"""
    return A @ B - B @ A


def anti_commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """反对易子 {A, B} = AB + BA"""
    return A @ B + B @ A


# ============================================================
# 期望值
# ============================================================

def expect(oper: np.ndarray, state: np.ndarray) -> complex:
    """期望值 ⟨O⟩ = ⟨ψ|O|ψ⟩ 或 Tr[ρ O]"""
    if state.ndim == 1:
        return np.conj(state) @ oper @ state
    return np.trace(state @ oper)


def variance(oper: np.ndarray, state: np.ndarray) -> float:
    """方差 Var(O) = ⟨O²⟩ - ⟨O⟩²"""
    e1 = expect(oper, state)
    e2 = expect(oper @ oper, state)
    return float(np.real(e2 - abs(e1)**2))


# ============================================================
# 光子统计
# ============================================================

def mean_photon(state: np.ndarray, fb: FockBasis = None) -> float:
    """平均光子数 ⟨â†â⟩"""
    if fb is None:
        fb = get_basis(len(state))
    return float(np.real(expect(fb.n_op, state)))


def g2(state: np.ndarray, fb: FockBasis = None) -> float:
    """二阶关联函数 g²(0) = ⟨â†â†ââ⟩ / ⟨â†â⟩²"""
    if fb is None:
        fb = get_basis(len(state))
    n_mean = mean_photon(state, fb)
    if n_mean < 1e-15:
        return 0.0
    a2 = fb.a @ fb.a
    ad2 = fb.a_dag @ fb.a_dag
    return float(np.real(expect(ad2 @ a2, state) / n_mean**2))


def mandel_q(state: np.ndarray, fb: FockBasis = None) -> float:
    """Mandel Q 参数 = ⟨Δn²⟩/⟨n⟩ - 1"""
    if fb is None:
        fb = get_basis(len(state))
    n_mean = mean_photon(state, fb)
    if n_mean < 1e-15:
        return 0.0
    n2 = expect(fb.n_op @ fb.n_op, state)
    var_n = float(np.real(n2)) - n_mean**2
    return var_n / n_mean - 1.0


# ============================================================
# 矩阵属性
# ============================================================

def is_hermitian(A: np.ndarray, tol: float = 1e-10) -> bool:
    """检查厄米性 A = A†"""
    return np.allclose(A, A.conj().T, atol=tol)


def is_unitary(U: np.ndarray, tol: float = 1e-10) -> bool:
    """检查幺正性 U†U = I"""
    return np.allclose(U @ U.conj().T, np.eye(U.shape[0]), atol=tol)

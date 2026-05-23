"""量子动力学 — 时间演化与主方程

两种演化方式:
    1. Schrödinger 方程:  iℏ d|ψ⟩/dt = H|ψ⟩
    2. Lindblad 主方程:  dρ/dt = -i/ℏ[H,ρ] + Σ D[L_k]ρ

方法:
    sesolve()   — Schrödinger 方程 (矩阵指数, 适用于小 N)
    mesolve()   — Lindblad 主方程 (RK4)
    steadystate() — 稳态求解
"""

import numpy as np
from typing import List, Optional, Callable


# ============================================================
# Schrödinger 方程
# ============================================================

def sesolve(H: np.ndarray, psi0: np.ndarray, tlist: np.ndarray,
            e_ops: List[np.ndarray] = None,
            hbar: float = 1.0) -> dict:
    """Schrödinger 方程求解 (矩阵指数法)

    |ψ(t)⟩ = exp(-iHt/ℏ) |ψ(0)⟩

    参数:
        H:     哈密顿量 (N, N)
        psi0:  初始态向量 (N,)
        tlist: 时间点数组
        e_ops: 期望值算符列表
        hbar:  ℏ

    返回:
        {'times': tlist, 'states': [...], 'expect': {idx: [...]}}
    """
    N = H.shape[0]
    # 对角化 H
    eigvals, eigvecs = np.linalg.eigh(H)

    result = {'times': tlist, 'states': []}
    if e_ops:
        result['expect'] = {i: np.zeros(len(tlist), dtype=complex)
                            for i in range(len(e_ops))}

    for idx, t in enumerate(tlist):
        # 时间演化
        U_diag = np.exp(-1j * eigvals * t / hbar)
        psi_t = eigvecs @ (U_diag * (eigvecs.conj().T @ psi0))
        result['states'].append(psi_t)

        if e_ops:
            for i, op in enumerate(e_ops):
                result['expect'][i][idx] = np.conj(psi_t) @ op @ psi_t

    return result


# ============================================================
# Lindblad 主方程
# ============================================================

def lindblad_rhs(H: np.ndarray, rho: np.ndarray,
                 c_ops: List[np.ndarray] = None,
                 hbar: float = 1.0) -> np.ndarray:
    """Lindblad 主方程右侧 dρ/dt

    dρ/dt = -i/ℏ[H, ρ] + Σ_k (L_k ρ L_k† - ½{L_k†L_k, ρ})
    """
    drho = -1j / hbar * (H @ rho - rho @ H)
    if c_ops:
        for L in c_ops:
            LdL = L.conj().T @ L
            drho += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return drho


def mesolve(H: np.ndarray, rho0: np.ndarray, tlist: np.ndarray,
            c_ops: List[np.ndarray] = None,
            e_ops: List[np.ndarray] = None,
            hbar: float = 1.0) -> dict:
    """Lindblad 主方程数值求解 (RK4)

    参数:
        H:      哈密顿量 (N,N)
        rho0:   初始密度矩阵 (N,N)
        tlist:  时间点数组
        c_ops:  坍缩算符列表 [L1, L2, ...]
        e_ops:  期望值算符列表
        hbar:   ℏ

    返回:
        {'times': tlist, 'states': [...], 'expect': {idx: [...]}}
    """
    n_steps = len(tlist)
    rho = rho0.copy()

    result = {'times': tlist, 'states': []}
    if e_ops:
        result['expect'] = {i: np.zeros(n_steps, dtype=complex)
                            for i in range(len(e_ops))}

    for idx in range(n_steps):
        result['states'].append(rho.copy())
        if e_ops:
            for i, op in enumerate(e_ops):
                result['expect'][i][idx] = np.trace(rho @ op)

        if idx < n_steps - 1:
            dt = tlist[idx + 1] - tlist[idx]
            rho = _rk4(lambda r: lindblad_rhs(H, r, c_ops, hbar), rho, dt)

    return result


def _rk4(f: Callable, y: np.ndarray, dt: float) -> np.ndarray:
    """4 阶 Runge-Kutta 步进"""
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


# ============================================================
# 稳态求解
# ============================================================

def steadystate(H: np.ndarray, c_ops: List[np.ndarray] = None,
                hbar: float = 1.0) -> np.ndarray:
    """Lindblad 主方程稳态: dρ/dt = 0

    直接求解 Liouville 超算符的线性系统。
    """
    N = H.shape[0]
    M = N * N
    I = np.eye(N)

    # Liouville 超算符 L: d/dt vec(ρ) = L vec(ρ)
    L = -1j / hbar * (np.kron(H, I) - np.kron(I, H.T))

    if c_ops:
        for L_op in c_ops:
            LdL = L_op.conj().T @ L_op
            L += (np.kron(L_op, L_op.conj())
                  - 0.5 * np.kron(LdL, I)
                  - 0.5 * np.kron(I, LdL.T))

    # 迹归一化约束
    trace_vec = np.zeros(M)
    for i in range(N):
        trace_vec[i * N + i] = 1.0
    L[-1, :] = trace_vec

    b = np.zeros(M)
    b[-1] = 1.0

    vec_rho = np.linalg.solve(L, b)
    rho_ss = vec_rho.reshape(N, N)

    # 确保物理性
    rho_ss = 0.5 * (rho_ss + rho_ss.conj().T)
    eigvals, eigvecs = np.linalg.eigh(rho_ss)
    eigvals[eigvals < 0] = 0
    rho_ss = eigvecs @ np.diag(eigvals) @ eigvecs.conj().T
    rho_ss /= np.trace(rho_ss)

    return rho_ss

"""Fock 基与算符 — QuTiP 风格的核心算符

在截断的 N 维 Fock 空间中定义:
    â   — 湮灭算符:  â|n⟩ = √n |n-1⟩
    â†  — 产生算符:  â†|n⟩ = √(n+1) |n+1⟩
    N̂   — 数算符:    N̂|n⟩ = n|n⟩
    x̂   — 坐标算符:  x̂ = √(ℏ/2mω) (â + â†)
    p̂   — 动量算符:  p̂ = i√(mℏω/2) (â† - â)

用法:
    from src.qm import FockBasis
    fb = FockBasis(20)       # 20 维截断
    a = fb.a                 # 湮灭算符矩阵
    H = fb.hamiltonian()     # 谐振子哈密顿量
"""

import numpy as np
from typing import Tuple


class FockBasis:
    """Fock 空间基与基本算符

    所有算符表示为 (N, N) numpy 复矩阵。

    参数:
        N:     Fock 空间截断维度
        hbar:  约化普朗克常数 (默认 1.0)
        mass:  粒子质量 (默认 1.0)
        omega: 参考频率 (默认 1.0)
    """

    def __init__(self, N: int = 50, hbar: float = 1.0,
                 mass: float = 1.0, omega: float = 1.0):
        self.N = N
        self.hbar = hbar
        self.mass = mass
        self.omega = omega

        # 特征标度
        self.x0 = np.sqrt(hbar / (mass * omega))   # 特征长度
        self.p0 = np.sqrt(mass * hbar * omega)      # 特征动量

        # 缓存
        self._cache = {}

    # ---------- 基础算符 ----------

    @property
    def a(self) -> np.ndarray:
        """湮灭算符 â"""
        if 'a' not in self._cache:
            a_mat = np.zeros((self.N, self.N), dtype=complex)
            for n in range(1, self.N):
                a_mat[n - 1, n] = np.sqrt(n)
            self._cache['a'] = a_mat
        return self._cache['a']

    @property
    def a_dag(self) -> np.ndarray:
        """产生算符 â†"""
        if 'a_dag' not in self._cache:
            self._cache['a_dag'] = self.a.conj().T
        return self._cache['a_dag']

    @property
    def n_op(self) -> np.ndarray:
        """数算符 N̂ = â†â"""
        if 'n_op' not in self._cache:
            self._cache['n_op'] = np.diag(np.arange(self.N, dtype=float))
        return self._cache['n_op']

    @property
    def x(self) -> np.ndarray:
        """坐标算符 x̂"""
        if 'x' not in self._cache:
            self._cache['x'] = self.x0 / np.sqrt(2) * (self.a + self.a_dag)
        return self._cache['x']

    @property
    def p(self) -> np.ndarray:
        """动量算符 p̂"""
        if 'p' not in self._cache:
            self._cache['p'] = 1j * self.p0 / np.sqrt(2) * (self.a_dag - self.a)
        return self._cache['p']

    @property
    def I(self) -> np.ndarray:
        """单位矩阵"""
        return np.eye(self.N)

    @property
    def parity(self) -> np.ndarray:
        """宇称算符 Π = (-1)^{N̂}"""
        if 'parity' not in self._cache:
            self._cache['parity'] = np.diag([(-1)**n for n in range(self.N)])
        return self._cache['parity']

    # ---------- 位移算符 ----------

    def displacement(self, alpha: complex) -> np.ndarray:
        """位移算符 D(α) = exp(α â† - α* â)

        通过对角化计算矩阵指数。
        """
        X_mat = alpha * self.a_dag - np.conj(alpha) * self.a
        eigvals, eigvecs = np.linalg.eig(X_mat)
        return eigvecs @ np.diag(np.exp(eigvals)) @ np.linalg.inv(eigvecs)

    # ---------- 哈密顿量 ----------

    def hamiltonian(self, omega: float = None) -> np.ndarray:
        """谐振子哈密顿量 Ĥ = ℏω (â†â + ½)"""
        w = omega if omega is not None else self.omega
        return self.hbar * w * (self.n_op + 0.5 * self.I)

    def make_hamiltonian(self, V_diag: np.ndarray = None) -> np.ndarray:
        """构建一般哈密顿量 Ĥ = p̂²/2m + V(x̂)

        参数:
            V_diag: 势能在坐标基本征值上的对角元 (可选)
        """
        T = self.p @ self.p / (2 * self.mass)
        if V_diag is None:
            return T
        # V(x̂) 通过对角化 x̂ 得到坐标基 → 变换回 Fock 基
        x_vals, x_vecs = np.linalg.eigh(self.x)
        V_mat = x_vecs @ np.diag(V_diag) @ x_vecs.conj().T
        return T + V_mat


# ============================================================
# 全局默认实例
# ============================================================

_default_basis = None


def get_basis(N: int = None) -> FockBasis:
    """获取或创建默认 FockBasis 实例"""
    global _default_basis
    if N is not None or _default_basis is None:
        _default_basis = FockBasis(N or 50)
    return _default_basis

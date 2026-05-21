"""符号矩阵力学 — 基于 sympy 的算符代数

提供符号量子力学计算能力:
    - 算符符号定义 (x̂, p̂, â, â†)
    - 对易子代数自动计算
    - 算符的矩阵表示 (在数态基)
    - 角动量代数 (J², J_z, J_±)

使用前需要安装 sympy:
    pip install sympy
"""

import numpy as np
from typing import Optional, Tuple, List, Dict, Union

try:
    import sympy as sp
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False
    sp = None  # type: ignore


class SymbolicQuantum:
    """符号量子力学计算引擎

    在 sympy 符号框架下处理量子算符和代数。
    如果 sympy 未安装，实例化时会抛出 ImportError。

    使用示例:
        sq = SymbolicQuantum()
        comm = sq.commutator_relation('x', 'p')  # → iℏ
        sq.ladder_algebra('a', 'a_dag')           # → 1
    """

    def __init__(self):
        if not HAS_SYMPY:
            raise ImportError(
                "SymbolicQuantum 需要 sympy。请安装: pip install sympy"
            )
        # 基本符号
        self.hbar_sym = sp.Symbol('hbar', real=True, positive=True)
        self.m_sym = sp.Symbol('m', real=True, positive=True)
        self.omega_sym = sp.Symbol('omega', real=True, positive=True)
        self.t_sym = sp.Symbol('t', real=True)

    # ---------- 对易子代数 ----------

    @staticmethod
    def commutator(A: 'sp.Expr', B: 'sp.Expr') -> 'sp.Expr':
        """对易子: [A, B] = A*B - B*A"""
        return A * B - B * A

    @staticmethod
    def simplify_commutator(expr: 'sp.Expr') -> 'sp.Expr':
        """简化对易子表达式"""
        return sp.simplify(expr)

    def canonical_commutation(self) -> 'sp.Expr':
        """正则对易关系: [x̂, p̂] = iℏ"""
        x, p = sp.symbols('x_hat p_hat', commutative=False)
        return sp.I * self.hbar_sym

    # ---------- 算符矩阵表示（截断数态基）----------

    def ladder_matrices(self, N: int) -> Tuple['sp.Matrix', 'sp.Matrix']:
        """产生/湮灭算符的矩阵表示（符号形式）

        â|n⟩ = √n |n-1⟩
        â†|n⟩ = √(n+1) |n+1⟩

        返回 (a, a_dag) 作为 sympy Matrix
        """
        a = sp.zeros(N, N)
        a_dag = sp.zeros(N, N)
        for n in range(1, N):
            a[n - 1, n] = sp.sqrt(n)
            a_dag[n, n - 1] = sp.sqrt(n)
        return a, a_dag

    def position_matrix(self, N: int) -> 'sp.Matrix':
        """坐标算符 x̂ = √(ℏ/2mω) (â + â†), 符号形式"""
        a, ad = self.ladder_matrices(N)
        x0 = sp.sqrt(self.hbar_sym / (2 * self.m_sym * self.omega_sym))
        return x0 * (a + ad)

    def momentum_matrix(self, N: int) -> 'sp.Matrix':
        """动量算符 p̂ = i √(mℏω/2) (â† - â), 符号形式"""
        a, ad = self.ladder_matrices(N)
        p0 = sp.sqrt(self.m_sym * self.hbar_sym * self.omega_sym / 2)
        return sp.I * p0 * (ad - a)

    def number_matrix(self, N: int) -> 'sp.Matrix':
        """数算符 N̂ = â†â, 符号形式"""
        a, ad = self.ladder_matrices(N)
        return ad * a

    def harmonic_hamiltonian_matrix(self, N: int) -> 'sp.Matrix':
        """谐振子哈密顿量 Ĥ = ℏω(â†â + ½)"""
        return self.hbar_sym * self.omega_sym * (self.number_matrix(N) +
                                                  sp.eye(N) / 2)

    # ---------- 符号本征值问题 ----------

    def harmonic_energies_symbolic(self, n_max: int = 5) -> List['sp.Expr']:
        """谐振子解析能级: E_n = ℏω(n + ½)"""
        energies = []
        for n in range(n_max + 1):
            E_n = self.hbar_sym * self.omega_sym * (n + sp.Rational(1, 2))
            energies.append(E_n)
        return energies

    # ---------- 常用恒等式 ----------

    @staticmethod
    def baker_campbell_hausdorff(A: 'sp.Expr', B: 'sp.Expr',
                                  commutator_AB: 'sp.Expr') -> 'sp.Expr':
        r"""BCH 公式 (截断到二阶):
        eᴬ eᴮ = exp(A + B + ½[A,B])  (当 [A,[A,B]] = [B,[A,B]] = 0)
        """
        return sp.exp(A + B + commutator_AB / 2)

    @staticmethod
    def hadamard_lemma(A: 'sp.Expr', B: 'sp.Expr', commutator_AB: 'sp.Expr',
                       t: 'sp.Expr' = None) -> 'sp.Expr':
        r"""Hadamard 引理: exp(tA) B exp(-tA) = B + t[A,B] + ..."""
        return B + commutator_AB

    # ---------- 角动量代数 ----------

    def angular_momentum_commutators(self) -> Dict[str, 'sp.Expr']:
        r"""角动量对易关系 (SO(3)):
        [J_x, J_y] = iℏ J_z, [J_y, J_z] = iℏ J_x, [J_z, J_x] = iℏ J_y
        """
        Jx, Jy, Jz = sp.symbols('J_x J_y J_z', commutative=False)
        return {
            '[Jx, Jy]': sp.I * self.hbar_sym * Jz,
            '[Jy, Jz]': sp.I * self.hbar_sym * Jx,
            '[Jz, Jx]': sp.I * self.hbar_sym * Jy,
        }

    def spin_matrices(self, s=None) -> Tuple['sp.Matrix', ...]:
        r"""自旋算符 S_x, S_y, S_z 的矩阵表示
        
        对于自旋 ½: S_i = (ℏ/2) σ_i (Pauli 矩阵)
        返回 (Sx, Sy, Sz, S2)
        """
        if s is None:
            s = sp.Rational(1, 2)
        sx = sp.Matrix([[0, 1], [1, 0]]) * self.hbar_sym / 2
        sy = sp.Matrix([[0, -sp.I], [sp.I, 0]]) * self.hbar_sym / 2
        sz = sp.Matrix([[1, 0], [0, -1]]) * self.hbar_sym / 2
        s2 = sx**2 + sy**2 + sz**2
        return sx, sy, sz, s2


# ============================================================
# 便捷封装：解析 + 数值混合
# ============================================================

class MatrixMechanics:
    """矩阵力学统一接口

    结合数值矩阵表示和符号推导。

    用法:
        mm = MatrixMechanics(n_basis=50, hbar=1.0, mass=1.0, omega=2.0)
    """

    def __init__(self, n_basis: int = 50, hbar: float = 1.0,
                 mass: float = 1.0, omega: float = 1.0):
        from .operators import NumericOperatorSystem
        self.numeric = NumericOperatorSystem(n_basis, hbar, mass, omega)
        self._hbar = hbar
        self._mass = mass
        self._omega = omega
        self._N = n_basis

    @property
    def x(self) -> np.ndarray:
        return self.numeric.x

    @property
    def p(self) -> np.ndarray:
        return self.numeric.p

    @property
    def a(self) -> np.ndarray:
        return self.numeric.a

    @property
    def a_dag(self) -> np.ndarray:
        return self.numeric.a_dag

    @property
    def N_matrix(self) -> np.ndarray:
        return self.numeric.number

    @property
    def H_harmonic(self) -> np.ndarray:
        return self.numeric.harmonic_hamiltonian()

    def eigensolve(self, H: np.ndarray = None, k: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        if H is None:
            H = self.H_harmonic
        return self.numeric.eigensolve(H, k)

    def check_commutation(self, A: np.ndarray, B: np.ndarray) -> Dict[str, float]:
        return self.numeric.check_commutation(A, B)

    def commutator(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        from .operators import commutator
        return commutator(A, B)

    def expectation(self, O: np.ndarray, state: np.ndarray) -> complex:
        return self.numeric.expectation(O, state)

    def report(self) -> str:
        lines = [
            "=" * 50,
            "  Matrix Mechanics Report",
            "=" * 50,
            f"  Basis:  N = {self._N}  (number state truncation)",
            f"  ℏ = {self._hbar}, m = {self._mass}, ω = {self._omega}",
            "",
            "  --- Commutation Relations ---",
        ]
        xp = self.check_commutation(self.x, self.p)
        lines.append(f"  [x̂, p̂] norm: {xp['frobenius_norm']:.2e}")
        lines.append(f"  Expected:    iℏI (trace = iℏ·N = {1j * self._hbar * self._N})")
        aa = self.check_commutation(self.a, self.a_dag)
        lines.append(f"  [â, â†] norm: {aa['frobenius_norm']:.2e}")
        lines.append(f"  Expected:    I (Frobenius norm = √N = {np.sqrt(self._N):.2f})")
        x_herm = self.numeric.check_hermiticity(self.x)
        p_herm = self.numeric.check_hermiticity(self.p)
        lines.append(f"  x̂ hermitian: {x_herm['is_hermitian']}")
        lines.append(f"  p̂ hermitian: {p_herm['is_hermitian']}")
        energies, _ = self.eigensolve(k=5)
        lines.append("")
        lines.append("  --- Harmonic Oscillator Energies (ℏω units) ---")
        for n, E in enumerate(energies):
            expected = self._hbar * self._omega * (n + 0.5)
            lines.append(f"  E_{n} = {E:.6f}  (expected: {expected:.6f})")
        lines.append("=" * 50)
        return "\n".join(lines)

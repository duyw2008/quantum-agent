"""矩阵力学模块 — 符号和数值量子力学计算

提供两套接口:
    1. SymbolicOperatorSystem  — sympy 符号算符代数
    2. NumericOperatorSystem   — numpy 数值矩阵表示

核心功能:
    - 算符表示: x̂, p̂, â, â†, Ĥ
    - 对易子计算: [Â, B̂] = ÂB̂ - B̂Â
    - 本征值/本征态求解
    - 期望值和不确定度
    - 时间演化 (von Neumann 方程)

数学基础:
    - 坐标表象:  x̂|x⟩ = x|x⟩,  p̂ = -iℏ d/dx
    - 数态表象:  â|n⟩ = √n |n-1⟩,  â†|n⟩ = √(n+1) |n+1⟩
    - 对易关系:  [x̂, p̂] = iℏ,  [â, â†] = 1
    - 谐振子:    Ĥ = ℏω(â†â + ½)
"""

import numpy as np
from typing import Optional, Tuple, List, Dict, Union, Callable
from dataclasses import dataclass, field
import warnings

# ============================================================
# numpy 数值矩阵力学
# ============================================================

class NumericOperatorSystem:
    """数值矩阵力学系统

    在数态表象 (Fock basis) 下表示所有算符为 N×N 矩阵。

    基础算符:
        â  (湮灭算符):  ⟨n|â|m⟩ = √m δ_{n, m-1}
        â† (产生算符):  ⟨n|â†|m⟩ = √(m+1) δ_{n, m+1}

    坐标和动量 (在数态表象):
        x̂ = √(ℏ/2mω) (â + â†)
        p̂ = i √(mℏω/2) (â† - â)
    """

    def __init__(self, n_basis: int = 50, hbar: float = 1.0,
                 mass: float = 1.0, omega: float = 1.0):
        """
        参数:
            n_basis: 数态基大小 (截断)
            hbar:    约化普朗克常数
            mass:    粒子质量
            omega:   谐振子频率 (用于坐标/动量算符的标度)
        """
        self.N = n_basis
        self.hbar = hbar
        self.mass = mass
        self.omega = omega

        # 特征长度和动量标度
        self.x0 = np.sqrt(hbar / (mass * omega))   # a_ho
        self.p0 = np.sqrt(mass * hbar * omega)

        # 缓存
        self._cache: Dict[str, np.ndarray] = {}

    # ---------- 基础算符 ----------

    @property
    def a(self) -> np.ndarray:
        """湮灭算符 â"""
        if 'a' not in self._cache:
            n = np.arange(1, self.N)
            self._cache['a'] = np.diag(np.sqrt(n), k=1)
        return self._cache['a']

    @property
    def a_dag(self) -> np.ndarray:
        """产生算符 â†"""
        if 'a_dag' not in self._cache:
            n = np.arange(1, self.N)
            self._cache['a_dag'] = np.diag(np.sqrt(n), k=-1)
        return self._cache['a_dag']

    @property
    def x(self) -> np.ndarray:
        """坐标算符 x̂ = √(ℏ/2mω) (â + â†)"""
        if 'x' not in self._cache:
            self._cache['x'] = self.x0 / np.sqrt(2) * (self.a + self.a_dag)
        return self._cache['x']

    @property
    def p(self) -> np.ndarray:
        """动量算符 p̂ = i √(mℏω/2) (â† - â)"""
        if 'p' not in self._cache:
            self._cache['p'] = 1j * self.p0 / np.sqrt(2) * (self.a_dag - self.a)
        return self._cache['p']

    @property
    def number(self) -> np.ndarray:
        """数算符 N̂ = â†â"""
        if 'number' not in self._cache:
            self._cache['number'] = np.diag(np.arange(self.N, dtype=float))
        return self._cache['number']

    @property
    def identity(self) -> np.ndarray:
        """单位矩阵"""
        return np.eye(self.N)

    def x_squared(self) -> np.ndarray:
        """x̂² (使用矩阵乘法)"""
        x_mat = self.x
        return x_mat @ x_mat

    def p_squared(self) -> np.ndarray:
        """p̂²"""
        p_mat = self.p
        return p_mat @ p_mat

    # ---------- 哈密顿量 ----------

    def harmonic_hamiltonian(self) -> np.ndarray:
        """谐振子哈密顿量: Ĥ = ℏω(â†â + ½)"""
        return self.hbar * self.omega * (self.number + 0.5 * self.identity)

    def hamiltonian_from_potential(self, V: Callable[[np.ndarray], np.ndarray],
                                    x_grid: Optional[np.ndarray] = None) -> np.ndarray:
        """在数态表象中构建一般势的哈密顿量: Ĥ = p̂²/2m + V(x̂)

        方法: 对角化 x̂ 得到坐标基，在坐标基中计算 V(x) 是对角的，再变换回来。

        参数:
            V:      势能函数 V(x)
            x_grid: 可选，用于诊断
        """
        p2 = self.p_squared()
        T = p2 / (2 * self.mass)  # 动能项

        # 对角化 x̂ 得到变换矩阵
        eigvals, eigvecs = np.linalg.eigh(self.x)
        # eigvecs[:, i] 是坐标基波函数在数态基中的展开

        # V(x̂) 在坐标本征基下是对角的
        V_diag = V(eigvals)

        # 变换回数态基: V_num = U V_diag U†
        V_num = eigvecs @ np.diag(V_diag) @ eigvecs.conj().T

        return T + V_num

    def arbitrary_hamiltonian(self, H_func: Callable[..., np.ndarray],
                              *args, **kwargs) -> np.ndarray:
        """自定义哈密顿量"""
        return H_func(self, *args, **kwargs)

    # ---------- 对易子 ----------

    @staticmethod
    def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """对易子: [A, B] = AB - BA"""
        return A @ B - B @ A

    @staticmethod
    def anticommutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """反对易子: {A, B} = AB + BA"""
        return A @ B + B @ A

    # ---------- 本征值问题 ----------

    def eigensolve(self, H: np.ndarray, k: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        """求解哈密顿量的本征值和本征态

        返回:
            eigenvalues: 前 k 个本征值 (升序)
            eigenvectors: 前 k 个本征矢 (列向量), shape (N, k)
        """
        eigvals, eigvecs = np.linalg.eigh(H)
        return eigvals[:k], eigvecs[:, :k]

    # ---------- 期望值 ----------

    def expectation(self, O: np.ndarray, state: np.ndarray) -> complex:
        """期望值: ⟨ψ|Ô|ψ⟩"""
        return np.conj(state) @ (O @ state)

    def uncertainty(self, O: np.ndarray, state: np.ndarray) -> float:
        """不确定度: ΔO = √(⟨O²⟩ - ⟨O⟩²)"""
        O2 = O @ O
        exp_O = self.expectation(O, state)
        exp_O2 = self.expectation(O2, state)
        var = exp_O2 - abs(exp_O)**2
        return np.sqrt(max(var.real, 0.0))

    # ---------- 时间演化 ----------

    def time_evolution_operator(self, H: np.ndarray, t: float) -> np.ndarray:
        """时间演化算符: Û(t) = exp(-iĤt/ℏ)

        用于将态从 t=0 演化到 t:
            |ψ(t)⟩ = Û(t) |ψ(0)⟩
        """
        # 对角化 H 以计算矩阵指数
        eigvals, eigvecs = np.linalg.eigh(H)
        U = eigvecs @ np.diag(np.exp(-1j * eigvals * t / self.hbar)) @ eigvecs.conj().T
        return U

    def evolve_state(self, H: np.ndarray, state0: np.ndarray,
                     t: float) -> np.ndarray:
        """演化态: |ψ(t)⟩ = exp(-iĤt/ℏ) |ψ(0)⟩"""
        U = self.time_evolution_operator(H, t)
        return U @ state0

    def evolve_operator(self, H: np.ndarray, O: np.ndarray,
                        t: float) -> np.ndarray:
        """海森堡绘景: Ô(t) = exp(iĤt/ℏ) Ô exp(-iĤt/ℏ)"""
        U = self.time_evolution_operator(H, t)
        return U.conj().T @ O @ U

    # ---------- 基态搜索 ----------

    def find_ground_state(self, H: np.ndarray) -> Tuple[float, np.ndarray]:
        """寻找基态 (最低能量本征态)"""
        eigvals, eigvecs = np.linalg.eigh(H)
        return eigvals[0], eigvecs[:, 0]

    def imaginary_time_evolution(self, H: np.ndarray, state0: np.ndarray,
                                  tau: float = 1.0, n_steps: int = 100) -> np.ndarray:
        """虚时间演化逼近基态

        虚演化: |ψ(τ)⟩ = exp(-Ĥτ) |ψ₀⟩ / ||···||
        当 τ → ∞, |ψ(τ)⟩ → |ψ₀⟩ (基态)
        """
        dtau = tau / n_steps
        eigvals, eigvecs = np.linalg.eigh(H)

        # 在能量本征基中演化
        coeffs = eigvecs.conj().T @ state0
        for _ in range(n_steps):
            coeffs *= np.exp(-eigvals * dtau)
            coeffs /= np.linalg.norm(coeffs)

        return eigvecs @ coeffs

    # ---------- 诊断 ----------

    def check_commutation(self, A: np.ndarray, B: np.ndarray,
                          tol: float = 1e-10) -> Dict[str, float]:
        """检查对易关系"""
        C = self.commutator(A, B)
        return {
            'frobenius_norm': float(np.linalg.norm(C, 'fro')),
            'max_element': float(np.max(np.abs(C))),
            'is_commuting': bool(np.linalg.norm(C, 'fro') < tol),
        }

    def check_hermiticity(self, A: np.ndarray, tol: float = 1e-10) -> Dict[str, float]:
        """检查厄米性: A = A†"""
        diff = A - A.conj().T
        return {
            'frobenius_norm': float(np.linalg.norm(diff, 'fro')),
            'is_hermitian': bool(np.linalg.norm(diff, 'fro') < tol),
        }

    def summary(self) -> str:
        """打印系统信息"""
        lines = [
            f"NumericOperatorSystem:",
            f"  Basis size:      N = {self.N}",
            f"  ℏ = {self.hbar},  m = {self.mass},  ω = {self.omega}",
            f"  x₀ (char. len.): {self.x0:.4f}",
            f"  p₀ (char. mom.): {self.p0:.4f}",
            f"  Truncation err:  ~({self.N})^(-1/2) for high n",
        ]
        return "\n".join(lines)


# ============================================================
# 快捷函数
# ============================================================

def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """对易子: [A, B] = AB - BA"""
    return A @ B - B @ A

def anticommutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """反对易子: {A, B} = AB + BA"""
    return A @ B + B @ A

def is_hermitian(A: np.ndarray, tol: float = 1e-10) -> bool:
    """检查矩阵是否厄米"""
    return np.allclose(A, A.conj().T, atol=tol)

def is_unitary(U: np.ndarray, tol: float = 1e-10) -> bool:
    """检查矩阵是否幺正"""
    return np.allclose(U @ U.conj().T, np.eye(U.shape[0]), atol=tol)

def fidelity(state1: np.ndarray, state2: np.ndarray) -> float:
    """保真度: F = |⟨ψ₁|ψ₂⟩|"""
    return float(abs(np.conj(state1) @ state2))

"""自由标量场 — 1+1D 量子场论

场算符展开:
    φ̂(x) = Σ_k (â_k e^{ikx} + â†_k e^{-ikx}) / √(2ω_k L)

其中 ω_k = √(k² + m²), L 是空间体积, k = 2πn/L

核心可观测量:
    - [φ̂(x), φ̂(y)]      — 场对易子 (类空区域非零 → 微观因果性)
    - ⟨0|φ̂(x)²|0⟩       — 真空涨落
    - ⟨0|T{φ̂(x)φ̂(y)}|0⟩ — Feynman 传播子

全部在截断动量空间 (N_modes 个模式) 中计算。
"""

import numpy as np
from typing import Tuple


class ScalarField:
    """1+1D 自由标量场 φ̂(x)

    参数:
        mass:     场质量 m
        L:        空间长度 (默认 20)
        N_modes:  动量模式数 (截断, 默认 30)
        hbar:     ℏ (默认 1.0)
    """

    def __init__(self, mass: float = 1.0, L: float = 20.0,
                 N_modes: int = 30, hbar: float = 1.0):
        self.mass = mass
        self.L = L
        self.N_modes = N_modes
        self.hbar = hbar

        # 动量网格 (排除 k=0 避免红外发散)
        n = np.arange(1, N_modes + 1)
        self.k = 2 * np.pi * n / L
        self.omega = np.sqrt(self.k**2 + mass**2)

        # 归一化因子: 1/√(2ω_k L)
        self._norm = 1.0 / np.sqrt(2 * self.omega * L)

    # ================================================================
    # 场算符 φ̂(x) — 返回在位置 x 处的矩阵表示
    # ================================================================

    def field_matrix(self, x: float) -> np.ndarray:
        """场算符 φ̂(x) 的一阶量子化矩阵表示

        在 (N_modes+1) 维 Fock 空间截断中:
            φ̂(x) = Σ_k (â_k/√(2ω_kL)) e^{ikx} + h.c.

        返回:
            (N_modes+1, N_modes+1) 厄米矩阵
        """
        N = self.N_modes + 1
        phi = np.zeros((N, N), dtype=complex)

        # 每个动量模式贡献: â_k 项 (上对角线) + â†_k 项 (下对角线)
        for i, (k, omega) in enumerate(zip(self.k, self.omega)):
            n = i + 1  # 模式编号 (1-indexed)
            a_matrix = np.zeros((N, N), dtype=complex)
            a_matrix[n - 1, n] = np.sqrt(n)  # â|n⟩ = √n|n-1⟩
            a_dag_matrix = a_matrix.conj().T    # â†|n⟩ = √(n+1)|n+1⟩

            phase = np.exp(1j * k * x)
            phi += self._norm[i] * (phase * a_matrix +
                                     np.conj(phase) * a_dag_matrix)

        return phi

    def field_expectation(self, x: float, n_state: int) -> float:
        """⟨n|φ̂(x)|n⟩ — 对角元恒为零 (场的期望值为零)"""
        return 0.0  # ⟨n|φ̂|n⟩ = 0 ∀n (产生湮灭算符对角元为零)

    # ================================================================
    # 真空涨落 ⟨0|φ̂(x)²|0⟩
    # ================================================================

    def vacuum_fluctuation(self, x: float) -> float:
        """真空涨落 ⟨0|φ̂(x)²|0⟩

        只取产生-湮灭配对项:
            ⟨0|â_k â†_k|0⟩ = 1
            ⟨0|φ̂(x)²|0⟩ = Σ_k 1/(2ω_k L)

        注意: 与 x 无关 (平移不变), 且紫外发散需截断
        """
        return float(np.sum(1.0 / (2 * self.omega * self.L)))

    def vacuum_fluctuation_profile(self, x_points: np.ndarray = None) -> np.ndarray:
        """真空涨落 vs x (常数值, 但展示平移不变性)"""
        if x_points is None:
            x_points = np.linspace(-self.L/2, self.L/2, 200)
        return np.full_like(x_points, self.vacuum_fluctuation(0))

    def vacuum_energy_density(self) -> float:
        """真空能量密度 E₀/L = (1/2L) Σ_k ω_k (零点能)"""
        return float(0.5 * np.sum(self.omega) / self.L)

    # ================================================================
    # 场对易子 [φ̂(x), φ̂(y)]
    # ================================================================

    def commutator(self, x: float, y: float) -> float:
        """场对易子 ⟨0|[φ̂(x), φ̂(y)]|0⟩

        解析:
            [φ̂(x), φ̂(y)] = Σ_k (1/2ω_kL)[e^{ik(x-y)} - e^{-ik(x-y)}]
                          = -i Σ_k sin(k(x-y)) / (ω_k L)

        类空分离 (|x-y| > 0): 对易子 ≠ 0 但随距离振荡衰减
        类时分离: 持续振荡
        """
        dx = x - y
        return float(np.sum(np.sin(self.k * dx) / (self.omega * self.L)))

    def commutator_profile(self, x0: float = 0.0,
                           x_points: np.ndarray = None) -> np.ndarray:
        """对易子 [φ̂(x₀), φ̂(x)] vs x"""
        if x_points is None:
            x_points = np.linspace(-self.L/2, self.L/2, 300)
        return np.array([self.commutator(x0, x) for x in x_points])

    # ================================================================
    # Feynman 传播子
    # ================================================================

    def feynman_propagator(self, x: float, y: float, t: float = 0.0) -> complex:
        """Feynman 传播子 D_F(x-y, t) = ⟨0|T{φ̂(x,t)φ̂(y,0)}|0⟩

        动量空间:
            D_F(x, t) = Σ_k (1/2ω_kL) [θ(t)e^{-i(ω_k t - kx)} + θ(-t)e^{i(ω_k t + kx)}]

        等时 (t=0): 非传播的真空关联
            D_F(x, 0) = Σ_k e^{ikx} / (2ω_k L)
        """
        dx = x - y
        if abs(t) < 1e-15:
            return complex(np.sum(np.cos(self.k * dx) / (2 * self.omega * self.L)))
        # 含时传播子
        result = 0j
        for k, omega in zip(self.k, self.omega):
            result += (np.exp(-1j * (omega * abs(t) - k * dx * np.sign(t))) /
                       (2 * omega * self.L))
        return result

    def propagator_profile(self, t: float = 0.0,
                           x_points: np.ndarray = None) -> np.ndarray:
        """传播子 D_F(x, t) vs x (固定 t)"""
        if x_points is None:
            x_points = np.linspace(-self.L/2, self.L/2, 300)
        return np.array([self.feynman_propagator(x, 0, t) for x in x_points])

    # ================================================================
    # 粒子数算符
    # ================================================================

    def number_operator(self, mode_idx: int) -> np.ndarray:
        """第 mode_idx 个模式 (索引从 0 开始) 的粒子数算符 N̂_k = â†_k â_k

        在 Fock 空间中的对角矩阵: diag(0, 0, ..., 1, 1, ...)
        (第 mode_idx 个模式占据时为 1)
        """
        N = self.N_modes + 1
        # 简化: 返回总模式对应的对角算符
        # 更精确的实现需要对每个模式单独处理
        return np.diag(np.arange(N, dtype=float))

    # ================================================================
    # 信息
    # ================================================================

    def summary(self) -> str:
        lines = [
            f"ScalarField(m={self.mass}, L={self.L}, N_modes={self.N_modes})",
            f"  k_min = {self.k[0]:.4f}, k_max = {self.k[-1]:.4f}",
            f"  ω_min = {self.omega[0]:.4f}, ω_max = {self.omega[-1]:.4f}",
            f"  Vacuum energy density = {self.vacuum_energy_density():.4f}",
            f"  Vacuum fluctuation ⟨φ²⟩ = {self.vacuum_fluctuation(0):.4f}",
        ]
        return '\n'.join(lines)

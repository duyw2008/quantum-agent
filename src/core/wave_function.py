"""波函数表示模块

WaveFunction 封装了一维波函数的所有信息:
    - ψ(x): 复值波函数
    - ρ(x) = |ψ(x)|²: 概率密度
    - φ(x) = arg(ψ(x)): 相位
    - 归一化、内积、期望值计算等
"""

import numpy as np
from typing import Optional, Tuple, Callable
from dataclasses import dataclass

# numpy 2.x compatibility
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid


@dataclass
class Grid:
    """空间网格"""
    x_min: float
    x_max: float
    n_points: int

    @property
    def dx(self) -> float:
        return (self.x_max - self.x_min) / (self.n_points - 1)

    @property
    def x(self) -> np.ndarray:
        return np.linspace(self.x_min, self.x_max, self.n_points)

    @property
    def k(self) -> np.ndarray:
        """动量空间网格 (FFT 顺序)"""
        return 2 * np.pi * np.fft.fftfreq(self.n_points, self.dx)

    @property
    def dk(self) -> float:
        return 2 * np.pi / (self.x_max - self.x_min)

    def __repr__(self):
        return f"Grid(x=[{self.x_min}, {self.x_max}], N={self.n_points}, dx={self.dx:.4f})"


class WaveFunction:
    """一维波函数 ψ(x)

    属性:
        psi:  复值波函数数组 (shape: n_points,)
        grid: 空间网格
        t:    当前时间

    使用:
        psi = WaveFunction(grid)
        psi.set_gaussian(x0=0, p0=5.0, sigma=0.5)  # 高斯波包
        psi.set_eigenstate(V, n=0)                   # 定态
    """

    def __init__(self, grid: Grid, name: str = "ψ"):
        self.grid = grid
        self.name = name
        self.t: float = 0.0
        self._psi: np.ndarray = np.zeros(grid.n_points, dtype=np.complex128)

    # ---------- 属性 ----------

    @property
    def psi(self) -> np.ndarray:
        return self._psi

    @psi.setter
    def psi(self, value: np.ndarray):
        if len(value) != self.grid.n_points:
            raise ValueError(f"波函数长度 ({len(value)}) 与网格点数 ({self.grid.n_points}) 不匹配")
        self._psi = np.asarray(value, dtype=np.complex128)

    @property
    def probability_density(self) -> np.ndarray:
        """概率密度 ρ(x) = |ψ(x)|²"""
        return np.abs(self._psi)**2

    @property
    def phase(self) -> np.ndarray:
        """相位 φ(x) = arg(ψ(x))"""
        return np.angle(self._psi)

    @property
    def norm(self) -> float:
        """L² 范数: ∫|ψ|² dx"""
        return np.trapz(self.probability_density, self.grid.x)

    @property
    def is_normalized(self, tol: float = 1e-12) -> bool:
        """检查波函数是否已归一化"""
        return abs(self.norm - 1.0) < tol

    @property
    def psi_k(self) -> np.ndarray:
        """动量空间波函数 ψ̃(k) = F[ψ(x)]"""
        return np.fft.fft(self._psi)

    @property
    def momentum_density(self) -> np.ndarray:
        """动量空间概率密度 |ψ̃(k)|²"""
        return np.abs(self.psi_k)**2

    # ---------- 初始化 ----------

    def set_gaussian(self, x0: float = 0.0, p0: float = 0.0,
                     sigma: float = 0.5):
        """初始化为高斯波包

        ψ(x) = (1/(2πσ²)^(1/4)) exp(-(x - x0)² / 4σ²) exp(i p₀ x / ℏ)

        这是最小不确定度波包: Δx Δp = ℏ/2
        """
        x = self.grid.x
        norm_factor = 1.0 / (2 * np.pi * sigma**2)**(0.25)
        self._psi = (norm_factor *
                     np.exp(-(x - x0)**2 / (4 * sigma**2)) *
                     np.exp(1j * p0 * x))  # 使用 ℏ=1 单位
        self.t = 0.0

    def set_plane_wave(self, k0: float = 0.0, envelope_sigma: float = None):
        """初始化为平面波（可选高斯包络）

        ψ(x) ∝ exp(i k₀ x)，若 envelope_sigma 不为 None 则加高斯包络
        """
        x = self.grid.x
        self._psi = np.exp(1j * k0 * x)
        if envelope_sigma is not None:
            self._psi *= np.exp(-x**2 / (4 * envelope_sigma**2))
        self.normalize()

    def set_eigenstate(self, V, n: int = 0, hbar: float = 1.0,
                       mass: float = 1.0):
        """用数值对角化计算并设置第 n 个本征态

        构建哈密顿量矩阵 H = -ℏ²/2m d²/dx² + V(x)
        对角化得到本征值和本征态
        """
        from scipy.linalg import eigh_tridiagonal
        from scipy.sparse.linalg import eigsh

        N = self.grid.n_points
        dx = self.grid.dx
        x = self.grid.x

        # 动能算符 (三对角近似，二阶中心差分)
        ke_diag = np.ones(N) * hbar**2 / (mass * dx**2)
        ke_off_diag = np.ones(N - 1) * (-hbar**2 / (2 * mass * dx**2))

        # 势能
        Vx = V(x) if callable(V) else V

        # 三对角哈密顿量
        diag = ke_diag + Vx

        # 计算前 (n+1) 个最低本征态
        n_eigs = min(n + 3, N - 2)
        eigenvalues, eigenvectors = eigh_tridiagonal(diag, ke_off_diag,
                                                      select='i',
                                                      select_range=(0, n_eigs - 1))

        self._psi = eigenvectors[:, n]
        self.normalize()
        self._eigenvalue = eigenvalues[n]
        return eigenvalues[n]

    # ---------- 操作 ----------

    def normalize(self):
        """归一化: ψ → ψ / √∫|ψ|²dx"""
        nrm = self.norm
        if nrm > 0:
            self._psi /= np.sqrt(nrm)

    def inner(self, other: 'WaveFunction') -> complex:
        """内积: ⟨ψ|φ⟩ = ∫ ψ*(x) φ(x) dx"""
        return np.trapz(np.conj(self._psi) * other._psi, self.grid.x)

    def overlap(self, other: 'WaveFunction') -> float:
        """重叠积分 (保真度): |⟨ψ|φ⟩|²"""
        return abs(self.inner(other))**2

    # ---------- 期望值 ----------

    def expectation_x(self) -> float:
        """⟨x⟩ = ∫ x |ψ(x)|² dx"""
        return np.trapz(self.grid.x * self.probability_density, self.grid.x)

    def expectation_x2(self) -> float:
        """⟨x²⟩"""
        return np.trapz(self.grid.x**2 * self.probability_density, self.grid.x)

    def expectation_p(self) -> float:
        """⟨p⟩ = -iℏ ∫ ψ* ∂ψ/∂x dx   (使用动量空间更快)"""
        k = self.grid.k
        psi_k = self.psi_k
        prob_k = np.abs(psi_k)**2
        p_expect = np.trapz(k * prob_k, k)
        return p_expect / np.trapz(prob_k, k)  # ℏ=1, p = k

    def expectation_p2(self) -> float:
        """⟨p²⟩ = ∫ k² |ψ̃(k)|² dk / ∫ |ψ̃(k)|² dk"""
        k = self.grid.k
        psi_k = self.psi_k
        prob_k = np.abs(psi_k)**2
        norm = np.trapz(prob_k, k)
        return np.trapz(k**2 * prob_k, k) / norm

    def expectation_kinetic_energy(self, mass: float = 1.0,
                                   hbar: float = 1.0) -> float:
        """动能期望值: ⟨T⟩ = ⟨p²⟩ / 2m"""
        return self.expectation_p2() * hbar**2 / (2 * mass)

    def expectation_potential_energy(self, V, mass: float = 1.0) -> float:
        """势能期望值: ⟨V⟩ = ∫ V(x) |ψ(x)|² dx"""
        Vx = V(self.grid.x) if callable(V) else V
        return np.trapz(Vx * self.probability_density, self.grid.x)

    def uncertainty_x(self) -> float:
        """位置不确定度: Δx = √(⟨x²⟩ - ⟨x⟩²)"""
        return np.sqrt(self.expectation_x2() - self.expectation_x()**2)

    def uncertainty_p(self) -> float:
        """动量不确定度: Δp = √(⟨p²⟩ - ⟨p⟩²)"""
        return np.sqrt(self.expectation_p2() - self.expectation_p()**2)

    def uncertainty_product(self) -> float:
        """不确定度乘积: Δx Δp"""
        return self.uncertainty_x() * self.uncertainty_p()

    # ---------- 保存与复制 ----------

    def copy(self) -> 'WaveFunction':
        """深拷贝波函数"""
        new_wf = WaveFunction(self.grid, self.name)
        new_wf._psi = self._psi.copy()
        new_wf.t = self.t
        return new_wf

    def __copy__(self):
        return self.copy()

    # ---------- 表示 ----------

    def __repr__(self):
        return (f"WaveFunction({self.name}, N={self.grid.n_points}, "
                f"t={self.t:.3f}, norm={self.norm:.6f})")

    def summary(self) -> str:
        """返回波函数的状态摘要"""
        lines = [
            f"WaveFunction: {self.name}",
            f"  Grid:       {self.grid}",
            f"  Time:       t = {self.t:.4f}",
            f"  Norm:       ∫|ψ|²dx = {self.norm:.10f}",
            f"  ⟨x⟩:        {self.expectation_x():.6f} ± {self.uncertainty_x():.6f}",
            f"  ⟨p⟩:        {self.expectation_p():.6f} ± {self.uncertainty_p():.6f}",
            f"  Δx·Δp:      {self.uncertainty_product():.6f} (≥ 0.5 by HUP)",
        ]
        return "\n".join(lines)

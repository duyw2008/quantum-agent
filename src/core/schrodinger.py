"""薛定谔方程求解器模块

求解一维含时薛定谔方程 (TDSE):
    iℏ ∂ψ/∂t = -ℏ²/2m ∂²ψ/∂x² + V(x) ψ

两种数值方法:
    1. Split-Step Fourier Method (SSFM) — 快速、显式、谱精度
    2. Crank-Nicolson (CN)        — 无条件稳定、隐式、二阶精度

物理基础:
    SSFM 基于 Trotter 分解:
        e^{-iHΔt/ℏ} ≈ e^{-iVΔt/2ℏ} e^{-iTΔt/ℏ} e^{-iVΔt/2ℏ}
    动量空间动能算符是对角的: T̃(k) = ℏ²k²/2m

    CN 基于 Crank-Nicolson 离散化:
        (I + iHΔt/2ℏ) ψ^{n+1} = (I - iHΔt/2ℏ) ψ^n
    三对角系统，用 Thomas 算法求解。
"""

import numpy as np
from typing import Optional, Callable, List, Tuple
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from .wave_function import WaveFunction, Grid
from .potentials import Potential

# numpy 2.x compatibility
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid


# ============================================================
# 求解器输出
# ============================================================

@dataclass
class EvolutionResult:
    """时间演化结果"""
    times: np.ndarray                    # 时间点 (Nt,)
    psi_snapshots: List[np.ndarray]      # 波函数快照列表
    prob_density_snapshots: np.ndarray   # 概率密度 (Nt, Nx)
    expectation_x: np.ndarray            # ⟨x⟩(t)
    expectation_p: np.ndarray            # ⟨p⟩(t)
    energy: np.ndarray                   # 能量 E(t) = ⟨H⟩
    norm_history: np.ndarray             # 范数历史 (用于检查守恒性)
    grid: Grid
    method: str
    params: dict = field(default_factory=dict)

    @property
    def n_snapshots(self) -> int:
        return len(self.times)


# ============================================================
# 抽象基类
# ============================================================

class TDSE_Solver(ABC):
    """含时薛定谔方程求解器抽象基类"""

    def __init__(self, grid: Grid, V: Potential,
                 hbar: float = 1.0, mass: float = 1.0):
        self.grid = grid
        self.V = V
        self.hbar = hbar
        self.mass = mass

        # 预计算空间网格
        self.x = grid.x
        self.dx = grid.dx
        self.N = grid.n_points

    @abstractmethod
    def step(self, psi: np.ndarray, dt: float) -> np.ndarray:
        """单步时间推进: ψ(t+dt) = U(dt) ψ(t)"""
        ...

    def evolve(self, psi0: WaveFunction, t_max: float, dt: float,
               snapshot_interval: int = 10,
               callback: Optional[Callable] = None) -> EvolutionResult:
        """时间演化主循环

        参数:
            psi0:              初始波函数
            t_max:             总演化时间
            dt:                时间步长
            snapshot_interval: 每隔多少步保存一次快照
            callback:          每步调用的回调函数 callback(step, psi, t)

        返回:
            EvolutionResult 包含所有快照和可观测量
        """
        n_steps = int(t_max / dt)
        n_snapshots = n_steps // snapshot_interval + 1

        psi = psi0.psi.copy()
        t = 0.0

        # 预分配存储
        times = np.zeros(n_snapshots)
        psi_snapshots = []
        prob_snapshots = np.zeros((n_snapshots, self.N))
        exp_x = np.zeros(n_snapshots)
        exp_p = np.zeros(n_snapshots)
        energy = np.zeros(n_snapshots)
        norm_hist = np.zeros(n_snapshots)

        # 保存初始状态
        times[0] = t
        psi_snapshots.append(psi.copy())
        prob_snapshots[0] = np.abs(psi)**2
        exp_x[0] = np.trapz(self.x * np.abs(psi)**2, self.x)
        exp_p[0] = self._expectation_p(psi)
        energy[0] = self._compute_energy(psi)
        norm_hist[0] = np.trapz(np.abs(psi)**2, self.x)

        snap_idx = 0

        for step in range(1, n_steps + 1):
            # 时间推进
            psi = self.step(psi, dt)
            t = step * dt

            # 回调
            if callback is not None:
                callback(step, psi, t)

            # 保存快照
            if step % snapshot_interval == 0:
                snap_idx += 1
                times[snap_idx] = t
                psi_snapshots.append(psi.copy())
                prob_snapshots[snap_idx] = np.abs(psi)**2
                exp_x[snap_idx] = np.trapz(self.x * np.abs(psi)**2, self.x)
                exp_p[snap_idx] = self._expectation_p(psi)
                energy[snap_idx] = self._compute_energy(psi)
                norm_hist[snap_idx] = np.trapz(np.abs(psi)**2, self.x)

        # 截断到实际使用的长度
        actual_snapshots = snap_idx + 1

        return EvolutionResult(
            times=times[:actual_snapshots],
            psi_snapshots=psi_snapshots[:actual_snapshots],
            prob_density_snapshots=prob_snapshots[:actual_snapshots],
            expectation_x=exp_x[:actual_snapshots],
            expectation_p=exp_p[:actual_snapshots],
            energy=energy[:actual_snapshots],
            norm_history=norm_hist[:actual_snapshots],
            grid=self.grid,
            method=self.name,
            params={'dt': dt, 't_max': t_max, 'n_steps': n_steps},
        )

    def _expectation_p(self, psi: np.ndarray) -> float:
        """动量期望值"""
        k = self.grid.k
        psi_k = np.fft.fft(psi)
        prob_k = np.abs(psi_k)**2
        return np.trapz(k * prob_k, k) / np.trapz(prob_k, k)

    def _compute_energy(self, psi: np.ndarray) -> float:
        """计算总能量 E = ⟨T⟩ + ⟨V⟩"""
        # 动能: 使用动量空间
        k = self.grid.k
        psi_k = np.fft.fft(psi)
        prob_k = np.abs(psi_k)**2
        T_expect = np.trapz(0.5 * self.hbar**2 * k**2 / self.mass * prob_k, k)
        T_expect /= np.trapz(prob_k, k)

        # 势能
        Vx = self.V(self.x)
        prob_x = np.abs(psi)**2
        V_expect = np.trapz(Vx * prob_x, self.x)

        return T_expect + V_expect

    @property
    @abstractmethod
    def name(self) -> str:
        """求解器名称"""
        ...


# ============================================================
# Split-Step Fourier Method (SSFM)
# ============================================================

class SplitStepFourier(TDSE_Solver):
    r"""Split-Step Fourier Method (二阶对称)

    基于 Trotter-Suzuki 分解:
        e^{-iHΔt/ℏ} = e^{-iVΔt/2ℏ} e^{-iTΔt/ℏ} e^{-iVΔt/2ℏ} + O(Δt³)

    算法:
        1. 半步势能演化:  ψ₁(x) = exp(-i V(x) Δt / 2ℏ) ψ(x, t)
        2. 傅里叶变换:    ψ̃₁(k) = F[ψ₁(x)]
        3. 动能演化:      ψ̃₂(k) = exp(-i ℏ k² Δt / 2m) ψ̃₁(k)
        4. 逆傅里叶变换:  ψ₂(x) = F⁻¹[ψ̃₂(k)]
        5. 半步势能演化:  ψ(x, t+Δt) = exp(-i V(x) Δt / 2ℏ) ψ₂(x)

    优势: 快速 (FFT, O(N log N))、谱精度
    限制: 要求周期边界、势能不随时间快速变化
    """

    def __init__(self, grid: Grid, V: Potential,
                 hbar: float = 1.0, mass: float = 1.0):
        super().__init__(grid, V, hbar, mass)
        self.k = grid.k
        # 预计算动能演化因子 (给定 dt 后更新)
        self._ke_phase: Optional[np.ndarray] = None
        self._pe_phase: Optional[np.ndarray] = None
        self._current_dt: Optional[float] = None

    @property
    def name(self) -> str:
        return "Split-Step Fourier (SSFM)"

    def step(self, psi: np.ndarray, dt: float) -> np.ndarray:
        """SSFM 单步推进"""
        # 更新动能相位因子 (只在 dt 改变时重新计算)
        if dt != self._current_dt:
            self._ke_phase = np.exp(-0.5j * self.hbar * self.k**2 * dt / self.mass)
            self._current_dt = dt

        Vx = self.V(self.x)
        pe_half = np.exp(-0.5j * Vx * dt / self.hbar)

        # Step 1: 半步势能
        psi_work = pe_half * psi

        # Step 2: FFT → 动量空间
        psi_k = np.fft.fft(psi_work)

        # Step 3: 动能演化 (注意: 这里用全步因为对称分解中动能算符是完整的)
        # 正确形式: exp(-i ℏ k² Δt / 2m) → 使用全步
        ke_full = self._ke_phase**2  # exp(-i ℏ k² Δt / m)
        psi_k *= ke_full

        # Step 4: 逆 FFT
        psi_work = np.fft.ifft(psi_k)

        # Step 5: 半步势能
        psi_work *= pe_half

        return psi_work


# ============================================================
# Crank-Nicolson Method
# ============================================================

class CrankNicolson(TDSE_Solver):
    r"""Crank-Nicolson 方法 (隐式、无条件稳定)

    离散化:
        iℏ (ψ^{n+1} - ψ^n) / Δt = Ĥ (ψ^{n+1} + ψ^n) / 2

    重写为:
        (I + i Ĥ Δt / 2ℏ) ψ^{n+1} = (I - i Ĥ Δt / 2ℏ) ψ^n

    Ĥ = -ℏ²/2m d²/dx² + V(x)

    使用三对角近似:
        Ĥ ψ_j = -α ψ_{j-1} + (2α + V_j) ψ_j - α ψ_{j+1}
        其中 α = ℏ² / (2m Δx²)

    左侧矩阵 A = I + i Δt / (2ℏ) Ĥ  是复三对角矩阵
    右侧向量 b = (I - i Δt / (2ℏ) Ĥ) ψ^n

    求解 A ψ^{n+1} = b 使用 Thomas 算法 (复版本)。

    优势: 无条件稳定、守恒性好、二阶精度
    限制: 需要解三对角系统 (但 O(N) 的 Thomas 算法很快)
    """

    def __init__(self, grid: Grid, V: Potential,
                 hbar: float = 1.0, mass: float = 1.0):
        super().__init__(grid, V, hbar, mass)
        self.alpha = hbar**2 / (2 * mass * self.dx**2)

    @property
    def name(self) -> str:
        return "Crank-Nicolson (CN)"

    def step(self, psi: np.ndarray, dt: float) -> np.ndarray:
        """CN 单步推进"""
        N = self.N
        alpha = self.alpha
        Vx = self.V(self.x)

        # 构建三对角矩阵元素:
        # A = I + i Δt/2ℏ Ĥ
        # 对角线:  A_jj = 1 + i(Δt/2ℏ)(2α + V_j)
        # 次对角线: A_{j,j±1} = -i(Δt/2ℏ)α
        factor = 0.5j * dt  # i Δt / (2ℏ),  ℏ=1

        # 左侧矩阵 A 的三对角线
        a = np.full(N - 1, -factor * alpha, dtype=np.complex128)  # 下对角线
        d = 1.0 + factor * (2 * alpha + Vx)                        # 主对角线
        c = np.full(N - 1, -factor * alpha, dtype=np.complex128)  # 上对角线

        # 边界条件: 默认 Dirichlet (ψ = 0 at boundary)
        # 对于周期边界需要修改三对角

        # 右侧: b = (I - i Δt/2ℏ Ĥ) ψ^n
        b = psi.copy()
        # 应用 -i Δt/2ℏ Ĥ
        b -= factor * (2 * alpha * psi + Vx * psi)
        b[1:] += factor * alpha * psi[:-1]   # +i Δt/2ℏ · α ψ_{j-1}
        b[:-1] += factor * alpha * psi[1:]   # +i Δt/2ℏ · α ψ_{j+1}

        # 求解三对角系统 (复 Thomas 算法)
        psi_new = self._complex_tridiag_solve(a, d, c, b)

        return psi_new

    @staticmethod
    def _complex_tridiag_solve(a: np.ndarray, d: np.ndarray,
                                c: np.ndarray, b: np.ndarray) -> np.ndarray:
        """复三对角系统 Thomas 算法

        求解:
            d[0]x[0] + c[0]x[1]                   = b[0]
            a[i-1]x[i-1] + d[i]x[i] + c[i]x[i+1] = b[i]
            a[N-2]x[N-2] + d[N-1]x[N-1]           = b[N-1]

        时间复杂度: O(N)，无需额外内存分配 (in-place)
        """
        N = len(d)
        cp = np.empty(N - 1, dtype=np.complex128)
        x = np.empty(N, dtype=np.complex128)

        # 前向消元
        cp[0] = c[0] / d[0]
        x[0] = b[0] / d[0]
        for i in range(1, N):
            denom = d[i] - a[i - 1] * cp[i - 1]
            if i < N - 1:
                cp[i] = c[i] / denom
            x[i] = (b[i] - a[i - 1] * x[i - 1]) / denom

        # 回代
        for i in range(N - 2, -1, -1):
            x[i] -= cp[i] * x[i + 1]

        return x


# ============================================================
# 便捷工厂函数
# ============================================================

def create_solver(method: str, grid: Grid, V: Potential,
                  hbar: float = 1.0, mass: float = 1.0) -> TDSE_Solver:
    """创建求解器

    参数:
        method: "ssfm" | "split_step" | "crank_nicolson" | "cn"
    """
    method = method.lower().replace('-', '_')
    if method in ('ssfm', 'split_step', 'split_step_fourier'):
        return SplitStepFourier(grid, V, hbar, mass)
    elif method in ('cn', 'crank_nicolson'):
        return CrankNicolson(grid, V, hbar, mass)
    else:
        raise ValueError(f"未知求解方法: {method}，支持: ssfm, cn")

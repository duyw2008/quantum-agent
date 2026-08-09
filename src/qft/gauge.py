"""U(1) 规范场 — 光子动量空间模式展开 + Ward 恒等式

A_μ 在 3+1D 中的动量空间展开:
    A_μ(x) = Σ_{k,λ} [ε_μ(k,λ) a_{k,λ} e^{-ikx} + ε_μ*(k,λ) a†_{k,λ} e^{ikx}] / √(2ω_k V)

偏振矢量 ε_μ(k,λ) 满足:
    k^μ ε_μ = 0     (横向条件, 对 λ=1,2)
    ε^μ*(k,λ) ε_μ(k,λ') = -δ_{λλ'}  (归一化)

Minkowski 度规: g^{μν} = diag(1, -1, -1, -1)

Feynman 规范光子传播子:
    D_{μν}(k) = -i g_{μν} / (k² + iε)

Ward 恒等式: k_μ M^μ = 0 (规范不变性要求)

参考文献:
    Peskin & Schroeder, 第 9 章
    Weinberg Vol.1, 第 5, 8 章
"""

import numpy as np
from typing import Tuple, Callable, Optional


# ================================================================
# Minkowski 度规
# ================================================================

def minkowski_metric(signature: str = 'west_coast') -> np.ndarray:
    """Minkowski 度规 g^{μν} (逆变形式)

    参数:
        signature: 'west_coast' → (+, -, -, -)  [Peskin 约定]
                   'east_coast' → (-, +, +, +)  [Weinberg 约定]

    返回:
        (4, 4) 对角矩阵

    使用方式:
        a_dot_b = Σ a_μ g^{μν} b_ν = a @ g @ b
        a² = a @ g @ a
    """
    g = np.diag([1.0, -1.0, -1.0, -1.0]) if signature == 'west_coast' \
        else np.diag([-1.0, 1.0, 1.0, 1.0])
    return g


G_MUNU = minkowski_metric()  # g_{μν} (协变 = 逆变, 因为对角)


def lorentz_dot(k1: np.ndarray, k2: np.ndarray) -> complex:
    """四动量内积 k1·k2 = k1^μ g_{μν} k2^ν

    参数:
        k1, k2: 4 分量四矢量 (E, px, py, pz)

    返回:
        标量 k1·k2 (实数或复数)
    """
    result = k1[0] * k2[0] - k1[1] * k2[1] - k1[2] * k2[2] - k1[3] * k2[3]
    return complex(result)


def k_squared(k_mu: np.ndarray) -> float:
    """k² = k^μ k_μ = E² - p²

    参数:
        k_mu: 4 分量四动量

    返回:
        k² (质量壳条件: k² = m²)
    """
    return float(lorentz_dot(k_mu, k_mu).real)


# ================================================================
# 偏振矢量
# ================================================================

def polarization_vectors(k_4vec: np.ndarray, lambda_idx: int) -> np.ndarray:
    """光子的横向偏振矢量 ε_μ(k, λ)

    在 Coulomb 规范中, 取 k = (ω, 0, 0, |k|) 方向:
        ε^μ(k, 1) = (0, 1, 0, 0)   // x 方向线偏振
        ε^μ(k, 2) = (0, 0, 1, 0)   // y 方向线偏振

    满足: k_μ ε^μ = 0, ε_μ* ε^μ = -1

    参数:
        k_4vec:     光子四动量 (ω, kx, ky, kz)
        lambda_idx: 偏振指标 1 或 2 (1-indexed)

    返回:
        (4,) 复数数组 ε_μ
    """
    if lambda_idx not in (1, 2):
        raise ValueError(f"偏振 λ 必须是 1 或 2, 得到 {lambda_idx}")

    k = k_4vec[1:]  # 空间部分
    k_mag = float(np.linalg.norm(k))

    if k_mag < 1e-15:
        # k=0 极限: ε_μ(0) = (0, e_x) 或 (0, e_y)
        eps = np.zeros(4, dtype=complex)
        eps[1 if lambda_idx == 1 else 2] = 1.0
        return eps

    # 构建与 k 正交的偏振基矢
    k_hat = k / k_mag

    # 构建与 k_hat 正交的矢量 (避开退化情况)
    if abs(k_hat[0]) < abs(k_hat[2]):
        # k 主要在 z 方向: 用 x 轴作参考
        e1 = np.array([0.0, -k_hat[2], k_hat[1]])
    else:
        # k 主要在 x 方向: 用 z 轴作参考
        e1 = np.array([-k_hat[1], k_hat[0], 0.0])

    e1_norm = np.linalg.norm(e1)
    if e1_norm < 1e-15:
        # 完全退化: k = (0,0,±1) → e1 = (1,0,0)
        e1 = np.array([1.0, 0.0, 0.0])
    else:
        e1 = e1 / e1_norm
    e2 = np.cross(k_hat, e1)

    eps_spatial = e1 if lambda_idx == 1 else e2
    eps = np.zeros(4, dtype=complex)
    eps[0] = 0.0  # 时间分量为零 (Coulomb 规范)
    eps[1:] = eps_spatial

    return eps


def check_polarization(k_4vec: np.ndarray) -> dict:
    """验证偏振矢量的横向条件与归一化

    返回:
        {'k_dot_eps1': ..., 'k_dot_eps2': ..., 'eps1²': ..., 'eps2²': ..., 'eps1·eps2': ...}
    """
    eps1 = polarization_vectors(k_4vec, 1)
    eps2 = polarization_vectors(k_4vec, 2)

    k_dot_eps1 = lorentz_dot(k_4vec, eps1)
    k_dot_eps2 = lorentz_dot(k_4vec, eps2)
    eps1_sq = lorentz_dot(eps1, eps1)
    eps2_sq = lorentz_dot(eps2, eps2)
    eps1_dot_eps2 = lorentz_dot(eps1, eps2)

    return {
        'k·ε₁': abs(k_dot_eps1),
        'k·ε₂': abs(k_dot_eps2),
        'ε₁²': eps1_sq,
        'ε₂²': eps2_sq,
        'ε₁·ε₂': abs(eps1_dot_eps2),
    }


# ================================================================
# 光子传播子
# ================================================================

def photon_propagator(k_sq: float, xi: float = 1.0,
                      epsilon: float = 1e-12) -> np.ndarray:
    """光子传播子 D_{μν}(k) (Feynman 规范)

    一般 R_ξ 规范:
        D_{μν}(k) = (-i/k²) [g_{μν} - (1-ξ) k_μ k_ν / k²]

    Feynman 规范 (ξ=1):
        D_{μν}(k) = -i g_{μν} / k²

    Landau 规范 (ξ=0):
        D_{μν}(k) = (-i/k²) [g_{μν} - k_μ k_ν / k²]

    参数:
        k_sq:    k² (动量平方)
        xi:      规范参数 (ξ=1 → Feynman, ξ=0 → Landau)
        epsilon: 红外截断 (避免 k²=0 发散)

    返回:
        (4, 4) 复数矩阵 D_{μν}
    """
    k2 = k_sq if abs(k_sq) > epsilon else np.sign(k_sq) * epsilon if abs(k_sq) > 0 else epsilon
    factor = -1j / k2
    return factor * G_MUNU.astype(complex)


def photon_propagator_landau(k_4vec: np.ndarray,
                              epsilon: float = 1e-12) -> np.ndarray:
    """光子传播子在 Landau 规范 (ξ=0)

    D_{μν}(k) = (-i/k²) [g_{μν} - k_μ k_ν / k²]

    参数:
        k_4vec:  光子四动量
        epsilon: 红外截断

    返回:
        (4, 4) 复数矩阵
    """
    k2 = k_squared(k_4vec)
    if abs(k2) < epsilon:
        k2 = epsilon if k2 >= 0 else -epsilon

    g = G_MUNU.astype(complex)
    k_mu = k_4vec.astype(complex)
    tensor = g - np.outer(k_mu, k_mu) / k2
    return -1j / k2 * tensor


# ================================================================
# 场强张量 (占位)
# ================================================================

def field_strength_tensor(A_mu: callable, x_mu: np.ndarray,
                           dx: float = 1e-6) -> np.ndarray:
    """场强张量 F_{μν} = ∂_μ A_ν - ∂_ν A_μ  [占位实现]

    通过有限差分计算:

        F_{μν}(x) ≈ [A_ν(x + h ê_μ) - A_ν(x - h ê_μ) - (μ↔ν)] / 2h

    其中 ê_μ 是第 μ 方向的单位矢量。

    参数:
        A_mu:  函数 A_μ(x), 接受 (4,) 数组返回 (4,) 数组
        x_mu:  位置 (4,)
        dx:    有限差分步长

    返回:
        (4, 4) 反对称矩阵 F_{μν}
    """
    F = np.zeros((4, 4), dtype=complex)
    for mu in range(4):
        for nu in range(4):
            x_plus_mu = x_mu.copy()
            x_plus_mu[mu] += dx
            x_minus_mu = x_mu.copy()
            x_minus_mu[mu] -= dx
            dmu_A_nu = (A_mu(x_plus_mu)[nu] - A_mu(x_minus_mu)[nu]) / (2 * dx)

            x_plus_nu = x_mu.copy()
            x_plus_nu[nu] += dx
            x_minus_nu = x_mu.copy()
            x_minus_nu[nu] -= dx
            dnu_A_mu = (A_mu(x_plus_nu)[mu] - A_mu(x_minus_nu)[mu]) / (2 * dx)

            F[mu, nu] = dmu_A_nu - dnu_A_mu

    return F


# ================================================================
# 规范变换 (占位)
# ================================================================

def gauge_transform(A_mu: np.ndarray, alpha: complex) -> np.ndarray:
    """U(1) 规范变换 A_μ → A_μ - ∂_μ α  [占位实现]

    简化: 对均匀规范变换返回相同的 A_μ。
    完整实现需要知道 ∂_μ α。

    参数:
        A_mu:  4 分量矢势 (4,)
        alpha: 规范参数 α(x)

    返回:
        变换后的 A_μ (4,)
    """
    return A_mu.astype(complex).copy()


# ================================================================
# Ward 恒等式验证
# ================================================================

def ward_identity_check(amplitude_fn: Callable[[np.ndarray], np.ndarray],
                         k_mu: np.ndarray,
                         epsilon: float = 0.01,
                         tol: float = 1e-8) -> dict:
    """验证 Ward 恒等式 k_μ M^μ = 0

    对于在壳光子, 振幅 M^μ 满足:
        k_μ M^μ = 0  (Ward 恒等式, 规范不变性的直接推论)

    测试方法: 替换 ε_μ → k_μ 后振幅应当消失。

    参数:
        amplitude_fn: 振幅函数 f(k_μ) → M^μ  (4 分量数组)
        k_mu:         光子四动量 (基值)
        epsilon:      动量微扰用于有限差分
        tol:          零判断容差

    返回:
        {'k_dot_M': 值, 'passes': bool, 'norm': |M^μ|, 'direction_checks': [...]}
    """
    # 主检查: k_μ M^μ(k)
    M = amplitude_fn(k_mu)
    M = np.asarray(M, dtype=complex)

    k_dot_M = lorentz_dot(k_mu, M)
    result = {
        'k·M': abs(k_dot_M),
        'passes': abs(k_dot_M) < tol * max(abs(k_mu).max(), 1.0),
        'norm': float(np.linalg.norm(M)),
    }

    # 方向检查: 在多个方向的动量上验证
    direction_checks = []
    for i in range(1, 4):
        k_perturbed = k_mu.copy().astype(complex)
        k_perturbed[i] += epsilon
        M_p = amplitude_fn(k_perturbed)
        M_p = np.asarray(M_p, dtype=complex)
        k_dot_M_p = abs(lorentz_dot(k_perturbed, M_p))
        direction_checks.append(float(k_dot_M_p))

    result['direction_checks'] = direction_checks
    return result


# ================================================================
# GaugeField 类 — 光子动量空间模式展开
# ================================================================

class GaugeField:
    """3+1D 自由 U(1) 规范场 A_μ(x) 的简单表示

    模式展开:
        A_μ(x) = Σ_{k,λ} [ε_μ(k,λ) a_{k,λ} e^{-ikx} + ε_μ*(k,λ) a†_{k,λ} e^{ikx}] / √(2ω_k V)

    参数:
        volume:     正则化体积 V (默认 1.0)
        k_modes:     离散动量列表 [(ω, kx, ky, kz), ...]
        polarizations: 每模式的偏振列表 (默认 [1, 2])
    """

    def __init__(self, volume: float = 1.0,
                 k_modes: list = None,
                 polarizations: list = None):
        self.volume = volume
        self.polarizations = polarizations if polarizations else [1, 2]

        if k_modes is None:
            # 默认: 几个代表性动量模式
            self.k_modes = [
                np.array([1.0, 0.0, 0.0, 1.0]),   # 沿 z 方向, ω=1
                np.array([2.0, 1.0, 1.0, 0.0]),   # 沿 (1,1,0) 方向, ω=2
                np.array([0.5, 0.0, 0.5, 0.0]),   # 沿 y 方向, ω=0.5
            ]
        else:
            self.k_modes = [np.asarray(k) for k in k_modes]

        self.n_modes = len(self.k_modes)

        # 预计算 ω_k
        self.omega = np.array([k[0] for k in self.k_modes])

        # 归一化因子: 1/√(2ω_k V)
        self._norm = np.where(self.omega > 1e-15,
                               1.0 / np.sqrt(2 * self.omega * volume), 0.0)

    # ================================================================
    # 场算符
    # ================================================================

    def field_at(self, x_mu: np.ndarray,
                 amplitudes: np.ndarray = None) -> np.ndarray:
        """A_μ(x) — 在场点 x 处的矢势

        使用给定的产生/湮灭振幅 (经典场替代) 计算。

        参数:
            x_mu:       四维坐标 x^μ = (t, x, y, z)
            amplitudes: 模式振幅 (n_modes, 2) — 每个模式每个偏振的复数系数

        返回:
            4 分量复数数组 A_μ(x)
        """
        x_mu = np.asarray(x_mu, dtype=float)
        A = np.zeros(4, dtype=complex)

        if amplitudes is None:
            # 默认: 单位振幅
            amplitudes = np.ones((self.n_modes, 2), dtype=complex)

        for i, k in enumerate(self.k_modes):
            kx = lorentz_dot(k, x_mu)
            phase = np.exp(-1j * kx)
            norm = self._norm[i]

            for pol_idx in self.polarizations:
                eps = polarization_vectors(k, pol_idx)
                amp = amplitudes[i, pol_idx - 1]
                A += norm * (amp * phase * eps + np.conj(amp * phase) * eps.conj())

        return A

    def field_profile(self, x_points: np.ndarray) -> np.ndarray:
        """A_μ(x) 沿着 x 轴的轮廓 (固定 t)

        参数:
            x_points: (N,) 空间位置列表

        返回:
            (N, 4) 复数数组
        """
        N = len(x_points)
        profile = np.zeros((N, 4), dtype=complex)
        for i, x in enumerate(x_points):
            x_mu = np.array([0.0, x, 0.0, 0.0])
            profile[i] = self.field_at(x_mu)
        return profile

    # ================================================================
    # 传播子
    # ================================================================

    def propagator_representation(self, k_sq: float,
                                   xi: float = 1.0) -> np.ndarray:
        """光子传播子 D_{μν}(k) 在动量空间中的矩阵表示"""
        return photon_propagator(k_sq, xi)

    # ================================================================
    # 信息
    # ================================================================

    def summary(self) -> str:
        lines = [
            f"GaugeField(V={self.volume}, n_modes={self.n_modes})",
            f"  Polarizations: {self.polarizations}",
        ]
        for i, k in enumerate(self.k_modes):
            lines.append(
                f"  Mode {i}: k^μ=({k[0]:.1f}, {k[1]:.1f}, {k[2]:.1f}, {k[3]:.1f}), "
                f"ω={self.omega[i]:.3f}"
            )
        return '\n'.join(lines)

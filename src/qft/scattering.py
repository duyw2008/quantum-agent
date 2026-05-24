"""散射与 Feynman 图 — φ⁴ 理论微扰展开

核心概念:
    1. Wick 定理 — 场算符时序积展开为所有可能的收缩
    2. Dyson 级数 — S = T exp(-i∫V(t)dt)
    3. Feynman 振幅 — 树图 2→2 散射
    4. 跃迁概率 — 数值计算 |⟨f|U(t)|i⟩|²

全部在截断 Fock 空间中计算, 与 LatticePhi4 一致。
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import permutations, combinations


# ================================================================
# Wick 收缩
# ================================================================

def wick_expand(operators: List[str]) -> List[List[Tuple[int, int]]]:
    """Wick 定理: 将场算符的时序积展开为所有可能的收缩对

    参数:
        operators: 算符列表, 如 ['φ₁', 'φ₂', 'φ₃', 'φ₄']

    返回:
        所有可能的收缩方式, 每种方式是一个 (i,j) 对列表

    示例:
        wick_expand(['φ₁','φ₂','φ₃','φ₄']) →
        [
            [(0,1),(2,3)],  # φ₁φ₂ 收缩, φ₃φ₄ 收缩
            [(0,2),(1,3)],  # φ₁φ₃ 收缩, φ₂φ₄ 收缩
            [(0,3),(1,2)],  # φ₁φ₄ 收缩, φ₂φ₃ 收缩
        ]
    """
    n = len(operators)
    if n % 2 != 0:
        return []  # 奇数个算符无完全收缩

    if n == 0:
        return [[]]
    if n == 2:
        return [[(0, 1)]]

    results = []
    # 取第一个算符, 与剩余每个配对
    for j in range(1, n):
        # 递归处理剩余 n-2 个算符
        remaining = [k for k in range(n) if k != 0 and k != j]
        # 重新映射索引
        sub_results = wick_expand(['_'] * (n - 2))
        if not sub_results:
            results.append([(0, j)])
        else:
            for sub in sub_results:
                remapped = []
                for a, b in sub:
                    ra = remaining[a]
                    rb = remaining[b]
                    remapped.append((ra, rb))
                results.append([(0, j)] + remapped)

    return results


def propagator(mass: float, dx: float = 0.0, dt: float = 0.0) -> complex:
    """自由标量场传播子 (动量空间, 1+1D)

    D_F(p) = i / (p² - m² + iε)

    位置空间 (大质量极限近似):
        D_F(x, t) ≈ (1/2m) e^{-m|x|} δ(t)  (瞬时近似)
    """
    if abs(dt) < 1e-15:
        return complex(1.0 / (2 * mass) * np.exp(-mass * abs(dx)))
    # 含时传播子 (简化)
    return complex(np.exp(-mass * abs(dx) - 1j * mass * abs(dt)) / (2 * mass))


# ================================================================
# Feynman 振幅 — φ⁴ 2→2 散射
# ================================================================

def feynman_amplitude_phi4_2to2(coupling: float, momenta: np.ndarray = None) -> complex:
    """φ⁴ 理论树图 2→2 散射振幅

    φ⁴ 相互作用: L_int = -(λ/4!) φ⁴

    树图振幅 (所有阶):
        iM = -iλ    (动量无关, 因为顶点无导数耦合)

    参数:
        coupling: λ
        momenta:  忽略 (φ⁴ 树图振幅与动量无关)

    返回:
        iM (复振幅)
    """
    return complex(0, -coupling)  # iM = -iλ


def differential_cross_section(coupling: float, s: float,
                                mass: float = 1.0) -> float:
    """微分散射截面 dσ/dΩ (质心系, 2→2)

    对于 φ⁴ 理论:
        dσ/dΩ = λ² / (64π² s)    (s ≫ m² 极限)

    参数:
        coupling: λ
        s:        Mandelstam s = E_cm²
        mass:     粒子质量

    返回:
        dσ/dΩ
    """
    amp2 = coupling**2  # |M|² = λ²
    # 相空间因子 (1+1D)
    if s <= 4 * mass**2:
        return 0.0
    p_cm = np.sqrt(s / 4 - mass**2)
    return float(amp2 / (16 * np.pi * s) * (p_cm / np.sqrt(s)))


# ================================================================
# Dyson 级数 — 时间演化算符的微扰展开
# ================================================================

def dyson_term(H0: np.ndarray, V: np.ndarray, t: float,
               order: int, hbar: float = 1.0) -> np.ndarray:
    """Dyson 级数第 n 阶项

    U^{(n)}(t) = (-i/ℏ)^n ∫₀^t dt₁ ... ∫₀^{t_{n-1}} dt_n
                  e^{-iH₀(t-t₁)/ℏ} V e^{-iH₀(t₁-t₂)/ℏ} V ... V e^{-iH₀t_n/ℏ}

    简化: 等时间步长近似
        U^{(n)}(t) ≈ (-iΔt/ℏ)^n Σ_{k₁<...<k_n}
                      e^{-iH₀(t-t_{k₁})/ℏ} V ... V e^{-iH₀t_{k_n}/ℏ}

    参数:
        H0:    自由哈密顿量
        V:     相互作用
        t:     总时间
        order: 阶数 n

    返回:
        U^{(n)} 矩阵
    """
    N_steps = max(20, order * 5)
    dt = t / N_steps

    # 自由演化算符
    eigvals, eigvecs = np.linalg.eigh(H0)

    def U0(tau):
        return eigvecs @ np.diag(np.exp(-1j * eigvals * tau / hbar)) @ eigvecs.conj().T

    # 相互作用绘景中的 V(t) = e^{iH₀t/ℏ} V e^{-iH₀t/ℏ}
    U_n = np.zeros_like(H0, dtype=complex)

    # 遍历所有可能的时间排序
    times = np.arange(1, N_steps + 1) * dt
    for indices in combinations(range(N_steps), order):
        term = np.eye(H0.shape[0], dtype=complex)
        t_prev = 0.0
        for idx in indices:
            ti = times[idx]
            term = U0(t - ti) @ V @ term  # simplified
            t_prev = ti

        factor = (-1j * dt / hbar) ** order
        U_n += factor * term

    return U_n


def transition_probability(H0: np.ndarray, V: np.ndarray,
                            initial: np.ndarray, final: np.ndarray,
                            t: float, max_order: int = 2,
                            hbar: float = 1.0) -> Dict[int, float]:
    """跃迁概率 |⟨f|U(t)|i⟩|² 逐阶计算

    参数:
        H0:        自由哈密顿量
        V:         相互作用
        initial:   初态 |i⟩
        final:     末态 |f⟩
        t:         演化时间
        max_order: 最大微扰阶数

    返回:
        {order: probability} 字典
    """
    # 对角化 H0
    eigvals, eigvecs = np.linalg.eigh(H0)

    def U0(tau):
        return eigvecs @ np.diag(np.exp(-1j * eigvals * tau / hbar)) @ eigvecs.conj().T

    # 相互作用绘景 V_I(t) = U0†(t) V U0(t)
    result = {}
    U_total = np.eye(H0.shape[0], dtype=complex)

    for n in range(max_order + 1):
        if n == 0:
            amp = final.conj() @ U0(t) @ initial
            result[0] = float(abs(amp)**2)
            continue

        # 第 n 阶: 时间排序积分
        N_int = 50
        dt_val = t / N_int
        U_n = np.zeros_like(H0, dtype=complex)

        # 用离散时间求和近似连续积分
        for indices in combinations(range(N_int), n):
            # 时间排序
            sorted_idx = sorted(indices)
            t_vals = [(idx + 0.5) * dt_val for idx in sorted_idx]

            # 构建每一阶的贡献
            contrib = np.eye(H0.shape[0], dtype=complex)
            for ti in reversed(t_vals):
                V_I = U0(ti).conj().T @ V @ U0(ti)
                contrib = (-1j / hbar) * V_I @ contrib

            U_n += contrib * dt_val**n

        U_n = U0(t) @ U_n
        amp = final.conj() @ U_n @ initial
        result[n] = float(abs(amp)**2)

    return result


# ================================================================
# Feynman 图可视化 (ASCII)
# ================================================================

def draw_feynman_phi4_2to2() -> str:
    """绘制 φ⁴ 2→2 树图 Feynman 图 (ASCII)"""
    return r"""
    p₁    p₃
     \   /
      \ /
       ●   = -iλ
      / \
     /   \
    p₂    p₄

    iM = -iλ
    |M|² = λ²
    dσ/dΩ = λ²/(64π²s)   (s ≫ m²)
    """


def draw_feynman_phi4_loop() -> str:
    """绘制 φ⁴ 单圈自能图"""
    return r"""
       ____
      /    \
     /  ⭕  \    单圈修正:
    ●        ●   Π(p²) ∝ λ ∫ d²k/(k²-m²)
     \      /
      \____/
    """

"""有效势与自发对称性破缺 — 单圈 Coleman-Weinberg 有效势

核心概念:
    1. 有效势 V_eff(φ_c) — 经典势 + 单圈量子修正
    2. Coleman-Weinberg 机制 — 无质量标量 QED 的辐射对称性破缺
    3. 对称性破缺 — 序参量、Goldstone 定理、Higgs 机制
    4. 数值最小值搜索 — find_minimum 在给定区间寻找 V_eff 极小值
    5. 可视化数据 — potential_plot_data 输出 (φ, V_eff, V_tree) 供绘图使用

全部基于 numpy, 无外部依赖。
"""

import numpy as np
from typing import Tuple, Callable, Optional


# ================================================================
# 单圈有效势
# ================================================================

class OneLoopEffectivePotential:
    """φ⁴ 理论单圈有效势

    有效势 = 树图势 + 单圈量子修正 (Coleman-Weinberg 势):

        V_eff(φ_c) = V_tree(φ_c) + V_1loop(φ_c)

        V_tree(φ_c) = ½m²φ_c² + (λ/4!)φ_c⁴

        V_1loop(φ_c) = (ħ/64π²) M⁴(φ_c) [ln(M²(φ_c)/μ²) - 3/2]

    其中 M²(φ_c) = m² + ½λφ_c² 是场依赖的有效质量平方。

    M²(φ_c) < 0 时, V_1loop 的实部来自 |M²| (取模),
    虚部对应不稳定模式的衰减率 Γ = Im(V_eff)/ħ。

    参数:
        mass_sq:  裸质量平方 m² (负值 → 树图自发对称性破缺, 零 → Coleman-Weinberg)
        coupling: 耦合常数 λ
        mu:       重整化标度 μ
        hbar:     ℏ (默认 1.0)
    """

    def __init__(self, mass_sq: float = 0.0, coupling: float = 0.1,
                 mu: float = 1.0, hbar: float = 1.0):
        self.mass_sq = mass_sq
        self.coupling = coupling
        self.mu = mu
        self.hbar = hbar

        # 验证: λ > 0 (真空稳定性)
        if coupling <= 0:
            raise ValueError(f"耦合常数 λ 必须为正 (真空稳定性), 得到 λ={coupling}")

    # ================================================================
    # 场依赖的有效质量 M²(φ_c)
    # ================================================================

    def effective_mass_sq(self, phi_c: float) -> float:
        """场依赖的有效质量平方

        M²(φ_c) = m² + ½λφ_c² = d²V_tree/dφ_c²

        这是单圈涨落算符的本征值。

        参数:
            phi_c: 经典背景场 φ_c

        返回:
            M²(φ_c)
        """
        return self.mass_sq + 0.5 * self.coupling * phi_c**2

    # ================================================================
    # 树图势
    # ================================================================

    def tree_potential(self, phi_c: float) -> float:
        """树图经典势

        V_tree(φ_c) = ½m²φ_c² + (λ/4!)φ_c⁴

        当 m² < 0 时, 势有两个非零极小值:
            φ_min = ±√(-6m²/λ)

        参数:
            phi_c: 经典背景场

        返回:
            V_tree(φ_c)
        """
        return 0.5 * self.mass_sq * phi_c**2 + (self.coupling / 24.0) * phi_c**4

    # ================================================================
    # 单圈量子修正
    # ================================================================

    def one_loop_correction(self, phi_c: float) -> float:
        """单圈 Coleman-Weinberg 量子修正

        V_1loop(φ_c) = (ħ/64π²) M⁴(φ_c) [ln(M²(φ_c)/μ²) - 3/2]

        当 M² < 0 时:
            ln(M²/μ²) = ln(|M²|/μ²) + iπ   (主值分支)
            Re[V_1loop] = (ħ/64π²) (M²)² [ln(|M²|/μ²) - 3/2]
            Im[V_1loop] = -ħ|M²|²/64π    ← 衰变率 Γ = |Im(V_eff)|/ħ

        当 M² → 0 (CW 在 φ_c → 0):
            M⁴ ln(M²) → 0  (极限为 0)

        参数:
            phi_c: 经典背景场

        返回:
            Re[V_1loop(φ_c)] (实数部分)
        """
        M2 = self.effective_mass_sq(phi_c)
        abs_M2 = abs(M2)

        # M² → 0 极限: M⁴ ln(M²) → 0
        if abs_M2 < 1e-30:
            return 0.0

        if M2 <= 0:
            # 不稳定区域, 取实部 (使用 |M²|)
            return float(self.hbar / (64.0 * np.pi**2) * abs_M2**2
                         * (np.log(abs_M2 / self.mu**2) - 1.5))
        else:
            return float(self.hbar / (64.0 * np.pi**2) * M2**2
                         * (np.log(M2 / self.mu**2) - 1.5))

    def one_loop_correction_complex(self, phi_c: float) -> complex:
        """单圈量子修正 — 完整复数形式 (包含虚部)

        虚部来自 M² < 0 时 ln 的支割线:
            Im[V_1loop] = -ħ |M²|² / (64π)

        参数:
            phi_c: 经典背景场

        返回:
            V_1loop(φ_c) (复数)
        """
        M2 = self.effective_mass_sq(phi_c)
        abs_M2 = abs(M2)
        if abs_M2 < 1e-30:
            return complex(0.0, 0.0)
        if M2 > 0:
            return complex(self.hbar / (64.0 * np.pi**2) * M2**2
                           * (np.log(M2 / self.mu**2) - 1.5))
        else:
            real_part = self.hbar / (64.0 * np.pi**2) * abs_M2**2 * (np.log(abs_M2 / self.mu**2) - 1.5)
            imag_part = -self.hbar * abs_M2**2 / (64.0 * np.pi)
            return complex(real_part, imag_part)

    # ================================================================
    # 完整有效势
    # ================================================================

    def effective_potential(self, phi_c: float) -> float:
        """完整单圈有效势 (实数部分)

        V_eff(φ_c) = V_tree(φ_c) + Re[V_1loop(φ_c)]

        参数:
            phi_c: 经典背景场

        返回:
            V_eff(φ_c)

        >>> pot = OneLoopEffectivePotential(mass=1.0, coupling=0.1, mu=1.0)
        >>> pot.effective_potential(0.0)
        0.0
        >>> pot.effective_potential(1.0)
        0.504...  # ½m² + λ/24 + 量子修正
        """
        return self.tree_potential(phi_c) + self.one_loop_correction(phi_c)

    def effective_potential_array(self, phi_vals: np.ndarray) -> np.ndarray:
        """有效势矢量化计算

        参数:
            phi_vals: 背景场值数组

        返回:
            V_eff(phi_vals) 数组
        """
        # 矢量化计算
        phi_sq = phi_vals**2
        phi_4 = phi_vals**4
        tree = 0.5 * self.mass_sq * phi_sq + (self.coupling / 24.0) * phi_4

        M2 = self.mass_sq + 0.5 * self.coupling * phi_sq
        abs_M2 = np.abs(M2)

        # 单圈修正 (矢量化, 处理 M² ≤ 0, 以及 M² → 0 极限)
        loop = np.zeros_like(phi_vals, dtype=float)
        tiny_mask = abs_M2 < 1e-30
        pos_mask = (M2 > 0) & ~tiny_mask
        neg_mask = (M2 <= 0) & ~tiny_mask

        if np.any(pos_mask):
            loop[pos_mask] = (self.hbar / (64.0 * np.pi**2) * M2[pos_mask]**2
                              * (np.log(M2[pos_mask] / self.mu**2) - 1.5))
        if np.any(neg_mask):
            loop[neg_mask] = (self.hbar / (64.0 * np.pi**2) * abs_M2[neg_mask]**2
                              * (np.log(abs_M2[neg_mask] / self.mu**2) - 1.5))

        return tree + loop

    # ================================================================
    # 信息
    # ================================================================

    def summary(self) -> str:
        """有效势信息摘要"""
        m2 = self.mass_sq
        lam = self.coupling
        lines = [
            f"OneLoopEffectivePotential(m²={m2}, λ={lam}, μ={self.mu}, ħ={self.hbar})",
            f"  V_tree = ½m²φ² + λφ⁴/24",
            f"  V_1loop = ħM⁴[ln(M²/μ²)-3/2]/(64π²),  M²=m²+½λφ²",
        ]
        if m2 < 0:
            v_tree = np.sqrt(-6.0 * m2 / lam)
            lines.append(f"  树图 SSB: φ_min = ±{v_tree:.4f}")
        elif abs(m2) < 1e-15:
            lines.append("  Coleman-Weinberg: m=0, 辐射 SSB")
        else:
            lines.append(f"  对称相: φ_min = 0 (m² > 0)")
        return '\n'.join(lines)


# ================================================================
# Coleman-Weinberg 势 (m=0 特殊情况)
# ================================================================

def coleman_weinberg_potential(phi_c: float, lam: float, mu: float,
                               hbar: float = 1.0) -> float:
    """Coleman-Weinberg 有效势 — 无质量标量 QED / λφ⁴ (m=0)

    在 m=0 时有效势简化为:

        V_CW(φ_c) = (λ/4!)φ_c⁴ + (ħλ²φ_c⁴/(256π²)) [ln(λφ_c²/(2μ²)) - 3/2]

    由于对数项在小的 φ_c 处为负, ϕ=0 不再是极小值点 →
    辐射修正诱导自发对称性破缺。

    极小值条件 dV/dφ_c = 0:
        φ_min = (2μ²/λ)^{1/2} exp(1/2 - 8π²/(3ħλ))

    参数:
        phi_c: 经典背景场
        lam:   耦合常数 λ
        mu:    重整化标度 μ
        hbar:  ℏ (默认 1.0)

    返回:
        V_CW(φ_c)

    >>> coleman_weinberg_potential(1.0, 0.1, 1.0)
    0.004...  # λ/24 + 量子修正
    """
    if abs(phi_c) < 1e-15:
        return 0.0

    phi2 = phi_c * phi_c
    phi4 = phi2 * phi2

    tree = (lam / 24.0) * phi4
    M2 = 0.5 * lam * phi2  # M²(φ_c) = ½λφ_c² (m=0)
    M4 = M2 * M2           # λ²φ_c⁴/4

    log_arg = max(M2 / mu**2, 1e-30)
    loop = hbar / (64.0 * np.pi**2) * M4 * (np.log(log_arg) - 1.5)
    return float(tree + loop)


def coleman_weinberg_minimum(lam: float, mu: float,
                             hbar: float = 1.0) -> float:
    """Coleman-Weinberg 势极小值的解析位置 (近似)

    dV_CW/dφ_c = 0 的解析解:

        φ_min² = (2μ²/λ) exp[1 - 32π²/(3ħλ)]
               = (2μ²/λ) exp[1] · exp[-32π²/(3ħλ)]

    当 ħλ/(32π²) ≪ 1 时, φ_min ≪ 1 → 微扰论可靠。

    参数:
        lam:  耦合 λ
        mu:   重整化标度 μ
        hbar: ℏ

    返回:
        φ_min > 0 (正分支)

    >>> coleman_weinberg_minimum(0.1, 1.0)
    7.8e-136...  # 极小, 因为 exp(-32π²/(3*0.1))
    """
    exponent = 1.0 - 32.0 * np.pi**2 / (3.0 * hbar * lam)
    phi_min_sq = (2.0 * mu**2 / lam) * np.exp(exponent)
    return float(np.sqrt(max(phi_min_sq, 0.0)))


# ================================================================
# 最小值搜索
# ================================================================

def find_minimum(veff_func: Callable[[float], float],
                 phi_range: Tuple[float, float],
                 n_points: int = 1000) -> Tuple[float, float]:
    """在给定区间内寻找有效势的全局极小值

    使用网格扫描 + 抛物插值精确定位。

    算法:
        1. 在 phi_range 上均匀采样 n_points 个点
        2. 找到 V_eff 最小的点作为粗估计
        3. 在该点附近用三点抛物插值精确定位极小值

    参数:
        veff_func:  V_eff(phi_c) 函数
        phi_range:  (phi_min, phi_max) 搜索区间
        n_points:   网格点数

    返回:
        (phi_min, V_min) — 极小值位置与势值

    >>> pot = OneLoopEffectivePotential(mass=-1.0, coupling=0.1)
    >>> phi_min, V_min = find_minimum(pot.effective_potential, (-20, 20))
    >>> abs(phi_min - np.sqrt(6/0.1)) < 0.1  # ≈ 7.746
    True
    """
    phi_lo, phi_hi = phi_range
    phi_vals = np.linspace(phi_lo, phi_hi, n_points)
    v_vals = np.array([veff_func(p) for p in phi_vals])

    # 找到全局最小值的索引
    idx_min = int(np.argmin(v_vals))

    # 抛物插值精确定位 (用相邻三点)
    if 0 < idx_min < n_points - 1:
        x = phi_vals[idx_min - 1:idx_min + 2]
        y = v_vals[idx_min - 1:idx_min + 2]
        # 拟合抛物 y = ax² + bx + c → 极小值在 x* = -b/(2a)
        # 使用二阶差分
        h = x[1] - x[0]
        a = (y[2] - 2 * y[1] + y[0]) / (2 * h**2)
        if abs(a) > 1e-15:
            b = (y[2] - y[0]) / (2 * h)
            phi_fine = x[1] - b / (2 * a)
            # 确保插值结果在区间内
            phi_fine = np.clip(phi_fine, x[0], x[2])
            v_fine = float(veff_func(phi_fine))
            if v_fine < v_vals[idx_min]:
                return phi_fine, v_fine

    return float(phi_vals[idx_min]), float(v_vals[idx_min])


def find_all_minima(veff_func: Callable[[float], float],
                    phi_range: Tuple[float, float],
                    n_points: int = 2000) -> list:
    """寻找有效势的所有局域极小值

    使用导数符号变化检测局域极小值。

    参数:
        veff_func:  V_eff(φ_c) 函数
        phi_range:  搜索区间
        n_points:   网格点数

    返回:
        [(φ_min, V_min), ...] 按 V 值升序排列
    """
    phi_lo, phi_hi = phi_range
    phi_vals = np.linspace(phi_lo, phi_hi, n_points)
    v_vals = np.array([veff_func(p) for p in phi_vals])

    # 数值导数
    dv = np.diff(v_vals) / np.diff(phi_vals)

    minima = []
    for i in range(1, len(dv)):
        if dv[i - 1] < 0 and dv[i] > 0:
            # 导数从负变正 → 局域极小值
            phi_min = phi_vals[i]
            v_min = v_vals[i]
            minima.append((float(phi_min), float(v_min)))
        elif abs(dv[i - 1]) < 1e-15 and abs(dv[i]) < 1e-15:
            # 平坦区域
            pass

    # 按 V 值升序排列
    minima.sort(key=lambda x: x[1])
    return minima


# ================================================================
# 自发对称性破缺
# ================================================================

class SymmetryBreaking:
    """自发对称性破缺分析

    对于 φ⁴ 理论, 当 m² < 0 或通过 Coleman-Weinberg 机制时,
    系统展现出 Z₂ → 无 或 U(1) → 无 的对称性破缺。

    关键可观测量:
        - 序参量 (order parameter): VEV ⟨φ⟩
        - Goldstone 模: 连续对称性破缺的无质量激发
        - Higgs 机制: 规范场获得质量

    参数:
        vev:         真空期望值 ⟨φ⟩ = v
        mass:        标量场质量 (破缺相)
        coupling:    自耦合 λ
        gauge_coupling: 规范耦合 g (Higgs 机制, 默认 None)
    """

    def __init__(self, vev: float, mass: float = 1.0,
                 coupling: float = 0.1,
                 gauge_coupling: Optional[float] = None):
        self.vev = vev
        self.mass = mass
        self.coupling = coupling
        self.gauge_coupling = gauge_coupling

    # ================================================================
    # 序参量
    # ================================================================

    @property
    def order_parameter(self) -> float:
        """序参量 — 真空期望值 ⟨φ⟩ = v

        在对称相中 v = 0, 在破缺相中 v ≠ 0。
        """
        return self.vev

    @property
    def is_broken(self) -> bool:
        """是否处于对称性破缺相"""
        return abs(self.vev) > 1e-15

    # ================================================================
    # 质量谱
    # ================================================================

    def scalar_mass(self) -> float:
        """破缺相的标量场质量 (Higgs 质量)

        在破缺相, 对 V(φ=v+η) 展开:
            V = const + ½ m_h² η² + ...

            m_h² = d²V/dφ²|_{φ=v} = m² + ½λv²

        对于树图 SSB (m² < 0):
            v² = -6m²/λ  →  m_h² = m² + ½λ(-6m²/λ) = -2m² > 0

        返回:
            m_h (正质量)
        """
        m2_eff = self.mass**2 + 0.5 * self.coupling * self.vev**2
        return float(np.sqrt(max(m2_eff, 0.0)))

    def goldstone_count(self) -> int:
        """Goldstone 玻色子数目

        对于 Z₂ 破缺 (离散对称性): 0 个 Goldstone
        对于 U(1) 破缺 (连续对称性): 1 个 Goldstone (当 gauge_coupling=0)

        返回:
            Goldstone 模数
        """
        # Z₂ → 无: 离散对称性破缺, 无 Goldstone
        if self.gauge_coupling is None:
            return 0  # 全局 Z₂ 破缺
        else:
            # U(1) 规范对称性: Goldstone 玻色子被吃掉
            return 0  # Higgs 机制: Goldstone → 规范玻色子纵向分量

    # ================================================================
    # Higgs 机制 (标量 QED / Abel Higgs)
    # ================================================================

    def higgs_mechanism(self) -> dict:
        """Higgs 机制 — Abel Higgs 模型中的质量生成

        复标量场 Φ = (1/√2)(v + h + iG) 与 U(1) 规范场 A_μ 耦合:

            L = |D_μ Φ|² - V(Φ) - ¼F_μν F^{μν}

            D_μ = ∂_μ + i g A_μ

        在幺正规范下, Goldstone G 被规范玻色子吸收:

            |D_μ Φ|² → ½(∂h)² + ½g²v² A_μ A^μ + g²v h A_μ A^μ + ...

        规范玻色子质量: m_A = g v
        Higgs 质量:     m_h = √(λ/3) v   (从 m² = -λv²/6)

        返回:
            {
                'gauge_boson_mass': m_A,
                'higgs_mass': m_h,
                'vev': v,
                'gauge_coupling': g,
                'goldstone_eaten': True,
            }
        """
        g = self.gauge_coupling if self.gauge_coupling else 0.0
        m_A = g * abs(self.vev)
        m_h = self.scalar_mass()

        return {
            'gauge_boson_mass': float(m_A),
            'higgs_mass': float(m_h),
            'vev': float(self.vev),
            'gauge_coupling': float(g),
            'goldstone_eaten': g > 1e-15,
            'mass_ratio': float(m_h / m_A) if m_A > 1e-15 else float('inf'),
        }

    def is_goldstone_mode(self, mass: float, tolerance: float = 1e-10) -> bool:
        """检查给定质量是否为 Goldstone 模 (近似零质量)

        在连续对称性自发破缺时, 存在零质量激发。

        参数:
            mass:      待检查的质量
            tolerance: 零质量容忍度

        返回:
            是否为零质量 Goldstone 模
        """
        return abs(mass) < tolerance

    # ================================================================
    # 有效势曲率 (质量矩阵)
    # ================================================================

    def curvature_at_vev(self, pot: OneLoopEffectivePotential) -> float:
        """在 VEV 处的有效势曲率 (质量平方)

        m_eff² = d²V_eff/dφ_c²|_{φ_c=v}  (数值差分)

        参数:
            pot: OneLoopEffectivePotential 实例

        返回:
            m_eff²
        """
        h = max(abs(self.vev) * 0.001, 1e-6)
        vp = self.vev + h
        vm = self.vev - h
        d2v = (pot.effective_potential(vp) - 2 * pot.effective_potential(self.vev)
               + pot.effective_potential(vm)) / h**2
        return float(d2v)

    # ================================================================
    # 信息
    # ================================================================

    def summary(self) -> str:
        lines = [
            f"SymmetryBreaking(v={self.vev:.4f}, m={self.mass}, λ={self.coupling})",
            f"  破缺相: {'是' if self.is_broken else '否'}",
            f"  序参量 ⟨φ⟩ = {self.vev:.4f}",
            f"  Higgs 质量 m_h = {self.scalar_mass():.4f}",
            f"  Goldstone 数 = {self.goldstone_count()}",
        ]
        if self.gauge_coupling:
            hm = self.higgs_mechanism()
            lines.append(f"  规范玻色子质量 m_A = {hm['gauge_boson_mass']:.4f}")
            lines.append(f"  质量比 m_h/m_A = {hm['mass_ratio']:.4f}")
        return '\n'.join(lines)


# ================================================================
# 可视化数据生成
# ================================================================

def potential_plot_data(pot: OneLoopEffectivePotential,
                        phi_range: Tuple[float, float] = (-5.0, 5.0),
                        n_points: int = 500) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """生成有效势绘图数据 — (φ, V_eff, V_tree)

    参数:
        pot:       OneLoopEffectivePotential 实例
        phi_range: φ 范围 (phi_min, phi_max)
        n_points:  采样点数

    返回:
        (phi_vals, V_eff_vals, V_tree_vals)

    示例:
        >>> pot = OneLoopEffectivePotential(mass=-1.0, coupling=0.1)
        >>> phi, Veff, Vtree = potential_plot_data(pot, (-10, 10))
        >>> len(phi)
        500
        >>> Veff.shape == phi.shape
        True
    """
    phi_lo, phi_hi = phi_range
    phi_vals = np.linspace(phi_lo, phi_hi, n_points)
    v_eff = pot.effective_potential_array(phi_vals)
    v_tree = np.array([pot.tree_potential(p) for p in phi_vals])
    return phi_vals, v_eff, v_tree


def symmetry_breaking_demo_data(pot: OneLoopEffectivePotential,
                                phi_range: Tuple[float, float] = (-5.0, 5.0),
                                n_points: int = 500) -> dict:
    """生成对称性破缺演示的完整数据集

    返回包含极值点标注的绘图数据。

    参数:
        pot:       OneLoopEffectivePotential 实例
        phi_range: φ 范围
        n_points:  采样点数

    返回:
        {
            'phi':        phi_vals,
            'V_eff':      V_eff_vals,
            'V_tree':     V_tree_vals,
            'minima':     [(phi, V), ...],
            'maximum':    (phi, V) at φ=0 (if SSB),
            'mass_sq':    m²,
            'coupling':   λ,
        }
    """
    phi, Veff, Vtree = potential_plot_data(pot, phi_range, n_points)

    minima = find_all_minima(pot.effective_potential, phi_range, n_points)

    # φ=0 处的势值
    v_at_zero = pot.effective_potential(0.0)

    # 判断 φ=0 是否为局域极大值 (SSB 的特征)
    is_max = False
    if len(minima) > 0 and all(abs(m[0]) > 0.1 for m in minima):
        v_min_val = minima[0][1] if minima else float('inf')
        is_max = v_at_zero > v_min_val + 1e-10

    return {
        'phi': phi,
        'V_eff': Veff,
        'V_tree': Vtree,
        'minima': minima,
        'maximum': (0.0, v_at_zero) if is_max else None,
        'mass_sq': pot.mass_sq,
        'coupling': pot.coupling,
        'mu': pot.mu,
    }

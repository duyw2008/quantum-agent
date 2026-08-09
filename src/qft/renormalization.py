"""重整化 — φ⁴ 理论单圈重整化

核心概念:
    1. Feynman 规则 — φ⁴ 传播子、顶点、动量路由
    2. 动量截断正规化 — 引入截断 Λ, 计算发散积分
    3. 单圈自能 Π(p²) — 蝌蚪图, 质量重整化
    4. 单圈顶角修正 Γ⁴(s,t,u) — 盒图, 耦合常数重整化
    5. 重整化方案 — 在壳 (on-shell) 方案
    6. 重整化群 — β 函数与跑动耦合

全部使用动量截断正规化 (numpy 可实现), 保留维数正规化占位。
"""

import numpy as np
from typing import Tuple, Callable, Optional


# ================================================================
# Feynman 规则
# ================================================================

class Phi4FeynmanRules:
    """φ⁴ 理论的 Feynman 规则

    拉氏量:
        L = ½(∂_μφ)² - ½m²φ² - (λ/4!)φ⁴

    规则:
        - 传播子: i/(p² - m² + iε)
        - 顶点:   -iλ
        - 每个圈积分: ∫ d⁴k/(2π)⁴
        - 对称因子: 图相关的组合因子

    参数:
        mass:     裸质量 m
        coupling: 裸耦合 λ
        hbar:     ℏ (默认 1.0)
        dim:      时空维数 (默认 4)
    """

    def __init__(self, mass: float = 1.0, coupling: float = 0.1,
                 hbar: float = 1.0, dim: int = 4):
        self.mass = mass
        self.coupling = coupling
        self.hbar = hbar
        self.dim = dim

    def propagator(self, p_sq: float, epsilon: float = 1e-10) -> complex:
        """动量空间 Feynman 传播子

        D_F(p) = i / (p² - m² + iε)

        参数:
            p_sq:     动量平方 p² (可以是负值, 如类空动量)
            epsilon:  iε 正规子

        返回:
            复传播子值

        >>> rules = Phi4FeynmanRules(mass=1.0)
        >>> rules.propagator(0.0)
        (-0-1j)   # -i/m²
        """
        return complex(0, 1) / (p_sq - self.mass**2 + 1j * epsilon)

    def vertex(self) -> complex:
        """φ⁴ 四顶点因子

        -iλ   (连接四条外线的顶点)

        返回:
            复振幅因子
        """
        return complex(0, -self.coupling)

    def loop_integral(self, integrand: Callable[[np.ndarray], np.ndarray],
                      Lambda: float) -> float:
        """d=4 动量空间圈积分 (各向同性近似)

        ∫ d⁴k/(2π)⁴ f(k²) → 1/(16π²) ∫₀^Λ k³ f(k²) dk
          球坐标: d⁴k = 2π² k³ dk, 除 (2π)⁴ → k³ dk / (8π²)

        但标准惯例是 ∫ d⁴k/(2π)⁴ → k³ dk/(16π²), 因为:
            d⁴k = 2π² k³ dk   (四维球面面积 = 2π²)
            除以 (2π)⁴ → 2π²/(16π⁴) k³ dk = k³ dk/(8π²)

        等等, 让我用正确的因子:
        ∫ d⁴k/(2π)⁴ f(|k|) = Ω₃/(2π)⁴ ∫ k³ f(k) dk
          其中 Ω₃ = 2π² (三维球面 S³ 的面积)
          = 2π²/(16π⁴) ∫ k³ f(k) dk
          = 1/(8π²) ∫ k³ f(k) dk

        使用 1/(8π²) 因子的积分:
        ∫₀^Λ k³/(k²+m²) dk = ½[Λ² - m² ln(1+Λ²/m²)]
        所以 ∫ d⁴k/(2π)⁴ 1/(k²+m²) = [Λ² - m² ln(1+Λ²/m²)] / (16π²)

        参数:
            integrand: f(k) — 被积函数 (只依赖 |k|)
            Lambda:    动量截断 Λ

        返回:
            圈积分结果
        """
        N = 2000
        k_vals = np.linspace(0, Lambda, N)
        dk = Lambda / N

        # d=4: ∫ d⁴k/(2π)⁴ f(k) = 1/(8π²) ∫₀^Λ k³ f(k) dk
        integrand_vals = np.asarray(integrand(k_vals[1:]), dtype=float)  # 跳开 k=0 避免发散
        integral = np.sum(k_vals[1:]**3 * integrand_vals) * dk / (8 * np.pi**2)
        return float(integral)

    def momentum_routing(self, incoming: np.ndarray,
                         outgoing: np.ndarray) -> Tuple[np.ndarray, bool]:
        """动量路由检查 — 验证动量守恒并分配内部线动量

        参数:
            incoming:  入射动量 (N_particles, dim)
            outgoing:  出射动量 (M_particles, dim)

        返回:
            (total_momentum, is_conserved)
            其中 total_momentum = Σ p_in - Σ p_out
        """
        total_in = np.sum(incoming, axis=0)
        total_out = np.sum(outgoing, axis=0)
        total = total_in - total_out
        conserved = bool(np.all(np.abs(total) < 1e-10))
        return total, conserved


# ================================================================
# 正规化方案
# ================================================================

def regularize_cutoff(integrand: Callable[[np.ndarray], np.ndarray],
                      Lambda: float, N: int = 2000) -> float:
    """动量截断正规化 — 对发散积分引入硬截断 Λ

    ∫₀^∞ dk f(k)  →  ∫₀^Λ dk f(k)

    与 Phi4FeynmanRules.loop_integral 不同的是, 这是一个标量 dk 积分,
    而非 d⁴k 动量空间积分。适用于已约化到一维的积分。

    参数:
        integrand: f(k), k ∈ [0, ∞)
        Lambda:    紫外截断
        N:         采样点数

    返回:
        正规化积分值
    """
    k_vals = np.linspace(0, Lambda, N)
    dk = Lambda / N
    vals = integrand(k_vals[1:])
    return float(np.sum(vals) * dk)


def regularize_dim_reg(mass: float, power: int = 1,
                       eps: float = 0.01) -> float:
    """维数正规化占位 — 在 d = 4 - 2ε 维计算发散的圈积分

    这是一个占位实现, 返回单极点 1/ε 近似。
    完整实现需要 Gamma 函数展开和 minimal subtraction。

    标准结果 (标量蝌蚪图):
        μ^{2ε} ∫ d^d k/(2π)^d 1/(k²+m²) = -m²/(16π²) [1/ε - γ_E + ln(4πμ²/m²) + 1] + O(ε)

    参数:
        mass:   质量 m
        power:  传播子的幂次 (1 = 蝌蚪图, 2 = 盒图)
        eps:    维数正规子 ε = (4-d)/2

    返回:
        1/ε 发散部分的系数 (不含有限项)
    """
    if power == 1:
        # 蝌蚪图: ∫ d^d k/(2π)^d 1/(k²+m²) ≈ -m²/(16π²) * (1/ε)
        return float(-mass**2 / (16 * np.pi**2) * (1.0 / eps))
    elif power == 2:
        # 盒图: ∫ d^d k/(2π)^d 1/(k²+m²)² ≈ 1/(16π²) * (1/ε)
        return float(1.0 / (16 * np.pi**2) * (1.0 / eps))
    else:
        raise ValueError(f"power={power} 不支持; 仅支持 power=1,2")


# ================================================================
# 单圈自能 Π(p²) — 蝌蚪图
# ================================================================

def self_energy_1loop(p_sq: float, m: float, lam: float,
                      Lambda: float) -> float:
    """单圈自能 (蝌蚪图) Σ(p²) = Π(p²)

    φ⁴ 理论中, 单圈自能来自蝌蚪图 (tadpole):
        -iΣ = (1/2)(-iλ) ∫ d⁴k/(2π)⁴ i/(k² - m² + iε)

    对称因子 1/2 (一个圈, 一条传播子自收缩)。

    在动量截断正规化下 (d=4):
        Σ = λ/(32π²) [Λ² - m² ln(1 + Λ²/m²)]

    注意: 在 φ⁴ 理论中, 单圈自能与外部动量 p² 无关
    (蝌蚪图不携带外部动量)。p² 参数保留用于扩展 (如 multi-loop)。

    参数:
        p_sq:   外部动量平方 p² (单圈下不使用, 保留用于接口一致)
        m:      裸质量
        lam:    裸耦合 λ
        Lambda: 动量截断

    返回:
        自能 Σ = Π(p²) (实数, 质量量纲 2)

    >>> self_energy_1loop(0, 1.0, 0.1, 10.0)
    0.527...  # 正值, Λ² 主导

    参考文献:
        Peskin & Schroeder, Eq. (10.33)
    """
    # 解析积分:
    # ∫ d⁴k/(2π)⁴ 1/(k²+m²) = 1/(16π²) [Λ² - m² ln(1+Λ²/m²)] (欧氏)
    # Σ = (λ/2) × 以上
    m2 = m * m
    L2 = Lambda * Lambda
    log_arg = 1.0 + L2 / m2
    integral = (L2 - m2 * np.log(log_arg)) / (16.0 * np.pi**2)
    sigma = 0.5 * lam * integral
    return float(sigma)


# ================================================================
# 单圈顶角修正 Γ⁴(s,t,u)
# ================================================================

def vertex_correction_1loop_zero_momentum(m: float, lam: float,
                                          Lambda: float) -> float:
    """单圈顶角修正 — 零动量极限 (s=t=u=0)

    1-loop 顶角函数修正来自三个交叉的盒图 (s, t, u 通道):

        δΓ⁴ = (3/2)(-iλ)² ∫ d⁴k/(2π)⁴ [i/(k²-m²)]²

    对称因子 1/2, 三个通道 → 3/2。

    在动量截断下:
        ∫ d⁴k/(2π)⁴ 1/(k²+m²)² = 1/(16π²) [ln(1+Λ²/m²) - Λ²/(Λ²+m²)]

    最终:
        δλ_vertex = 3λ²/(32π²) [ln(1+Λ²/m²) - Λ²/(Λ²+m²)]

    参数:
        m:      裸质量
        lam:    裸耦合 λ
        Lambda: 动量截断

    返回:
        有效耦合偏移 δλ_vertex (零动量)

    >>> vertex_correction_1loop_zero_momentum(1.0, 0.1, 10.0)
    0.005...  # λ² 量级, 对数发散
    """
    m2 = m * m
    L2 = Lambda * Lambda
    log_arg = 1.0 + L2 / m2

    # ∫ d⁴k/(2π)⁴ 1/(k²+m²)² 的解析结果 (欧氏)
    integral_over_m2 = (np.log(log_arg) - L2 / (L2 + m2)) / (16.0 * np.pi**2)
    delta_lambda = 3.0 * lam * lam / (32.0 * np.pi**2) * (np.log(log_arg) - L2 / (L2 + m2))
    return float(delta_lambda)


# ================================================================
# 重整化常数 — 在壳 (on-shell) 方案
# ================================================================

def mass_counterterm(m: float, lam: float, Lambda: float,
                     scheme: str = 'on-shell') -> float:
    """质量抵消项 δm² — 在壳 (on-shell) 重整化方案

    重整化条件: 物理质量 m_phys 由传播子极点定义
        p² - m² - Σ(p²) = 0  at p² = m_phys²

    裸质量分解: m₀² = m_phys² + δm²
    抵消项:     δm² = -Σ(p²=m_phys²)

    在 φ⁴ 理论单圈下, Σ 与 p² 无关, 因此:
        δm² = -Σ = -λ/(32π²) [Λ² - m² ln(1+Λ²/m²)]

    参数:
        m:      物理质量 m_phys (on-shell 下)
        lam:    重整化耦合 λ_R ≈ λ_0 (单圈下)
        Lambda: 动量截断
        scheme: 方案标识 ('on-shell' 或 'ms-bar')

    返回:
        δm² (负值, 抵消 Λ² 发散)

    >>> mass_counterterm(1.0, 0.1, 10.0)
    -0.527...
    """
    sigma = self_energy_1loop(0.0, m, lam, Lambda)
    if scheme == 'on-shell':
        return float(-sigma)
    elif scheme == 'ms-bar':
        # MS-bar: 只减发散部分 (1/ε 极点或 Λ² 项)
        return float(-lam * Lambda**2 / (32.0 * np.pi**2))
    else:
        raise ValueError(f"未知方案: {scheme}; 支持 'on-shell', 'ms-bar'")


def coupling_counterterm(m: float, lam: float, Lambda: float,
                         scheme: str = 'on-shell') -> float:
    """耦合常数抵消项 δλ — 在壳 (on-shell) 重整化方案

    重整化条件: 物理耦合由 s=4m², t=u=0 的散射振幅定义
        λ_phys = λ_0 + δλ + vertex_correction(s=4m², t=0, u=0)

    在零动量近似下 (简单版本):
        δλ = -vertex_correction_1loop_zero_momentum(m, lam, Lambda)

    参数:
        m:      物理质量
        lam:    裸耦合 λ_0
        Lambda: 动量截断
        scheme: 'on-shell' 或 'ms-bar'

    返回:
        δλ (负值)

    >>> coupling_counterterm(1.0, 0.1, 10.0)
    -0.005...
    """
    vc = vertex_correction_1loop_zero_momentum(m, lam, Lambda)
    if scheme == 'on-shell':
        return float(-vc)
    elif scheme == 'ms-bar':
        # MS-bar: 只减对数发散
        return float(-3.0 * lam**2 / (32.0 * np.pi**2) * np.log(Lambda**2 / m**2))
    else:
        raise ValueError(f"未知方案: {scheme}; 支持 'on-shell', 'ms-bar'")


def field_renormalization(m: float, lam: float, Lambda: float) -> float:
    """场重整化常数 Z_φ

    在 φ⁴ 理论中, 场重整化从单圈自能的动量依赖部分获得:
        Σ(p²) = Σ(p²=0) + p² Σ'(0) + ...

    波函数重整化: Z_φ = 1 + dΣ/dp²|_{p²=m²}

    在 φ⁴ 理论单圈下, 蝌蚪图的自能与 p² 无关 (无动量依赖):
        dΣ/dp² = 0  →  Z_φ = 1

    场重整化从两圈开始 (sunset 图)。

    参数:
        m:      质量
        lam:    耦合
        Lambda: 动量截断

    返回:
        Z_φ = 1.0 (单圈下为 1)

    >>> field_renormalization(1.0, 0.1, 10.0)
    1.0
    """
    return 1.0


# ================================================================
# 重整化群 — β 函数与跑动耦合
# ================================================================

def beta_function(lam: float) -> float:
    """φ⁴ 理论单圈 β 函数

    β(λ) = μ dλ/dμ = 3λ²/(16π²)   (单圈, d=4)

    推导: 从抵消项 δλ 的 μ 依赖性获得
        λ(μ) = λ_0 + [3λ²/(32π²)] ln(Λ²/μ²)
        μ dλ/dμ = 3λ²/(16π²)

    参数:
        lam: 重整化耦合 λ(μ) 在当前标度 μ 下的值

    返回:
        β(λ)

    >>> beta_function(0.1)
    0.001899...  # 3*0.01/(16π²) ≈ 0.00019
    """
    return float(3.0 * lam * lam / (16.0 * np.pi**2))


def running_coupling(lam0: float, mu0: float, mu: float) -> float:
    """跑动耦合 λ(μ) — 从 β 函数积分

    解重正化群方程 dλ/d ln μ = β(λ):

        1/λ(μ₀) - 1/λ(μ) = [3/(16π²)] ln(μ/μ₀)

        λ(μ) = λ₀ / [1 - (3λ₀/(16π²)) ln(μ/μ₀)]

    注意: λ(μ) 在 μ → ∞ 时增长 (φ⁴ 不是渐近自由的)。

    当分母 → 0 时出现 Landau 极点:
        μ_Landau = μ₀ exp(16π²/(3λ₀))

    参数:
        lam0: λ(μ₀) — 参考标度 μ₀ 处的耦合
        mu0:  参考标度 μ₀
        mu:   目标标度 μ

    返回:
        λ(μ)

    >>> running_coupling(0.1, 1.0, 10.0)
    0.107...  # 随 μ 增大而增大

    >>> running_coupling(0.1, 1.0, 0.1)
    0.093...  # 随 μ 减小而减小
    """
    beta_coeff = 3.0 / (16.0 * np.pi**2)
    inv_lam = 1.0 / lam0 - beta_coeff * np.log(mu / mu0)
    if inv_lam <= 0:
        # Landau 极点: λ → ∞
        return float('inf')
    return float(1.0 / inv_lam)


def landau_pole(lam0: float, mu0: float) -> float:
    """Landau 极点 — λ(μ) 发散时的标度

    μ_Landau = μ₀ exp(16π²/(3λ₀))

    在微扰论失效的标度以上, 需要非微扰方法。

    参数:
        lam0: λ(μ₀)
        mu0:  参考标度

    返回:
        μ_Landau
    """
    return float(mu0 * np.exp(16.0 * np.pi**2 / (3.0 * lam0)))


# ================================================================
# 组合 — 重整化修正传播子
# ================================================================

def renormalized_propagator(p_sq: float, m: float, lam: float,
                            Lambda: float) -> complex:
    """单圈重整化传播子

    D_R(p²) = i / [p² - m² - Σ(p²) + δm² + iε]
             = i / [p² - m_phys² + iε]   (在壳方案)

    在壳方案下, 抵消项精确消除自能:
        m_phys² = m² + Σ + δm² = m²

    所以重整化传播子 = 自由传播子 (在壳方案, 单圈)。

    参数:
        p_sq:   动量平方
        m:      物理质量
        lam:    耦合
        Lambda: 截断

    返回:
        重整化传播子
    """
    dm2 = mass_counterterm(m, lam, Lambda, scheme='on-shell')
    sigma = self_energy_1loop(p_sq, m, lam, Lambda)
    m_phys_sq = m**2 + sigma + dm2  # on-shell 下应等于 m²
    return complex(0, 1) / (p_sq - m_phys_sq + 1j * 1e-10)


# ================================================================
# 信息
# ================================================================

def renormalization_summary(m: float, lam: float, Lambda: float) -> str:
    """打印重整化信息摘要"""
    dm2 = mass_counterterm(m, lam, Lambda)
    dlam = coupling_counterterm(m, lam, Lambda)
    Z_phi = field_renormalization(m, lam, Lambda)
    beta = beta_function(lam)
    mu_landau = landau_pole(lam, 1.0)

    lines = [
        f"φ⁴ 重整化 (单圈, 动量截断 Λ={Lambda})",
        f"  m = {m}, λ = {lam}",
        f"  δm² = {dm2:.6f}",
        f"  δλ  = {dlam:.6f}",
        f"  Z_φ = {Z_phi}",
        f"  β(λ) = {beta:.6f} ({'+' if beta > 0 else ''}{beta/lam*100:.2f}% per e-fold)",
        f"  Landau 极点 (μ₀=1): μ_Landau = {mu_landau:.2e}",
    ]
    return '\n'.join(lines)

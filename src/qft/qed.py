"""QED 散射过程 — Compton, 对湮灭, Møller 散射

精细结构常数: α = e²/(4π) ≈ 1/137.036

Mandelstam 不变量 (2→2 散射, p₁ + p₂ → p₃ + p₄):
    s = (p₁ + p₂)² = 4(p_cm² + m²)
    t = (p₁ - p₃)² = -2p_cm²(1 - cosθ)
    u = (p₁ - p₄)² = -2p_cm²(1 + cosθ)
    满足: s + t + u = 4m²

截面公式:
    dσ/dΩ = |M|² / (64π² s) · |p_f|/|p_i|   (质心系, 等质量)

包含的过程:
    1. Compton 散射:  γ + e⁻ → γ + e⁻  (Klein-Nishina 公式)
    2. 对湮灭:        e⁺ + e⁻ → μ⁺ + μ⁻  (s 道)
    3. Møller 散射:   e⁻ + e⁻ → e⁻ + e⁻  (t+u 道)

参考文献:
    Peskin & Schroeder, 第 5 章
    Halzen & Martin, 第 7 章
    Klein-Nishina: Z. Phys. 52, 853 (1929)
"""

import numpy as np
from typing import Tuple


# ================================================================
# 常数
# ================================================================

ALPHA_QED = 1.0 / 137.035999084  # α = e²/4π
ALPHA_SQ = ALPHA_QED ** 2

M_E = 0.511  # 电子质量 MeV/c²
M_MU = 105.658  # μ 子质量 MeV/c²


# ================================================================
# Mandelstam 不变量
# ================================================================

def mandelstam(p1_4vec: np.ndarray, p2_4vec: np.ndarray,
               p3_4vec: np.ndarray, p4_4vec: np.ndarray,
               mass: float = 0.0) -> Tuple[float, float, float]:
    """计算 2→2 散射的 Mandelstam 不变量

        s = (p₁ + p₂)²
        t = (p₁ - p₃)²
        u = (p₁ - p₄)²

    参数:
        p1, p2, p3, p4: 四动量 (E, px, py, pz)
        mass:            粒子质量 (用于验证 s+t+u = 4m²)

    返回:
        (s, t, u) 标量值
    """
    # Lorentz 内积: p·q = p₀q₀ - p₁q₁ - p₂q₂ - p₃q₃
    def dot(p, q):
        return p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]

    s = dot(p1_4vec + p2_4vec, p1_4vec + p2_4vec)
    t = dot(p1_4vec - p3_4vec, p1_4vec - p3_4vec)
    u = dot(p1_4vec - p4_4vec, p1_4vec - p4_4vec)

    return float(s), float(t), float(u)


def mandelstam_from_cm(s: float, theta_cm: float,
                        mass: float = 0.0) -> Tuple[float, float, float]:
    """从质心能量平方 s 和散射角 θ 计算 Mandelstam 不变量

    质心系 2→2 散射 (等质量):
        s = 4(E_cm²) = 4(p² + m²)
        t = -2p²(1 - cosθ) = m² - s/2 + (s/2)√(1 - 4m²/s) cosθ   [简化]
           ≈ -s/2 (1 - cosθ)  (无质量极限)
        u = -2p²(1 + cosθ) ≈ -s/2 (1 + cosθ)

    参数:
        s:        质心能量平方
        theta_cm: 质心散射角 (弧度)
        mass:     粒子质量

    返回:
        (s, t, u)
    """
    s_min = 4 * mass**2
    if s <= s_min:
        s = s_min * 1.001

    p_sq = s / 4 - mass**2  # 质心动量平方

    if p_sq < 0:
        p_sq = 0.0

    cos_theta = np.cos(theta_cm)
    t = mass**2 + mass**2 - s / 2 + 2 * p_sq * cos_theta
    u = mass**2 + mass**2 - s / 2 - 2 * p_sq * cos_theta

    # 化简形式:
    # t = 2*mass**2 - s/2 + s/2 * beta**2 * cos_theta
    # u = 2*mass**2 - s/2 - s/2 * beta**2 * cos_theta
    # 其中 beta² = 1 - 4m²/s

    return float(s), float(t), float(u)


# ================================================================
# 相空间因子
# ================================================================

def phase_space_factor(s: float, m3: float, m4: float,
                       m1: float = 0.0, m2: float = 0.0) -> float:
    """2→2 散射相空间因子 |p_f| / (16π |p_i| s)

    质心系:
        dσ/dΩ = (1/64π² s) · (|p_f|/|p_i|) · |M|²

    返回:
        |p_f| / (64π² s |p_i|)  (等质量时 |p_f| = |p_i| → 1/(64π² s))
    """
    s_min = (m1 + m2)**2
    if s <= s_min:
        return 0.0

    # 初态质心动量
    p_i_sq = (s - (m1 + m2)**2) * (s - (m1 - m2)**2) / (4 * s)
    p_i_sq = max(p_i_sq, 0.0)

    # 末态质心动量
    p_f_sq = (s - (m3 + m4)**2) * (s - (m3 - m4)**2) / (4 * s)
    p_f_sq = max(p_f_sq, 0.0)

    if p_i_sq < 1e-30:
        return 0.0

    factor = np.sqrt(p_f_sq / p_i_sq) / (64 * np.pi ** 2 * s)
    return float(factor)


# ================================================================
# Compton 散射 (Klein-Nishina 公式)
# ================================================================

def compton_amplitude_squared(s: float, u: float,
                               m_e: float = M_E) -> float:
    """Compton 散射自旋平均振幅平方 |M|² (γ + e⁻ → γ + e⁻)

    Klein-Nishina 振幅:
        |M|² = 2e⁴ [ -u/s - s/u + 4m²(1/s + 1/u) - 4m⁴(1/s + 1/u)² ]
              = 32π² α² [ -u/s - s/u + 4m²(1/s + 1/u) - 4m⁴(1/s + 1/u)² ]

    注意: 使用 Peskin 约定, 其中 e² = 4πα

    参数:
        s:    Mandelstam s
        u:    Mandelstam u
        m_e:  电子质量

    返回:
        |M|² (无量纲, 含 e⁴ 耦合常数)
    """
    if abs(s) < 1e-15 or abs(u) < 1e-15:
        return 0.0

    inv_s = 1.0 / s
    inv_u = 1.0 / u

    # Peskin (5.87): |M|² = 2e⁴ [p·k' / p·k + p·k / p·k' + 2m²(1/p·k - 1/p·k') + m⁴(1/p·k - 1/p·k')²]
    # 其中 p·k = -(s - m²)/2, p·k' = -(u - m²)/2
    # 等价形式:
    term1 = -u / s - s / u
    term2 = 4 * m_e**2 * (inv_s + inv_u)
    term3 = -4 * m_e**4 * (inv_s + inv_u)**2

    # e⁴ = (4πα)² = 16π² α²
    e4 = (4 * np.pi * ALPHA_QED) ** 2
    result = 2 * e4 * (term1 + term2 + term3)
    return float(result)


def compton_cross_section(omega_lab: float, theta_lab: float,
                           m_e: float = M_E) -> float:
    """Compton 散射微分截面 dσ/dΩ (Klein-Nishina 公式)

    实验室系 (电子静止):
        dσ/dΩ = α²/(2m²) · (ω'/ω)² · [ω/ω' + ω'/ω - sin²θ]

    其中 Compton 公式:
        ω' = ω / [1 + (ω/m)(1 - cosθ)]

    参数:
        omega_lab: 入射光子能量 (实验室系)
        theta_lab: 散射角 (实验室系, 弧度)
        m_e:       电子质量

    返回:
        dσ/dΩ (自然单位, 量纲 [E]^{-2})

    参考文献:
        Klein & Nishina, Z. Phys. 52, 853 (1929)
    """
    # 散射光子能量 (Compton 公式)
    ratio = 1.0 + (omega_lab / m_e) * (1.0 - np.cos(theta_lab))
    omega_prime = omega_lab / ratio

    # Klein-Nishina 公式
    factor = (ALPHA_QED**2) / (2 * m_e**2)
    cross_section_ratio = (omega_prime / omega_lab) ** 2
    angular = omega_lab / omega_prime + omega_prime / omega_lab - np.sin(theta_lab)**2

    return float(factor * cross_section_ratio * angular)


def compton_total_cross_section(omega_lab: float,
                                  m_e: float = M_E) -> float:
    """Compton 散射总截面 (数值积分)

    参数:
        omega_lab: 入射光子能量
        m_e:       电子质量

    返回:
        σ_total (自然单位)
    """
    n_theta = 200
    theta_vals = np.linspace(0, np.pi, n_theta)

    dsigma = np.array([compton_cross_section(omega_lab, t, m_e)
                        for t in theta_vals])
    integrand = dsigma * 2 * np.pi * np.sin(theta_vals)

    return float(np.trapezoid(integrand, theta_vals))


# ================================================================
# 对湮灭 e⁺e⁻ → μ⁺μ⁻
# ================================================================

def pair_annihilation_amplitude_squared(s: float, t: float, u: float,
                                          m_e: float = M_E,
                                          m_mu: float = M_MU) -> float:
    """e⁺e⁻ → μ⁺μ⁻ 自旋平均振幅平方 |M|²

    纯 QED, s 道光子交换:
        |M|² = (e⁴/s²) · 4 · [(p₁·p₃)(p₂·p₄) + (p₁·p₄)(p₂·p₃)
                              + m_e²(m_μ² - p₂·p₃ - ...) ... ]
    简化 (忽略 m_e):
        |M|² = 2e⁴ (t² + u²) / s²    (m_e → 0)

    参数:
        s, t, u: Mandelstam 不变量
        m_e:     电子质量
        m_mu:     μ 子质量

    返回:
        |M|²  (含 e⁴)
    """
    e4 = (4 * np.pi * ALPHA_QED) ** 2

    # 阈值检查: s < 4*m_mu² → 无法产生 μ⁺μ⁻
    s_threshold = 4 * m_mu ** 2
    if s < s_threshold:
        return 0.0

    if m_e < 1e-10 and m_mu < 1e-10:
        # 无质量极限: |M|² = 2e⁴ (t² + u²) / s²
        if abs(s) < 1e-15:
            return 0.0
        return float(2 * e4 * (t**2 + u**2) / s**2)

    # 带质量完整公式 (Peskin 5.12):
    # |M|² = e⁴/s² · [2(t² + u² + 4s(m_e²+m_μ²) - 2(m_e²+m_μ²)²)]
    if abs(s) < 1e-15:
        return 0.0

    term = 2 * (t**2 + u**2 + 4 * s * (m_e**2 + m_mu**2) -
                2 * (m_e**2 + m_mu**2)**2)
    return float(e4 * term / s**2)


def pair_annihilation_cross_section(s: float, theta_cm: float,
                                      m_e: float = M_E,
                                      m_mu: float = M_MU) -> float:
    """e⁺e⁻ → μ⁺μ⁻ 微分散射截面 dσ/dΩ

    在 s ≫ m_e², m_μ² 时:
        dσ/dΩ = (α²/4s) (1 + cos²θ)

    参数:
        s:        质心能量平方
        theta_cm: 质心散射角 (弧度)
        m_e:      电子质量
        m_mu:     μ 子质量

    返回:
        dσ/dΩ  (自然单位)

    参考文献:
        Peskin & Schroeder, (5.13)
    """
    # 相空间因子
    ps = phase_space_factor(s, m_mu, m_mu, m_e, m_e)
    if ps < 1e-30:
        return 0.0

    # Mandelstam 不变量
    _, t, u = mandelstam_from_cm(s, theta_cm, mass=m_e)

    # 振幅平方
    amp2 = pair_annihilation_amplitude_squared(s, t, u, m_e, m_mu)

    # dσ/dΩ = ps × |M|²  (注意 ps 已含 1/(64π²s) 因子)
    # 实际: ps = (|p_f|/|p_i|) / (64π²s)
    # dσ/dΩ = ps × |M|²
    return float(amp2 / (64 * np.pi**2 * s))


def pair_annihilation_total_cross_section(s: float,
                                            m_e: float = M_E,
                                            m_mu: float = M_MU) -> float:
    """e⁺e⁻ → μ⁺μ⁻ 总截面 (数值积分)

    参数:
        s:    质心能量平方
        m_e:  电子质量
        m_mu: μ 子质量

    返回:
        σ_total
    """
    n_theta = 200
    theta_vals = np.linspace(0, np.pi, n_theta)
    dsigma = np.array([pair_annihilation_cross_section(s, t, m_e, m_mu)
                        for t in theta_vals])
    integrand = dsigma * 2 * np.pi * np.sin(theta_vals)
    return float(np.trapezoid(integrand, theta_vals))


# ================================================================
# Møller 散射 e⁻e⁻ → e⁻e⁻
# ================================================================

def moller_amplitude_squared(s: float, t: float, u: float,
                               m_e: float = M_E) -> float:
    """Møller 散射自旋平均振幅平方 |M|² (e⁻e⁻ → e⁻e⁻)

    包含 t 道和 u 道 (交换图), 以及 t-u 干涉:
        |M|² = e⁴ [ s²+u²    s²+t²    2s²       4m²s     8m⁴s²  ]
             = e⁴ [ ────── + ────── - ──── + ──────── + ──────── ]
                  [   t²       u²       tu    (t u)     t² u²    ]

    无质量极限 (m_e → 0):
        |M|² = 2e⁴ (s² + u²)/t² + 2e⁴ (s² + t²)/u² + 4e⁴ s²/(t u)

    参数:
        s, t, u: Mandelstam 不变量
        m_e:     电子质量

    返回:
        |M|²

    参考文献:
        Møller, Ann. Phys. 14, 531 (1932)
        Halzen & Martin, (6.60)
    """
    e4 = (4 * np.pi * ALPHA_QED) ** 2

    # 避免除零
    if abs(t) < 1e-15 or abs(u) < 1e-15:
        return 0.0

    # 无质量极限
    if m_e < 1e-10:
        term_t = 2 * (s**2 + u**2) / t**2
        term_u = 2 * (s**2 + t**2) / u**2
        term_interference = 4 * s**2 / (t * u)
        return float(e4 * (term_t + term_u + term_interference) / 4)

    # 带质量完整公式 (Peskin 习题 5.3)
    m2 = m_e**2
    term_t = (s**2 + u**2 + 4 * m2 * (t - m2)) / t**2
    term_u = (s**2 + t**2 + 4 * m2 * (u - m2)) / u**2
    term_interference = (2 * s**2 - 8 * m2 * s + 12 * m2**2) / (t * u)

    return float(e4 * (term_t + term_u + term_interference))


def moller_cross_section(s: float, theta_cm: float,
                           m_e: float = M_E) -> float:
    """Møller 散射微分散射截面 dσ/dΩ

    无质量极限:
        dσ/dΩ = α² (3 + cos²θ)² / (2s sin⁴θ)

    带质量 (s ≫ m²):
        dσ/dΩ = α²/(4m²) · [用 |M|² 完整公式]

    参数:
        s:        质心能量平方
        theta_cm: 质心散射角
        m_e:      电子质量

    返回:
        dσ/dΩ

    参考文献:
        Halzen & Martin, 第 6 章
    """
    # 相空间因子
    ps = phase_space_factor(s, m_e, m_e, m_e, m_e)
    if ps < 1e-30:
        return 0.0

    # Mandelstam 不变量
    _, t, u = mandelstam_from_cm(s, theta_cm, mass=m_e)

    # 振幅平方
    amp2 = moller_amplitude_squared(s, t, u, m_e)

    return float(amp2 / (64 * np.pi**2 * s))


def moller_total_cross_section(s: float, m_e: float = M_E,
                                 theta_min: float = 0.01) -> float:
    """Møller 散射总截面 (数值积分, 避开前向发散)

    参数:
        s:         质心能量平方
        m_e:       电子质量
        theta_min: 最小散射角 (规避 t=0 发散)

    返回:
        σ_total
    """
    n_theta = 200
    theta_vals = np.linspace(theta_min, np.pi - theta_min, n_theta)
    dsigma = np.array([moller_cross_section(s, t, m_e)
                        for t in theta_vals])
    integrand = dsigma * 2 * np.pi * np.sin(theta_vals)
    return float(np.trapezoid(integrand, theta_vals))


# ================================================================
# 信息
# ================================================================

def qed_summary() -> str:
    """QED 模块概览"""
    lines = [
        "QED Scattering Module",
        f"  α = 1/{1/ALPHA_QED:.0f}",
        f"  m_e = {M_E} MeV/c²",
        f"  m_μ = {M_MU} MeV/c²",
        "",
        "  Processes:",
        "    compton_cross_section(ω_lab, θ_lab)     → dσ/dΩ (Klein-Nishina)",
        "    pair_annihilation_cross_section(s, θ_cm) → dσ/dΩ (e⁺e⁻→μ⁺μ⁻)",
        "    moller_cross_section(s, θ_cm)            → dσ/dΩ (e⁻e⁻→e⁻e⁻)",
        "",
        "  Mandelstam:",
        "    mandelstam(p1,p2,p3,p4) → (s, t, u)",
        "    mandelstam_from_cm(s, θ_cm) → (s, t, u)",
    ]
    return '\n'.join(lines)

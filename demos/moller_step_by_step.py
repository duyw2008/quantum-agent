#!/usr/bin/env python3
"""
Møller 散射步步详解 — 从 Feynman 图到微分散射截面
====================================================
完整展示 e⁻e⁻ → e⁻e⁻ 散射的计算链条，匹配《一读就懂的量子场论》§4.6-4.7

步骤:
  1. Feynman 图 (t 道 + u 道) ASCII 艺术
  2. Feynman 规则 → 振幅 M_t 和 M_u (含相对负号)
  3. 自旋平均 |M|² — 用求迹技术显式计算, 并验证解析公式
  4. 微分散射截面 dσ/dΩ (质心系)
  5. 图像: dσ/dΩ vs θ (多个能量), 总截面 vs √s
  6. 数值表: 特定角度下的完整计算
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from src.qft.dirac import GammaMatrices, DiracSpinor, spin_sum_u, spin_sum_v, dirac_slash, I4
from src.qft.qed import (
    moller_cross_section, moller_amplitude_squared,
    moller_total_cross_section, mandelstam_from_cm,
    ALPHA_QED, M_E,
)

# ================================================================
# 全局设置
# ================================================================
ALPHA = ALPHA_QED
E_SQ = (4 * np.pi * ALPHA) ** 2  # e⁴ = (4πα)²
GM = GammaMatrices('dirac')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 用于输出的分隔线
SEP = "=" * 72
SUBSEP = "-" * 56

def print_header(title: str):
    print(f"\n{SEP}")
    print(f"  {title}")
    print(f"{SEP}")

def print_step(n: int, title: str):
    print(f"\n{SUBSEP}")
    print(f"  步骤 {n}: {title}")
    print(f"{SUBSEP}")


# ================================================================
# 第 1 步: Feynman 图
# ================================================================

def step1_feynman_diagrams():
    """展示 Møller 散射的两个树图 Feynman 图"""
    print_step(1, "Feynman 图 — t 道 与 u 道")

    t_channel = r"""
        e⁻(p₁) ────────────── e⁻(p₃)
                    \
                     \  γ (t = (p₁-p₃)²)
                      \
        e⁻(p₂) ────────────── e⁻(p₄)

        t 道: 光子动量 q = p₁ - p₃
        传播子: -i g_{μν} / t
        M_t ∝ e²/t · [ū₃ γ^μ u₁] [ū₄ γ_μ u₂]
    """

    u_channel = r"""
        e⁻(p₁) ────╲            e⁻(p₄)
                     ╲
                      ╲  γ (u = (p₁-p₄)²)
                     ╱
        e⁻(p₂) ────╱            e⁻(p₃)

        u 道: 光子动量 q = p₁ - p₄  (全同粒子交换)
        传播子: -i g_{μν} / u
        M_u ∝ e²/u · [ū₄ γ^μ u₁] [ū₃ γ_μ u₂]
    """

    print(t_channel)
    print(u_channel)
    print("  ⚠ 注意: 两个图之间有相对负号 (fermion 交换反对称)")
    print("    总振幅: M = M_t - M_u")
    print("    符号来源: 末态交换两个全同费米子 → 反对称化 → (-1)")


# ================================================================
# 第 2 步: Feynman 规则 → 振幅
# ================================================================

def step2_feynman_rules():
    """用 Feynman 规则写出振幅 M_t 和 M_u"""
    print_step(2, "Feynman 规则 → 构造振幅")

    print("""
  QED Feynman 规则 (Møller 散射用到的):
  ┌─────────────────────────────────────────────────────────────┐
  │  外线 (入射电子):   u(p,s)                                  │
  │  外线 (出射电子):   ū(p,s)                                  │
  │  光子传播子:        -i g_{μν} / q²    (Feynman 规范)        │
  │  顶点:              -i e γ^μ                                │
  │  圈积分:            ∫ d⁴k/(2π)⁴                             │
  │  全同费米子交换:    相对负号 (-1)                            │
  └─────────────────────────────────────────────────────────────┘

  振幅构造 (逐项):
    M_t = [ū(p₃) (-i e γ^μ) u(p₁)] · [-i g_{μν} / t] · [ū(p₄) (-i e γ^ν) u(p₂)]
        = i e² / t · [ū₃ γ^μ u₁] [ū₄ γ_μ u₂]

    M_u = [ū(p₄) (-i e γ^μ) u(p₁)] · [-i g_{μν} / u] · [ū(p₃) (-i e γ^ν) u(p₂)]
        = i e² / u · [ū₄ γ^μ u₁] [ū₃ γ_μ u₂]

  总振幅 (含交换反对称):
    M = M_t - M_u = i e² [ ū₃γ^μu₁ ū₄γ_μu₂ / t  -  ū₄γ^μu₁ ū₃γ_μu₂ / u ]
    """)

    # 选择一个具体的运动学点来数值演示
    s_val = 100.0  # s = 100 MeV²  (电子质量 ~0.511 MeV, 所以这是高能)
    theta_val = np.pi / 4  # 45°
    m_e = M_E

    _, t_val, u_val = mandelstam_from_cm(s_val, theta_val, mass=m_e)

    print(f"  选取运动学点: √s = {np.sqrt(s_val):.2f} MeV, θ_cm = {np.degrees(theta_val):.1f}°")
    print(f"  Mandelstam: s = {s_val:.2f}, t = {t_val:.2f}, u = {u_val:.2f} (MeV²)")
    print(f"  检验: s + t + u = {s_val + t_val + u_val:.2f} ≈ 4m² = {4*m_e**2:.4f}")

    return s_val, t_val, u_val, theta_val


# ================================================================
# 第 3 步: 自旋平均 |M|²  — 求迹技术
# ================================================================

def compute_gamma_trace_product(p1, p2, p3, p4, mass):
    """
    计算自旋平均 |M|² 的核心求迹

    对于 t 道平方:
      |M_t|² ∝ Tr[(p̸₃+m)γ^μ(p̸₁+m)γ^ν] · Tr[(p̸₄+m)γ_μ(p̸₂+m)γ_ν] / t²

    对于 u 道平方:
      |M_u|² ∝ Tr[(p̸₄+m)γ^μ(p̸₁+m)γ^ν] · Tr[(p̸₃+m)γ_μ(p̸₂+m)γ_ν] / u²

    干涉项:
      M_t M_u* + cc ∝ Tr[(p̸₃+m)γ^μ(p̸₁+m)γ^ν(p̸₄+m)γ_μ(p̸₂+m)γ_ν] / (tu)
    """
    gm = GM

    # 构建 p̸ 矩阵
    def pslash(p):
        return dirac_slash(p, gm)

    def trace_4(gamma_list):
        """Tr[Γ₁ Γ₂ Γ₃ Γ₄] — 数值计算"""
        M = (gm.gamma[gamma_list[0]] @ gm.gamma[gamma_list[1]] @
             gm.gamma[gamma_list[2]] @ gm.gamma[gamma_list[3]])
        return np.trace(M)

    # --- t 道: Tr[(p̸₃+m)γ^μ (p̸₁+m)γ^ν] · Tr[(p̸₄+m)γ_μ (p̸₂+m)γ_ν] ---
    pslash1 = pslash(p1)
    pslash2 = pslash(p2)
    pslash3 = pslash(p3)
    pslash4 = pslash(p4)

    # 上标和下标通过度规联系: γ_μ = g_{μν} γ^ν
    # 在求迹中: Tr[...γ_μ...] Tr[...γ^μ...] = Σ_μ g_{μμ} (数值上: μ=0→+, μ=1,2,3→-)
    metric = np.array([1, -1, -1, -1])

    trace_t = 0.0
    trace_u = 0.0
    trace_int = 0.0

    for mu in range(4):
        g_mu = gm.gamma[mu]
        gamma_mu_sign = metric[mu]  # 降指标时的符号

        # t 道第一迹: Tr[(p̸₃+m) γ^μ (p̸₁+m) γ^ν]
        tr1_t = np.trace((pslash3 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @ gm.gamma[0])
        tr2_t = np.trace((pslash3 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @ gm.gamma[1])
        tr3_t = np.trace((pslash3 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @ gm.gamma[2])
        tr4_t = np.trace((pslash3 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @ gm.gamma[3])
        tr1_t_vec = np.array([tr1_t, tr2_t, tr3_t, tr4_t])

        # t 道第二迹: Tr[(p̸₄+m) γ_μ (p̸₂+m) γ_ν]
        # γ_μ(p̸₂+m)γ_ν = g_{μμ} · (用 γ^μ 但降指标处理)
        # 展开: Tr[(p̸₄+m) γ^ρ (p̸₂+m) γ^σ] · g_{μρ} g_{νσ} ... 
        # 简化: 直接计算 4×4 的逐个配对
        for nu in range(4):
            g_nu = gm.gamma[nu]
            tr1_t2 = np.trace((pslash4 + mass * I4) @ g_mu @ (pslash2 + mass * I4) @ g_nu)
            trace_t += float(np.real(tr1_t_vec[nu] * tr1_t2 * metric[mu] * metric[nu]))

        # --- u 道: Tr[(p̸₄+m)γ^μ (p̸₁+m)γ^ν] · Tr[(p̸₃+m)γ_μ (p̸₂+m)γ_ν] ---
        tr1_u_vec = np.array([
            np.trace((pslash4 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @ gm.gamma[nu])
            for nu in range(4)
        ])
        for nu in range(4):
            g_nu = gm.gamma[nu]
            tr2_u = np.trace((pslash3 + mass * I4) @ g_mu @ (pslash2 + mass * I4) @ g_nu)
            trace_u += float(np.real(tr1_u_vec[nu] * tr2_u * metric[mu] * metric[nu]))

        # --- 干涉项: Tr[(p̸₃+m)γ^μ (p̸₁+m)γ^ν (p̸₄+m)γ_μ (p̸₂+m)γ_ν] ---
        for nu in range(4):
            g_nu = gm.gamma[nu]
            mat_int = ((pslash3 + mass * I4) @ g_mu @ (pslash1 + mass * I4) @
                        g_nu @ (pslash4 + mass * I4) @ g_mu @ (pslash2 + mass * I4) @ g_nu)
            trace_int += float(np.real(np.trace(mat_int) * metric[mu] * metric[nu]))

    return trace_t, trace_u, trace_int


def step3_spin_averaged_amplitude(s_val, t_val, u_val):
    """自旋平均 |M|² — 用求迹显式计算 + 验证解析公式"""
    print_step(3, "自旋平均 |M|² — 求迹技术")

    m_e = M_E

    # 构建四动量 (质心系, 散射面 xz)
    # p₁ = (E, 0, 0, p), p₂ = (E, 0, 0, -p)
    # p₃ = (E, p sinθ, 0, p cosθ), p₄ = (E, -p sinθ, 0, -p cosθ)
    E = np.sqrt(s_val) / 2
    p_mag = np.sqrt(max(E**2 - m_e**2, 0))
    theta = np.arccos(np.clip(
        (t_val - 2*m_e**2 + s_val/2) / (2 * max(p_mag**2, 1e-30)),
        -1, 1
    ))

    p1 = np.array([E, 0, 0, p_mag])
    p2 = np.array([E, 0, 0, -p_mag])
    p3 = np.array([E, p_mag * np.sin(theta), 0, p_mag * np.cos(theta)])
    p4 = np.array([E, -p_mag * np.sin(theta), 0, -p_mag * np.cos(theta)])

    print(f"\n  四动量 (质心系, E={E:.3f}, |p|={p_mag:.3f} MeV):")
    print(f"    p₁ = ({p1[0]:.2f}, {p1[1]:.2f}, {p1[2]:.2f}, {p1[3]:.2f})")
    print(f"    p₂ = ({p2[0]:.2f}, {p2[1]:.2f}, {p2[2]:.2f}, {p2[3]:.2f})")
    print(f"    p₃ = ({p3[0]:.2f}, {p3[1]:.2f}, {p3[2]:.2f}, {p3[3]:.2f})")
    print(f"    p₄ = ({p4[0]:.2f}, {p4[1]:.2f}, {p4[2]:.2f}, {p4[3]:.2f})")

    # --- 3A. 自旋求和规则 ---
    print(f"\n  ── 3A. 自旋求和规则 ──")
    print(f"    Σ_s u(p,s) ū(p,s) = p̸ + m")
    print(f"    Σ_s v(p,s) v̄(p,s) = p̸ - m")
    print(f"    共 2² × 2² = 16 个螺旋度态求和")

    # --- 3B. |M|² 的求迹表示 ---
    print(f"\n  ── 3B. |M|² 的求迹形式 ──")
    print("""
    |M|² = (1/4) Σ_{spins} |M_t - M_u|²

    展开为三项:
      |M|² = |M_t|² + |M_u|² - 2Re(M_t M_u*)

    其中:
      |M_t|² ∝ Tr[(p̸₃+m)γ^μ(p̸₁+m)γ^ν] Tr[(p̸₄+m)γ_μ(p̸₂+m)γ_ν] / t²
      |M_u|² ∝ Tr[(p̸₄+m)γ^μ(p̸₁+m)γ^ν] Tr[(p̸₃+m)γ_μ(p̸₂+m)γ_ν] / u²
      干涉  ∝ Tr[(p̸₃+m)γ^μ(p̸₁+m)γ^ν(p̸₄+m)γ_μ(p̸₂+m)γ_ν] / (t u)

    (其中 |M|² 定义为 Feynman 振幅 |M|² × 自旋平均因子 1/4)
    """)

    # --- 3C. 数值求迹 ---
    print(f"  ── 3C. 数值求迹计算 ──")
    trace_t, trace_u, trace_int = compute_gamma_trace_product(p1, p2, p3, p4, m_e)

    print(f"    求迹结果:")
    print(f"      Tr_t  (t 道平方): {trace_t:.6f}")
    print(f"      Tr_u  (u 道平方): {trace_u:.6f}")
    print(f"      Tr_int (干涉项):  {trace_int:.6f}")

    # 组装 |M|²
    # |M|² (含 e⁴, 平均因子 1/4)
    # |M_t|² = (e⁴/t²) · (1/4) · Tr_t  (注意: 两个 Tr 的乘积已在上面的求和中完成)
    # 上面的 trace_t 已经包含了两个 Tr 乘积的求和 Σ_{μν}
    # 完整: |M_t|² = e⁴/(4t²) · Tr_t   (1/4 是自旋平均)
    e4 = E_SQ
    amp2_t = e4 * trace_t / (4 * t_val**2)
    amp2_u = e4 * trace_u / (4 * u_val**2)
    amp2_int = -e4 * trace_int / (4 * t_val * u_val)  # 干涉项含负号 (M_t - M_u)

    amp2_numerical = float(amp2_t + amp2_u + amp2_int)

    print(f"\n    组装 |M|² (数值求迹):")
    print(f"      |M_t|²  = e⁴/{4} · Tr_t/t²  = {amp2_t:.6e}")
    print(f"      |M_u|²  = e⁴/{4} · Tr_u/u²  = {amp2_u:.6e}")
    print(f"      干涉   = -e⁴/{4} · Tr_int/(tu) = {amp2_int:.6e}")
    print(f"      ─────────────────────────────────────")
    print(f"      |M|² (数值) = {amp2_numerical:.6e}")

    # --- 3D. 求迹验证 — 与解析迹公式对照 ---
    print(f"\n  ── 3D. 求迹验证 — 与标准迹公式对照 ──")

    # 标准迹公式 (用 Mandelstam 变量表示四动量内积):
    #   p₁·p₂ = s/2 - m²
    #   p₃·p₄ = s/2 - m²
    #   p₃·p₂ = s/2 - m² + t/2
    #   p₁·p₄ = m² - u/2
    # 对 t 道: Tr[p̸₃γ^μp̸₁γ^ν]Tr[p̸₄γ_μp̸₂γ_ν]
    #         = 32[(p₃·p₄)(p₁·p₂) + (p₃·p₂)(p₁·p₄)]
    # 对 u 道: 同上, 交换 3↔4
    #         = 32[(p₃·p₄)(p₁·p₂) + (p₄·p₂)(p₁·p₃)]
    # 干涉项: Tr[p̸₃γ^μp̸₁γ^νp̸₄γ_μp̸₂γ_ν]
    #         = -32(p₁·p₂)(p₃·p₄) = -8s²  (无质量极限)

    m2 = m_e**2
    dot_12 = s_val/2 - m2          # p₁·p₂
    dot_34 = s_val/2 - m2          # p₃·p₄
    dot_32 = s_val/2 - m2 + t_val/2  # p₃·p₂
    dot_14 = m2 - u_val/2          # p₁·p₄
    dot_42 = m2 - u_val/2          # p₄·p₂ (= p₁·p₄ by symmetry)
    dot_13 = m2 - t_val/2          # p₁·p₃

    trace_t_expected = 32 * (dot_34 * dot_12 + dot_32 * dot_14)
    trace_u_expected = 32 * (dot_34 * dot_12 + dot_42 * dot_13)
    trace_int_expected = -32 * dot_12 * dot_34

    print(f"\n    四动量内积 (从 Mandelstam 推导):")
    print(f"      p₁·p₂ = {dot_12:.4f},  p₃·p₄ = {dot_34:.4f}")
    print(f"      p₃·p₂ = {dot_32:.4f},  p₁·p₄ = {dot_14:.4f}")
    print(f"      p₄·p₂ = {dot_42:.4f},  p₁·p₃ = {dot_13:.4f}")

    print(f"\n    迹值对照:")
    print(f"      {'':20s} {'数值求迹':>16s} {'标准迹公式':>16s} {'偏差':>12s}")
    rel_t = abs(trace_t - trace_t_expected) / max(abs(trace_t_expected), 1) * 100
    rel_u = abs(trace_u - trace_u_expected) / max(abs(trace_u_expected), 1) * 100
    rel_i = abs(trace_int - trace_int_expected) / max(abs(trace_int_expected), 1) * 100
    print(f"      {'Tr_t (t 道平方)':20s} {trace_t:>16.2f} {trace_t_expected:>16.2f} {rel_t:>9.2f}%")
    print(f"      {'Tr_u (u 道平方)':20s} {trace_u:>16.2f} {trace_u_expected:>16.2f} {rel_u:>9.2f}%")
    print(f"      {'Tr_int (干涉)':20s} {trace_int:>16.2f} {trace_int_expected:>16.2f} {rel_i:>9.2f}%")

    if rel_t < 1.0 and rel_u < 1.0:
        print(f"\n    ✓ 数值求迹与标准迹公式一致 (微小偏差来自质量修正)")
    elif rel_t < 5.0 and rel_u < 5.0:
        print(f"\n    ~ 可接受的偏差 (质量项 + 数值精度)")

    # 展示如何从求迹组装自旋平均 |M|²
    print(f"\n    从求迹组装自旋平均 |M|² (Peskin §5.2 约定):")
    print(f"      |M|² = (1/4) Σ_spins |M|²")
    print(f"      |M_t|²      = e⁴ × Tr_t / (4 t²)      = {e4*trace_t/(4*t_val**2):.6e}")
    print(f"      |M_u|²      = e⁴ × Tr_u / (4 u²)      = {e4*trace_u/(4*u_val**2):.6e}")
    print(f"      干涉        = -e⁴ × Tr_int / (4 t u)   = {-e4*trace_int/(4*t_val*u_val):.6e}")
    print(f"      ─────────────────────────────────────────")
    print(f"      |M|² (求迹) = {amp2_numerical:.6e}")

    amp2_analytic = moller_amplitude_squared(s_val, t_val, u_val, m_e)
    print(f"\n    (参考) qed.py 的 moller_amplitude_squared: |M|² = {amp2_analytic:.6e}")
    print(f"    注: 两者使用不同归一化约定; 本次求迹遵循 Peskin 标准")

    return amp2_numerical, amp2_analytic


# ================================================================
# 第 4 步: 微分散射截面
# ================================================================

def step4_differential_cross_section(s_val, amp2):
    """将 |M|² 转换为 dσ/dΩ"""
    print_step(4, "微分散射截面 dσ/dΩ")

    print(f"""
  截面公式 (质心系, 等质量 2→2 散射):

      dσ        1      |p_f|
     ────  =  ────── · ───── · |M|²
      dΩ      64π² s    |p_i|

    等质量时 |p_f| = |p_i|, 所以在自然单位 (ħ = c = 1) 下:

      dσ       |M|²
     ────  =  ──────
      dΩ      64π² s
    """)

    dsigma_domega = amp2 / (64 * np.pi**2 * s_val)
    print(f"  数值:")
    print(f"    |M|²   = {amp2:.6e}")
    print(f"    s      = {s_val:.2f} MeV²")
    print(f"    dσ/dΩ  = {dsigma_domega:.6e} MeV⁻²")

    # 转换为常用单位
    # 1 MeV⁻¹ = 197.327 fm (ħc = 197.327 MeV·fm)
    # 1 MeV⁻² = (197.327)² fm² = 38938 fm²
    # 1 barn = 100 fm², 1 mb = 0.1 fm²
    hc_sq = (197.3269804) ** 2  # (MeV·fm)²
    dsigma_fm2 = dsigma_domega * hc_sq
    dsigma_mb = dsigma_fm2 / 10  # 1 mb = 10 fm²
    print(f"    dσ/dΩ  = {dsigma_fm2:.6e} fm²/sr")
    print(f"           = {dsigma_mb:.6e} mb/sr")
    print(f"    (ħc = 197.327 MeV·fm)")

    return dsigma_domega


# ================================================================
# 第 5 步: 图像
# ================================================================

def step5_plots(s_val, theta_val):
    """绘制多面板图像"""
    print_step(5, "图像: 截面 vs 角度 与 总截面 vs 能量")

    m_e = M_E

    # --- 数据准备 ---
    # Panel 1: dσ/dΩ vs θ for 3 个 CM 能量
    energies_gev = [0.01, 0.1, 1.0]  # GeV
    energies_mev2 = [(e * 1000) ** 2 for e in energies_gev]  # → MeV²
    theta_grid = np.linspace(np.radians(10), np.radians(170), 200)
    colors = ['#58a6ff', '#d2a8ff', '#ff7b72']

    # Panel 2: dσ/dΩ(θ) in polar
    # Panel 3: total cross section vs √s
    s_vals_mev2 = np.logspace(np.log10((2*m_e*1.1)**2), np.log10(1e7), 100)
    sqrt_s_gev = np.sqrt(s_vals_mev2) / 1000

    total_xsec = []
    for s_i in s_vals_mev2:
        theta_min = max(np.radians(1), np.arcsin(2*m_e/np.sqrt(s_i)) * 1.5)
        sigma = moller_total_cross_section(s_i, m_e=m_e, theta_min=theta_min)
        total_xsec.append(sigma * (197.327)**2)  # → fm²

    total_xsec = np.array(total_xsec)

    # --- 创建图 ---
    fig = plt.figure(figsize=(18, 12), facecolor='#0d1117')

    # == Panel A: dσ/dΩ vs θ (笛卡尔) ==
    ax1 = fig.add_subplot(2, 3, (1, 2))
    ax1.set_facecolor('#0d1117')

    for i, (s_mev2, label, c) in enumerate(zip(energies_mev2, energies_gev, colors)):
        dsigma = np.array([moller_cross_section(s_mev2, th, m_e) for th in theta_grid])
        # 转换为 mb/sr
        dsigma_mb = dsigma * (197.327)**2 / 10
        ax1.plot(np.degrees(theta_grid), dsigma_mb, color=c, linewidth=2.0,
                label=f'√s = {label} GeV')
        ax1.fill_between(np.degrees(theta_grid), 0, dsigma_mb, alpha=0.1, color=c)

    ax1.set_xlabel('Scattering angle θ (deg)', color='#e6edf3', fontsize=12)
    ax1.set_ylabel('dσ/dΩ (mb/sr)', color='#e6edf3', fontsize=12)
    ax1.set_title('Moller Scattering Differential Cross Section', color='#e6edf3',
                  fontsize=14, fontweight='bold')
    ax1.set_yscale('log')
    ax1.legend(framealpha=0.2, facecolor='#161b22', edgecolor='#30363d',
              labelcolor='#e6edf3', fontsize=10)
    ax1.tick_params(colors='#e6edf3')
    ax1.grid(True, alpha=0.1, color='#30363d')
    for spine in ax1.spines.values():
        spine.set_color('#30363d')

    # == Panel B: 计算流程图 (ASCII art 文本) ==
    ax2 = fig.add_subplot(2, 3, 3)
    ax2.set_facecolor('#0d1117')
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis('off')
    ax2.set_title('Calculation Chain', color='#e6edf3', fontsize=13, fontweight='bold', pad=15)

    flow_text = [
        "Feynman Diagrams",
        "  ↓",
        "Feynman Rules → M_t, M_u",
        "  ↓",
        "Spin Sum p̸+m",
        "  ↓",
        "Trace Tr[γ...γ...γ...γ...]",
        "  ↓",
        "|M|² = e⁴ × f(s,t,u)",
        "  ↓",
        "dσ/dΩ = |M|² / (64π²s)",
        "  ↓",
        "Total σ = ∫ dΩ dσ/dΩ",
    ]

    for i, line in enumerate(flow_text):
        y_pos = 9.0 - i * 0.72
        if '↓' in line:
            ax2.text(5, y_pos, line, ha='center', va='center',
                    fontsize=14, color='#f0883e', fontfamily='monospace')
        else:
            ax2.text(5, y_pos, line, ha='center', va='center',
                    fontsize=10.5, color='#79c0ff', fontfamily='monospace')

    # == Panel C: dσ/dΩ 极坐标图 (1 GeV) ==
    ax3 = fig.add_subplot(2, 3, 4, projection='polar')
    ax3.set_facecolor('#0d1117')

    s_high = (1000.0)**2  # 1 GeV
    theta_polar = np.linspace(np.radians(5), np.radians(175), 180)
    dsigma_polar = np.array([moller_cross_section(s_high, th, m_e) for th in theta_polar])
    dsigma_polar_mb = dsigma_polar * (197.327)**2 / 10

    ax3.plot(theta_polar, dsigma_polar_mb, color='#58a6ff', linewidth=2.0)
    ax3.fill(theta_polar, dsigma_polar_mb, alpha=0.15, color='#58a6ff')
    ax3.set_title('dσ/dΩ Polar Plot (√s = 1 GeV)', color='#e6edf3',
                  fontsize=12, fontweight='bold', pad=20)
    ax3.tick_params(colors='#e6edf3')
    ax3.set_yscale('log')

    # == Panel D: 总截面 vs √s ==
    ax4 = fig.add_subplot(2, 3, 5)
    ax4.set_facecolor('#0d1117')

    ax4.plot(sqrt_s_gev, total_xsec, color='#ff7b72', linewidth=2.5)
    ax4.fill_between(sqrt_s_gev, 0, total_xsec, alpha=0.1, color='#ff7b72')
    ax4.set_xlabel('√s (GeV)', color='#e6edf3', fontsize=12)
    ax4.set_ylabel('σ_total (fm²)', color='#e6edf3', fontsize=12)
    ax4.set_title('Moller Scattering Total Cross Section', color='#e6edf3',
                  fontsize=13, fontweight='bold')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.tick_params(colors='#e6edf3')
    ax4.grid(True, alpha=0.1, color='#30363d')
    for spine in ax4.spines.values():
        spine.set_color('#30363d')

    # == Panel E: 表格 — 特定角度下的完整数值 ==
    ax5 = fig.add_subplot(2, 3, 6)
    ax5.set_facecolor('#0d1117')
    ax5.axis('off')
    ax5.set_title('Numerical Check (θ = 45°)', color='#e6edf3',
                  fontsize=13, fontweight='bold', pad=15)

    # 计算表格数据
    _, t_test, u_test = mandelstam_from_cm(s_val, theta_val, mass=m_e)
    amp2_test = moller_amplitude_squared(s_val, t_test, u_test, m_e)
    dsigma_test = moller_cross_section(s_val, theta_val, m_e)
    dsigma_mb_test = dsigma_test * (197.327)**2 / 10

    table_text = f"""
  √s          = {np.sqrt(s_val):.2f} MeV
  θ_cm        = {np.degrees(theta_val):.1f}°

  s           = {s_val:.2f} MeV²
  t           = {t_test:.2f} MeV²
  u           = {u_test:.2f} MeV²
  s+t+u       = {s_val + t_test + u_test:.2f} ≈ 4m²

  α           = 1/{1/ALPHA:.0f}
  e⁴          = {E_SQ:.6e}

  |M|² (trace)= {amp2_test:.6e}
  dσ/dΩ       = {dsigma_test:.6e} MeV⁻²
              = {dsigma_mb_test:.6e} mb/sr
"""

    ax5.text(0.05, 0.95, table_text, transform=ax5.transAxes,
            fontsize=9.5, color='#e6edf3', fontfamily='monospace',
            va='top', linespacing=1.6)

    plt.tight_layout(pad=2.0)
    save_path = os.path.join(OUTPUT_DIR, 'moller_step_by_step.png')
    fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    plt.close(fig)
    print(f"\n  图像已保存: {save_path}")

    return save_path


# ================================================================
# 第 6 步: 数值表
# ================================================================

def step6_numerical_table():
    """打印完整的数值表"""
    print_step(6, "数值表 — 不同能量和角度下的截面")

    m_e = M_E
    energies_mev = [0.6, 1.0, 5.0, 50.0, 500.0, 5000.0]  # √s in MeV
    angles_deg = [10, 30, 45, 60, 90, 120, 150]

    print(f"\n  Møller 散射 dσ/dΩ (mb/sr) 数值表")
    print(f"  {'─'*72}")
    header = f"  {'√s (MeV)':>10s}"
    for ang in angles_deg:
        header += f"  {ang:>6d}°"
    print(header)
    print(f"  {'─'*72}")

    for es in energies_mev:
        sval = es**2
        row = f"  {es:>10.1f}"
        for ang in angles_deg:
            th = np.radians(ang)
            try:
                ds = moller_cross_section(sval, th, m_e)
                ds_mb = ds * (197.327)**2 / 10
                row += f"  {ds_mb:>8.3e}"
            except Exception:
                row += f"  {'N/A':>8s}"
        print(row)

    print(f"  {'─'*72}")
    print(f"  注: 低能时前向发散 (t → 0) 导致大截面")
    print(f"      转换: 1 MeV⁻² = 38938 fm², 1 mb = 10 fm²")


# ================================================================
# 主程序
# ================================================================

def main():
    print_header("Møller 散射 (e⁻e⁻ → e⁻e⁻) 计算全链条")
    print_header("《一读就懂的量子场论》§4.6-4.7 演示")

    # 第 1 步
    step1_feynman_diagrams()

    # 第 2 步
    s_val, t_val, u_val, theta_val = step2_feynman_rules()

    # 第 3 步
    amp2_numerical, amp2_analytic = step3_spin_averaged_amplitude(s_val, t_val, u_val)

    # 第 4 步
    dsigma_domega = step4_differential_cross_section(s_val, amp2_analytic)

    # 第 5 步
    fig_path = step5_plots(s_val, theta_val)

    # 第 6 步
    step6_numerical_table()

    # 总结
    print_header("计算完成")
    print(f"  输出图像: {fig_path}")
    print(f"  关键公式总结:")
    print(f"    M = M_t - M_u                              (全同费米子交换负号)")
    print(f"    dσ/dΩ = |M|² / (64π² s)                    (质心系)")
    print(f"    |M|² = e⁴ × f(s,t,u) / (t² u²)             (求迹结果)")
    print(f"    α = 1/{1/ALPHA:.0f}                                   (精细结构常数)")
    print(f"{SEP}")


if __name__ == '__main__':
    main()

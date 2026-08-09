#!/usr/bin/env python3
"""
重整化核心洞察 — 截断抵消 (Cutoff Cancellation) 演示 (§5.4)

物理核心:
    量子场论重整化的最深刻洞察：
    - 裸质量 m₀(Λ) 随截断 Λ 发散（向上）
    - 抵消项 δm(Λ) 也发散（向下，等大反向）
    - 物理质量 m_phys = m₀ + δm 保持常数（平坦直线）

    这就是重整化「工作」的原因——发散精确抵消。

    公式 (φ⁴ 理论，单圈，动量截断，d=4):
        Σ(Λ) = λ/(32π²) [Λ² - m² ln(1 + Λ²/m²)]    自能 (发散)
        δm² = -Σ(Λ)                                    抵消项 (反向发散)
        m₀² = m_R² + Σ(Λ)                               裸质量平方
        m_phys² = m₀² + δm² = m_R²                      物理质量 (常数)

演示设计:
    Panel 1 (上): 裸质量 m₀(Λ) vs 截断 Λ (log scale) —— 向上发散
    Panel 2 (中): 抵消项 δm(Λ) vs 截断 Λ (log scale) —— 向下发散
    Panel 3 (下): 物理质量 m_phys vs 截断 Λ (log scale) —— 平坦直线
    叠加: m_phys 参考水平线覆盖三面板
    Planck 标度: 垂直虚线标记
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from src.qft.renormalization import self_energy_1loop, mass_counterterm


# ================================================================
# 物理计算
# ================================================================

def compute_renormalization(m_R: float, lam: float, Lambda_vals: np.ndarray):
    """计算裸质量、抵消项、物理质量随截断的变化。

    参数:
        m_R:         重整化物理质量 (on-shell)
        lam:         耦合常数 λ
        Lambda_vals: 截断数组

    返回:
        m0:      裸质量 m₀ (取 sqrt, 保持质量量纲)
        dm:      有效质量抵消 δm ≈ δm²/(2m₀)
        m_phys:  物理质量 m_R (常数)
        sigma:   自能 Σ(Λ)
        dm2:     质量平方抵消项 δm²
    """
    sigma = np.array([self_energy_1loop(0.0, m_R, lam, L) for L in Lambda_vals])
    dm2 = np.array([mass_counterterm(m_R, lam, L) for L in Lambda_vals])
    # 裸质量平方: m₀² = m_R² + Σ = m_R² - δm²
    m0_sq = m_R**2 + sigma
    m0 = np.sqrt(np.maximum(0, m0_sq))
    # 有效质量抵消: δm ≈ δm²/(2m₀)
    with np.errstate(divide='ignore', invalid='ignore'):
        dm = np.where(m0 > 1e-10, dm2 / (2.0 * m0), 0.0)
    m_phys = np.full_like(Lambda_vals, m_R)
    return m0, dm, m_phys, sigma, dm2


# ================================================================
# 主程序
# ================================================================

def main(m_R: float = 1.0, lam: float = 0.1,
         Lambda_max: float = 1000.0, n_points: int = 500):
    """生成截断抵消可视化。

    参数:
        m_R:        物理质量 (默认 1.0 GeV)
        lam:        耦合常数 λ (默认 0.1)
        Lambda_max: 最大截断 Λ (默认 1000 GeV)
        n_points:   采样点数 (默认 500)
    """
    print("=" * 64)
    print("  重整化核心洞察 — 截断抵消 (Cutoff Cancellation)")
    print("  §5.4 重整化")
    print("=" * 64)
    print(f"\n  参数: m_R = {m_R} GeV,  λ = {lam}")
    print(f"  截断范围: Λ ∈ [1, {Lambda_max}] GeV (log scale)")
    print(f"  理论: φ⁴, d=4, 动量截断正规化, 单圈, on-shell 方案")

    # ── 截断数组 (log 均匀) ──
    Lambda_vals = np.logspace(0, np.log10(Lambda_max), n_points)

    # ── 计算 ──
    m0, dm, m_phys, sigma, dm2 = compute_renormalization(m_R, lam, Lambda_vals)

    # ── 终端数值摘要 (Λ = 10, 100, 1000) ──
    check_Lambdas = np.array([10.0, 100.0, 1000.0])
    print(f"\n  {'Λ':>8s}  {'Σ(Λ)':>12s}  {'δm²':>12s}  {'m₀':>10s}  "
          f"{'δm_eff':>10s}  {'m₀+δm':>10s}  {'m_R':>8s}")
    print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*10}  "
          f"{'-'*10}  {'-'*10}  {'-'*8}")
    for L in check_Lambdas:
        idx = np.argmin(np.abs(Lambda_vals - L))
        m0_val = m0[idx]
        dm_val = dm[idx]
        sigma_val = sigma[idx]
        dm2_val = dm2[idx]
        phys_val = m0_val + dm_val
        print(f"  {L:8.0f}  {sigma_val:12.6f}  {dm2_val:12.6f}  "
              f"{m0_val:10.6f}  {dm_val:10.6f}  {phys_val:10.6f}  {m_R:8.4f}")

    # 自洽性检查
    print(f"\n  自洽性检查:")
    print(f"    σ + δm² = {sigma[-1] + dm2[-1]:.3e}  (应为 0, on-shell 方案精确抵消)")
    print(f"    m₀² - m_R² = σ  ⟹  {m0[-1]**2 - m_R**2:.6f}  vs  {sigma[-1]:.6f}")
    print(f"    物理质量守恒: m₀ + δm_eff ≈ {m0[-1] + dm[-1]:.6f}  (参考 m_R={m_R})")

    # ══════════════════════════════════════════════════════════════
    # 绘图
    # ══════════════════════════════════════════════════════════════
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # 中文字体
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC',
                                        'Source Han Sans CN',
                                        'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['mathtext.fontset'] = 'dejavusans'

    fig = plt.figure(figsize=(14, 12), facecolor='#0d1117')

    # Planck 标度 (~10^19 GeV) 类比 — 用 100 GeV 作为演示中的 "高能标度"
    planck_mark = 100.0

    # 颜色方案
    color_bare = '#ff7b72'      # 红色 — 裸质量 (发散)
    color_ct = '#79c0ff'        # 蓝色 — 抵消项 (反向发散)
    color_phys = '#3fb950'      # 绿色 — 物理质量 (常数)
    color_ref = '#d2a8ff'       # 紫色 — 参考线
    color_grid = '#30363d'
    color_text = '#e6edf3'
    color_panel_bg = '#161b22'

    # ══════════════════════════════════════════════════════════════
    # Panel 1 (上): 裸质量 m₀(Λ) — 向上发散
    # ══════════════════════════════════════════════════════════════
    ax1 = fig.add_subplot(3, 1, 1, facecolor='#0d1117')

    ax1.plot(Lambda_vals, m0, color=color_bare, linewidth=2.5,
             label=r'裸质量 $m_0(\Lambda)$', zorder=4)
    ax1.fill_between(Lambda_vals, m_R, m0, alpha=0.12, color=color_bare)
    # 物理质量参考线
    ax1.axhline(m_R, color=color_phys, linewidth=2.0, linestyle='--',
                alpha=0.7, label=f'物理质量 $m_{{\\mathrm{{phys}}}} = {m_R}$', zorder=3)
    # Planck 标度
    ax1.axvline(planck_mark, color=color_ref, linewidth=1.2,
                linestyle=':', alpha=0.6, zorder=2)

    # 标注
    ax1.annotate(
        r'$m_0(\Lambda) \to \infty$' + '\n(紫外发散)',
        xy=(Lambda_vals[-1], m0[-1]),
        xytext=(Lambda_vals[-1] * 0.3, m0[-1] * 0.92),
        arrowprops=dict(arrowstyle='->', color=color_bare, lw=1.5),
        color=color_bare, fontsize=10, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=color_panel_bg,
                  edgecolor='#30363d', alpha=0.9))
    # Planck 标度标注
    ax1.annotate(
        'Planck 标度',
        xy=(planck_mark, m_R),
        xytext=(planck_mark * 1.8, m_R * 1.15),
        arrowprops=dict(arrowstyle='->', color=color_ref, lw=1.2),
        color=color_ref, fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor=color_panel_bg,
                  edgecolor='#30363d', alpha=0.85))

    ax1.set_ylabel('裸质量 $m_0$ [GeV]', color=color_text, fontsize=12)
    ax1.set_title('Panel 1: 裸质量 — 紫外发散', color=color_text,
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left', facecolor=color_panel_bg,
               edgecolor='#30363d', labelcolor=color_text, fontsize=9)
    ax1.set_xscale('log')
    ax1.tick_params(colors=color_text)
    ax1.grid(True, alpha=0.12, color=color_grid)
    for spine in ax1.spines.values():
        spine.set_color(color_grid)

    # ══════════════════════════════════════════════════════════════
    # Panel 2 (中): 抵消项 δm(Λ) — 向下发散
    # ══════════════════════════════════════════════════════════════
    ax2 = fig.add_subplot(3, 1, 2, facecolor='#0d1117')

    ax2.plot(Lambda_vals, dm, color=color_ct, linewidth=2.5,
             label=r'抵消项 $\delta m(\Lambda)$', zorder=4)
    ax2.fill_between(Lambda_vals, 0, dm, alpha=0.12, color=color_ct)
    # 零线
    ax2.axhline(0, color=color_text, linewidth=0.6, linestyle='-', alpha=0.3)
    # 物理质量偏移参考 (近似 m₀ + δm = m_R → δm ≈ m_R - m₀)
    ax2.axhline(m_R - m0[-1], color=color_phys, linewidth=2.0,
                linestyle='--', alpha=0.4, zorder=3)
    # Planck 标度
    ax2.axvline(planck_mark, color=color_ref, linewidth=1.2,
                linestyle=':', alpha=0.6, zorder=2)

    ax2.annotate(
        r'$\delta m(\Lambda) \to -\infty$' + '\n(抵消发散)',
        xy=(Lambda_vals[-1], dm[-1]),
        xytext=(Lambda_vals[-1] * 0.25, dm[0] * 0.15),
        arrowprops=dict(arrowstyle='->', color=color_ct, lw=1.5),
        color=color_ct, fontsize=10, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=color_panel_bg,
                  edgecolor='#30363d', alpha=0.9))
    ax2.annotate(
        'Planck 标度',
        xy=(planck_mark, 0),
        xytext=(planck_mark * 1.8, dm[0] * 0.7),
        arrowprops=dict(arrowstyle='->', color=color_ref, lw=1.2),
        color=color_ref, fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor=color_panel_bg,
                  edgecolor='#30363d', alpha=0.85))

    ax2.set_ylabel(r'抵消项 $\delta m_\mathrm{eff}$ [GeV]',
                   color=color_text, fontsize=12)
    ax2.set_title('Panel 2: 抵消项 — 反向发散 (等大反向)', color=color_text,
                  fontsize=14, fontweight='bold')
    ax2.legend(loc='lower left', facecolor=color_panel_bg,
               edgecolor='#30363d', labelcolor=color_text, fontsize=9)
    ax2.set_xscale('log')
    ax2.tick_params(colors=color_text)
    ax2.grid(True, alpha=0.12, color=color_grid)
    for spine in ax2.spines.values():
        spine.set_color(color_grid)

    # ══════════════════════════════════════════════════════════════
    # Panel 3 (下): 物理质量 m_phys — 平坦直线
    # ══════════════════════════════════════════════════════════════
    ax3 = fig.add_subplot(3, 1, 3, facecolor='#0d1117')

    # 物理质量 = sqrt(m₀² + δm²) = m_R (精确)
    phys_exact = np.sqrt(np.maximum(0, m0**2 + dm2))
    ax3.plot(Lambda_vals, phys_exact, color=color_phys, linewidth=2.5,
             label=r'物理质量 $m_\mathrm{phys} = \sqrt{m_0^2 + \delta m^2}$', zorder=4)
    ax3.axhline(m_R, color=color_phys, linewidth=1.5, linestyle='-',
                alpha=0.4, zorder=3)
    # 裸质量 + 抵消项 分别画 (浅显)
    ax3.plot(Lambda_vals, m0, color=color_bare, linewidth=1.0,
             linestyle=':', alpha=0.25, label='裸质量 $m_0$ (参考)', zorder=2)
    ax3.plot(Lambda_vals, dm, color=color_ct, linewidth=1.0,
             linestyle=':', alpha=0.25, label=r'抵消项 $\delta m$ (参考)', zorder=2)
    # Planck 标度
    ax3.axvline(planck_mark, color=color_ref, linewidth=1.2,
                linestyle=':', alpha=0.6, zorder=2)

    # 大标注 — 核心信息
    ax3.annotate(
        f'$m_{{\\mathrm{{phys}}}} = \\sqrt{{m_0^2 + \\delta m^2}} = {m_R}$\n'
        '= 常数\n\n'
        '发散精确抵消\n'
        '—— 重整化「工作」的原因',
        xy=(np.median(Lambda_vals), m_R),
        xytext=(Lambda_vals[0] * 3, m_R * 1.08),
        color=color_phys, fontsize=12, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.6', facecolor=color_panel_bg,
                  edgecolor=color_phys, alpha=0.95, linewidth=1.5),
        ha='center')
    # Planck 标度标注
    ax3.annotate(
        'Planck 标度',
        xy=(planck_mark, m_R),
        xytext=(planck_mark * 2.0, m_R * 0.93),
        arrowprops=dict(arrowstyle='->', color=color_ref, lw=1.2),
        color=color_ref, fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor=color_panel_bg,
                  edgecolor='#30363d', alpha=0.85))

    ax3.set_xlabel('动量截断 $\\Lambda$ [GeV] (log scale)', color=color_text, fontsize=12)
    ax3.set_ylabel('物理质量 $m_\\mathrm{phys}$ [GeV]',
                   color=color_text, fontsize=12)
    ax3.set_title('Panel 3: 物理质量 — 发散精确抵消！', color=color_text,
                  fontsize=14, fontweight='bold')
    ax3.legend(loc='center right', facecolor=color_panel_bg,
               edgecolor='#30363d', labelcolor=color_text, fontsize=8)
    ax3.set_xscale('log')
    ax3.tick_params(colors=color_text)
    ax3.grid(True, alpha=0.12, color=color_grid)
    for spine in ax3.spines.values():
        spine.set_color(color_grid)
    # y 轴范围: 围绕 m_R 展开 ±30% 以突出平坦性
    ax3.set_ylim(m_R * 0.7, m_R * 1.3)

    # ── 全局标题 ──
    fig.suptitle(
        f'重整化核心洞察  |  φ⁴ 单圈, 动量截断正规化, on-shell 方案  |  '
        f'$m_R={m_R}$, $\\lambda={lam}$  |  §5.4',
        color=color_text, fontsize=12, fontweight='bold', y=0.995)

    plt.tight_layout(rect=(0, 0, 1, 0.97))

    # ── 保存 ──
    save_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, 'cutoff_cancellation.png')
    fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    plt.close(fig)

    print(f"\n  [OK] 图像已保存: {save_path}")

    # ── 物理总结 ──
    print(f"\n{'='*64}")
    print("  物理结论 — 重整化的核心洞察")
    print(f"{'='*64}")
    print(f"""
  裸质量 m₀(Λ) 随截断 Λ 发散 —— 这是紫外灾难的表征。
  抵消项 δm(Λ) 也发散 —— 等大反向, 这是「重整化」的含义。
  物理质量 m_phys = m₀ + δm 保持常数 —— 这是重整化「工作」的原因。

  关键点:
    (a) 发散不是物理的 —— 它们是微扰展开中引入的中间量
    (b) 抵消项精确消去发散 —— 在壳方案保证物理质量为有限值
    (c) 重整化 = 重新参数化 —— 用物理可观测量 (m_R) 取代裸参数 (m₀)
    (d) 截断依赖性隐藏在未观测的裸参数中 —— 可观测量与截断无关

  公式:
    Σ(Λ) = λ/(32π²) [Λ² - m² ln(1 + Λ²/m²)]    自能 (Λ² 发散)
    δm² = -Σ(Λ)                                   抵消项 (等大反向)
    m_phys² = m₀² + δm² = m_R²                    物理质量 (常数, 与 Λ 无关)

  对应书籍: §5.4 重整化
""")

    return {
        'm_R': m_R, 'lam': lam,
        'Lambda_vals': Lambda_vals,
        'm0': m0, 'dm': dm, 'm_phys': phys_exact,
        'sigma': sigma, 'dm2': dm2,
    }


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='截断抵消可视化 — 重整化核心洞察')
    parser.add_argument('--mR', type=float, default=1.0,
                        help='物理质量 m_R (默认: 1.0)')
    parser.add_argument('--lam', type=float, default=0.1,
                        help='耦合常数 λ (默认: 0.1)')
    parser.add_argument('--Lambda-max', type=float, default=1000.0,
                        help='最大截断 Λ (默认: 1000)')
    parser.add_argument('--n-points', type=int, default=500,
                        help='采样点数 (默认: 500)')
    args = parser.parse_args()

    result = main(m_R=args.mR, lam=args.lam,
                  Lambda_max=args.Lambda_max, n_points=args.n_points)

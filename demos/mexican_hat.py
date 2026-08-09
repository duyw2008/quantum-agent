#!/usr/bin/env python3
"""
Mexican Hat 势可视化 — 自发对称性破缺的经典图示 (§6.4-6.5)

展示:
  1. 2D 剖面: V(φ) 随 m² 从正到负的相变过程
     - m² > 0:  单阱 (对称相, 极小值在 φ=0)
     - m² = 0:  平坦底 (临界点)
     - m² < 0:  双阱 (破缺相, 极小值在 ±√(-6m²/λ))
  2. 3D 曲面: 复标量场 V(φ₁, φ₂) 的 Mexican hat 形状
     - 径向激发 = Higgs 玻色子 (有质量)
     - 角向激发 = Goldstone 玻色子 (无质量)
     - 谷底圆环 = 无穷简并真空

物理:
  V(φ) = ½m²|φ|² + (λ/4!)|φ|⁴

  当 m² < 0 时, φ=0 是局域极大值, 真实真空在 |φ| = v = √(-6m²/λ)
  展开 φ = (v + h) e^{iθ/v}:
    - h (径向/ Higgs) 质量 m_h = √(-2m²) = √(λv²/3)
    - θ (角向/ Goldstone) 质量 m_G = 0 (连续对称性破缺的 Nambu-Goldstone 定理)

输出: output/mexican_hat.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import FancyBboxPatch
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FuncFormatter

# 尝试导入项目模块 (可选: 用于一致性验证)
try:
    from src.qft.effective_potential import OneLoopEffectivePotential, find_minimum
    HAS_EFFECTIVE_POT = True
except ImportError:
    HAS_EFFECTIVE_POT = False

# ================================================================
# 全局样式
# ================================================================

BG = '#0d1117'
FG = '#e6edf3'
GRID = '#21262d'
ACCENT_BLUE = '#58a6ff'
ACCENT_ORANGE = '#f0883e'
ACCENT_GREEN = '#3fb950'
ACCENT_RED = '#f85149'
ACCENT_YELLOW = '#d29922'
ACCENT_PURPLE = '#bc8cff'

plt.rcParams.update({
    'figure.facecolor': BG,
    'axes.facecolor': BG,
    'axes.edgecolor': GRID,
    'axes.labelcolor': FG,
    'axes.titlecolor': FG,
    'xtick.color': FG,
    'ytick.color': FG,
    'text.color': FG,
    'grid.color': GRID,
    'grid.alpha': 0.4,
    'legend.facecolor': BG,
    'legend.edgecolor': GRID,
    'legend.labelcolor': FG,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Source Han Sans CN', 'Noto Sans CJK SC', 'DejaVu Sans'],
    'font.size': 11,
})

# ================================================================
# 势函数
# ================================================================

def tree_potential_real(phi, mass_sq, coupling):
    """实标量场树图势: V(φ) = ½m²φ² + (λ/4!)φ⁴"""
    return 0.5 * mass_sq * phi**2 + (coupling / 24.0) * phi**4


def tree_potential_complex(phi1, phi2, mass_sq, coupling):
    """复标量场树图势: V(φ) = ½m²|φ|² + (λ/4!)|φ|⁴"""
    r2 = phi1**2 + phi2**2
    return 0.5 * mass_sq * r2 + (coupling / 24.0) * r2**2


def minima_positions(mass_sq, coupling):
    """计算树图势极小值位置 φ = ±√(-6m²/λ) (仅 m²<0 时有效)"""
    if mass_sq >= 0:
        return [0.0]
    v = np.sqrt(-6.0 * mass_sq / coupling)
    return [-v, 0.0, v]  # 包含 φ=0 (局域极大值)


# ================================================================
# 主函数
# ================================================================

def main():
    coupling = 0.5   # λ (2D 剖面使用的耦合常数)

    # m²: 对称 → 临界 → 破缺 → 深破缺
    mass_sq_list = [+2.0, 0.0, -2.0, -5.0]
    m2_3d = -2.0       # 3D Mexican hat: 破缺相
    lam_3d = 1.0       # 3D Mexican hat: λ=1

    # --- 验证 (利用 effective_potential 模块) ---
    if HAS_EFFECTIVE_POT:
        pot_check = OneLoopEffectivePotential(mass_sq=-1.0, coupling=coupling)
        v_check = np.sqrt(-6.0 * (-1.0) / coupling)
        print(f"  [模块验证] OneLoopEffectivePotential 可用")
        print(f"    m²=-1, λ={coupling}: φ_min (解析) = ±{v_check:.4f}")
        phi_min, V_min = find_minimum(
            pot_check.tree_potential, (-v_check*2, v_check*2))
        print(f"    φ_min (数值搜索) = ±{phi_min:.4f}")
        print(f"    V_min = {V_min:.4f}")

    # ================================================================
    # 创建图形: 左 2D 剖面 + 右 3D Mexican Hat
    # ================================================================
    fig = plt.figure(figsize=(20, 9), facecolor=BG)

    # --- 左侧: 2D 剖面 (4条曲线) ---
    ax_2d = fig.add_axes((0.05, 0.12, 0.42, 0.82))
    ax_2d.set_facecolor(BG)

    phi_range = 6.0
    phi_vals = np.linspace(-phi_range, phi_range, 800)

    # 颜色梯度 (从对称到高度破缺)
    colors_2d = ['#58a6ff', '#3fb950', '#f0883e', '#f85149']
    labels_2d = [
        r'$m^2 = +2$ (对称相)',
        r'$m^2 = 0$ (临界点)',
        r'$m^2 = -2$ (破缺相)',
        r'$m^2 = -5$ (深破缺)',
    ]

    for m2, color, label in zip(mass_sq_list, colors_2d, labels_2d):
        V_vals = tree_potential_real(phi_vals, m2, coupling)
        # 将底部抬高以分离曲线 (便于观察形状)
        shift = 0.0
        ax_2d.plot(phi_vals, V_vals + shift, color=color, linewidth=2.2, label=label, zorder=3)

        # 标记极小值
        minima = minima_positions(m2, coupling)
        for m in minima:
            Vm = tree_potential_real(m, m2, coupling) + shift
            if m2 < 0 and abs(m) < 0.01:
                # φ=0 在破缺相是局域极大值
                ax_2d.plot(m, Vm, 'o', color=color, markersize=8,
                          markeredgecolor='white', markeredgewidth=1.2,
                          markerfacecolor='none', zorder=5)
            else:
                ax_2d.plot(m, Vm, 'o', color=color, markersize=10,
                          markeredgecolor='white', markeredgewidth=1.5,
                          zorder=5)

    # 标注破缺相的极小值和 Goldstone/Higgs
    m2_broken = -2.0
    v_val = np.sqrt(-6.0 * m2_broken / coupling)
    V_at_v = tree_potential_real(v_val, m2_broken, coupling)

    # 在 m²=-2 曲线上标注物理
    ax_2d.annotate(
        r'$\phi = \pm v = \pm\sqrt{-6m^2/\lambda}$',
        xy=(v_val, V_at_v), xytext=(v_val + 0.5, V_at_v + 3),
        fontsize=10, color=ACCENT_ORANGE,
        arrowprops=dict(arrowstyle='->', color=ACCENT_ORANGE, lw=1.5,
                       connectionstyle='arc3,rad=0.15'),
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG, edgecolor=ACCENT_ORANGE, alpha=0.85),
        zorder=10)

    ax_2d.annotate(
        r'$\phi = -\sqrt{-6m^2/\lambda}$',
        xy=(-v_val, V_at_v), xytext=(-v_val - 2.2, V_at_v + 2.5),
        fontsize=10, color=ACCENT_ORANGE,
        arrowprops=dict(arrowstyle='->', color=ACCENT_ORANGE, lw=1.5,
                       connectionstyle='arc3,rad=-0.15'),
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG, edgecolor=ACCENT_ORANGE, alpha=0.85),
        zorder=10)

    # Higgs 质量标注
    m_h_sq = -2.0 * m2_broken
    ax_2d.annotate(
        r'径向曲率 → $m_h^2 = -2m^2 = {:.1f}$'.format(m_h_sq),
        xy=(v_val + 0.3, V_at_v + 1.0), xytext=(v_val + 1.0, V_at_v + 6),
        fontsize=9.5, color=ACCENT_GREEN,
        arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=1.2),
        zorder=9)

    # φ=0 标注 (局域极大)
    V0_broken = tree_potential_real(0, m2_broken, coupling)
    ax_2d.annotate(
        r'$\phi = 0$ (不稳定平衡)',
        xy=(0, V0_broken), xytext=(0.8, V0_broken + 1.5),
        fontsize=9.5, color=ACCENT_RED,
        arrowprops=dict(arrowstyle='->', color=ACCENT_RED, lw=1.2),
        zorder=9)

    ax_2d.axhline(0, color=GRID, linewidth=0.8, linestyle='--', zorder=1)
    ax_2d.axvline(0, color=GRID, linewidth=0.8, linestyle='--', zorder=1)

    ax_2d.set_xlabel(r'$\phi$ (标量场)', fontsize=13, color=FG)
    ax_2d.set_ylabel(r'$V(\phi)$ (势能)', fontsize=13, color=FG)
    ax_2d.set_title(
        r'自发对称性破缺: $V(\phi)=\frac{1}{2}m^2\phi^2+\frac{\lambda}{4!}\phi^4$  '
        + r'($\lambda={:.1f}$)'.format(coupling),
        fontsize=14, fontweight='bold', color=FG, pad=12)

    ax_2d.set_xlim(-phi_range, phi_range)
    ax_2d.set_ylim(-10, 40)
    ax_2d.grid(True, alpha=0.25, color=GRID)
    ax_2d.legend(loc='upper center', fontsize=10, framealpha=0.9,
                 ncol=2, columnspacing=0.8)

    for spine in ax_2d.spines.values():
        spine.set_color(GRID)

    # --- 右侧: 3D Mexican Hat ---
    ax_3d = fig.add_axes((0.53, 0.05, 0.46, 0.92), projection='3d')
    ax_3d.set_facecolor(BG)
    # 使 3D 背景面板透明
    ax_3d.xaxis.pane.fill = False
    ax_3d.yaxis.pane.fill = False
    ax_3d.zaxis.pane.fill = False
    ax_3d.xaxis.pane.set_edgecolor(GRID)
    ax_3d.yaxis.pane.set_edgecolor(GRID)
    ax_3d.zaxis.pane.set_edgecolor(GRID)

    # 3D 曲面的参数 (已在上方定义: m2_3d, lam_3d)
    v_3d = np.sqrt(-6.0 * m2_3d / lam_3d)

    # 创建网格
    grid_size = 120
    r_max = v_3d * 2.2
    phi1 = np.linspace(-r_max, r_max, grid_size)
    phi2 = np.linspace(-r_max, r_max, grid_size)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)

    V_surf = tree_potential_complex(PHI1, PHI2, m2_3d, lam_3d)

    # 截断太大值以便可视化
    V_max_plot = V_surf.max() * 0.6
    V_surf_clipped = np.clip(V_surf, V_surf.min(), V_max_plot)

    # 自定义 colormap — 深底到亮顶
    colors_surf = [
        (0.0, '#0d1117'),     # 底部 (势阱) — 深色
        (0.15, '#1a3a5c'),    # 谷底
        (0.35, '#215070'),    # 上升
        (0.55, '#58a6ff'),    # 中部
        (0.75, '#bc8cff'),    # 高处
        (0.9, '#f0883e'),     # 顶部
        (1.0, '#f85149'),     # 顶点
    ]
    from matplotlib.colors import LinearSegmentedColormap
    mexican_cmap = LinearSegmentedColormap.from_list('mexican_hat', colors_surf)

    # 绘制曲面
    surf = ax_3d.plot_surface(
        PHI1, PHI2, V_surf_clipped,
        cmap=mexican_cmap,
        linewidth=0,
        antialiased=True,
        alpha=0.92,
        shade=True,
        lightsource=plt.matplotlib.colors.LightSource(azdeg=315, altdeg=45),
        zorder=1,
    )

    # 绘制谷底圆环 (极小值环)
    theta_circle = np.linspace(0, 2*np.pi, 200)
    circle_x = v_3d * np.cos(theta_circle)
    circle_y = v_3d * np.sin(theta_circle)
    circle_z = np.full_like(theta_circle, V_surf.min())
    ax_3d.plot(circle_x, circle_y, circle_z,
              color='#f0883e', linewidth=2.8, zorder=10, alpha=0.9,
              label=r'真空流形: $|\phi|=v$')

    # 内圈 (更小半径) — 显示势阱的曲率
    r_inner = v_3d * 0.5
    circle_x_in = r_inner * np.cos(theta_circle)
    circle_y_in = r_inner * np.sin(theta_circle)
    V_inner = tree_potential_complex(circle_x_in, circle_y_in, m2_3d, lam_3d)
    ax_3d.plot(circle_x_in, circle_y_in, V_inner,
              color='#3fb950', linewidth=1.5, zorder=10, alpha=0.6,
              linestyle='--')

    # 原点标记 (鞍点/局域极大)
    V_origin = tree_potential_complex(0, 0, m2_3d, lam_3d)
    ax_3d.scatter([0], [0], [V_origin], color='#f85149', s=80, zorder=11,
                 edgecolors='white', linewidths=1.5)

    # === 关键标注 ===

    # 径向模式 = Higgs 玻色子
    r_arrow_start = 0.3
    r_arrow_end = v_3d * 1.6
    ax_3d.annotate(
        '', xy=(r_arrow_end, 0), xytext=(r_arrow_start, 0),
        arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=3,
                       connectionstyle='arc3,rad=0'),
        zorder=20)
    # 标注文本 (使用 3D 坐标投影标注)
    ax_3d.text(
        v_3d * 1.25, -0.3, V_origin * 0.05,
        '径向模式\n= Higgs 玻色子\n(有质量)',
        color=ACCENT_GREEN, fontsize=12, fontweight='bold',
        ha='center', va='bottom', zorder=20,
        bbox=dict(boxstyle='round,pad=0.4', facecolor=BG,
                  edgecolor=ACCENT_GREEN, alpha=0.85))

    # 角向模式 = Goldstone 玻色子
    theta_arrow = np.pi / 4
    arc_r = v_3d * 1.05
    arc_theta = np.linspace(0.2, np.pi/2 - 0.1, 40)
    arc_x = arc_r * np.cos(arc_theta)
    arc_y = arc_r * np.sin(arc_theta)
    arc_z = np.full_like(arc_theta, V_surf.min())
    ax_3d.plot(arc_x, arc_y, arc_z,
              color=ACCENT_YELLOW, linewidth=3.5, zorder=15)
    # 箭头尖端
    ax_3d.annotate(
        '', xy=(arc_x[-1], arc_y[-1]),
        xytext=(arc_x[-2], arc_y[-2]),
        arrowprops=dict(arrowstyle='->', color=ACCENT_YELLOW, lw=3),
        zorder=20)

    ax_3d.text(
        v_3d * 1.0, v_3d * 1.15, V_surf.min() - 0.5,
        '角向模式\n= Goldstone 玻色子\n(无质量)',
        color=ACCENT_YELLOW, fontsize=12, fontweight='bold',
        ha='center', va='top', zorder=20,
        bbox=dict(boxstyle='round,pad=0.4', facecolor=BG,
                  edgecolor=ACCENT_YELLOW, alpha=0.85))

    # 原点标注
    ax_3d.text(
        0, -v_3d*0.25, V_origin + V_origin*0.08,
        r'$\phi=0$ (不稳定)',
        color=ACCENT_RED, fontsize=10, fontweight='bold',
        ha='center', zorder=20)

    # 坐标轴设置
    ax_3d.set_xlabel(r'$\phi_1$ (Re $\phi$)', fontsize=12, color=FG, labelpad=8)
    ax_3d.set_ylabel(r'$\phi_2$ (Im $\phi$)', fontsize=12, color=FG, labelpad=8)
    ax_3d.set_zlabel(r'$V(\phi_1,\phi_2)$', fontsize=12, color=FG, labelpad=8)
    ax_3d.set_title(
        '墨西哥帽势: ' + r'$m^2={:.0f} < 0,\; \lambda={:.1f}$'.format(m2_3d, lam_3d),
        fontsize=14, fontweight='bold', color=FG, pad=20)

    # 调整 3D 视图
    ax_3d.view_init(elev=28, azim=-55)

    # 设置 3D 轴颜色
    ax_3d.xaxis.line.set_color(GRID)
    ax_3d.yaxis.line.set_color(GRID)
    ax_3d.zaxis.line.set_color(GRID)
    ax_3d.tick_params(colors=FG, labelsize=9)

    # 限制 z 轴范围
    ax_3d.set_zlim(V_surf.min() - 1.0, V_max_plot)

    # --- 全局标题与说明 ---
    fig.suptitle(
        '墨西哥帽势 — 自发对称性破缺 (§6.4-6.5)',
        fontsize=17, fontweight='bold', color=FG, y=0.985)

    # 底部说明文字
    fig.text(
        0.5, 0.008,
        r'$V(\phi)=\frac{1}{2}m^2|\phi|^2+\frac{\lambda}{4!}|\phi|^4$'
        r'  |  $m^2<0$: 对称性自发破缺, '
        r'真空期望值 $|\langle\phi\rangle|=v=\sqrt{-6m^2/\lambda}$'
        r'  |  径向激发 = Higgs (有质量 $m_h=\sqrt{-2m^2}$), '
        r'角向激发 = Goldstone (无质量)',
        fontsize=9.5, color='#8b949e', ha='center', fontstyle='italic')

    # ================================================================
    # 保存
    # ================================================================
    save_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, 'mexican_hat.png')

    fig.savefig(save_path, dpi=200, bbox_inches='tight',
                facecolor=BG, edgecolor='none')
    plt.close(fig)

    print(f"\n  ✓ 墨西哥帽可视化已保存: {save_path}")
    print(f"    2D 剖面 (4 条曲线): m² = +2, 0, -2, -5")
    print(f"    3D 曲面: m² = {m2_3d}, λ = {lam_3d}")
    print(f"    真空期望值 v = {v_3d:.3f}")
    print(f"    Higgs 质量 m_h = {np.sqrt(-2*m2_3d):.3f}")
    print(f"    Goldstone 质量 m_G = 0")

    return save_path


if __name__ == '__main__':
    main()

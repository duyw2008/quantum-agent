#!/usr/bin/env python3
"""
费曼路径积分可视化 — 一维谐振子 "对所有路径求和"
==================================================

物理核心:
    量子传播子 K(x_f, t_f; x_i, t_i) 来自于对 ALL possible paths
    加权 exp(iS[x]/ħ) 的求和。在经典极限 ħ→0, 只有使作用量
    取极值 (δS=0) 的经典路径幸存 — 其它路径因相位干涉相消。

    K(x_f, x_i; T) = ∫ D[x(t)] exp(iS[x]/ħ)

演示设计:
    Panel 1 (左): 30 条随机路径 (细彩线) + 经典路径 (粗白线)
    Panel 2 (右): 所有 200 条路径的相位 exp(iS/ħ) 在复平面单位圆上的分布
        — 靠近经典作用量 S_cl 的路径相位聚集在一起 (稳定相位!)
        — 远离经典的路径相位随机散布, 相互抵消

对应书籍 §7.1-7.2
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb


# ================================================================
# 物理参数
# ================================================================

MASS = 1.0          # 粒子质量 m
OMEGA = 1.0         # 谐振子频率 ω
HBAR = 1.0          # 约化普朗克常数 ħ
XI = 0.8            # 初位置 x_i
XF = -0.5           # 末位置 x_f
TI = 0.0            # 初时间
TF = 2.0            # 末时间 T = tf - ti (避免 T=nπ 的共轭点)
N_T = 120           # 时间离散点数
N_PATHS = 200       # 总路径数
N_SHOW = 30         # 在 Panel 1 中显示的路径数
N_MODES = 15        # 随机路径的 Fourier 模数
FOURIER_SIGMA = 0.3  # Fourier 系数标准差 (控制量子涨落幅度, 带 1/n 衰减)


def classical_path(t, xi, xf, ti, tf, omega):
    """谐振子经典路径: 满足 ẍ + ω²x = 0, x(ti)=xi, x(tf)=xf

    x_cl(t) = A cos(ωt) + B sin(ωt)
    系数由边界条件确定。
    """
    T = tf - ti
    sin_wT = np.sin(omega * T)
    cos_wT = np.cos(omega * T)

    # 平移时间原点至 ti 简化计算
    tau = t - ti
    # x(τ) = A cos(ωτ) + B sin(ωτ)
    # x(0) = A = xi
    # x(T) = xi·cos(ωT) + B·sin(ωT) = xf  →  B = (xf - xi·cos(ωT)) / sin(ωT)
    A = xi
    B = (xf - xi * cos_wT) / sin_wT
    return A * np.cos(omega * tau) + B * np.sin(omega * tau)


def classical_action(xi, xf, ti, tf, mass, omega):
    """谐振子经典作用量 (解析公式)

    S_cl = (mω / 2sin(ωT)) [(xi² + xf²)cos(ωT) - 2xi xf]
    """
    T = tf - ti
    sin_wT = np.sin(omega * T)
    cos_wT = np.cos(omega * T)
    return (mass * omega / (2.0 * sin_wT)) * (
        (xi**2 + xf**2) * cos_wT - 2.0 * xi * xf
    )


# ================================================================
# 随机路径生成
# ================================================================

def generate_random_path(t_grid, xi, xf, ti, tf, n_modes, sigma):
    """生成一条满足端点条件的随机路径。

    使用 Fourier 正弦级数叠加于经典路径之上:
        x(t) = x_cl(t) + Σ_{n=1}^{N_modes} c_n · sin(nπ(t-ti)/T)

    正弦基自动满足在端点处为零, 因此边界条件 x(ti)=xi, x(tf)=xf 自然满足。
    """
    T = tf - ti
    tau = t_grid - ti  # [0, T]

    # 经典路径
    x_cl = classical_path(t_grid, xi, xf, ti, tf, OMEGA)

    # 随机 Fourier 系数 (高斯分布, 带 1/n 衰减以抑制高频动能)
    coeffs = np.random.randn(n_modes) * sigma
    # 按 1/n 衰减系数: 高频模对应更大的动能, 物理上应被抑制
    n_vals = np.arange(1, n_modes + 1)
    coeffs = coeffs / n_vals

    # 叠加正弦模
    fourier_sum = np.zeros_like(t_grid)
    for n in range(1, n_modes + 1):
        fourier_sum += coeffs[n - 1] * np.sin(n * np.pi * tau / T)

    return x_cl + fourier_sum, x_cl


# ================================================================
# 作用量计算
# ================================================================

def compute_action(path, t_grid, mass, hbar, omega):
    """离散化计算经典作用量 S[x]

    S[x] = ∫_{ti}^{tf} L(x, ẋ) dt
         = ∫ [½m ẋ² - ½mω²x²] dt

    离散化 (中点法则):
        S ≈ Σ_j [½m ((x_{j+1} - x_j)/Δt)² - ½mω² x_j²] · Δt
    """
    dt = t_grid[1] - t_grid[0]
    dx = path[1:] - path[:-1]
    x_mid = 0.5 * (path[1:] + path[:-1])

    kinetic = 0.5 * mass * (dx / dt)**2
    potential = 0.5 * mass * omega**2 * x_mid**2
    lagrangian = kinetic - potential

    return np.sum(lagrangian) * dt


# ================================================================
# 路径积分传播子 (离散化近似)
# ================================================================

def discretized_propagator(paths, actions, mass, hbar, dt):
    """离散路径积分的归一化传播子近似。

    K ≈ (m / 2πiħΔt)^{N/2} · (1/N_paths) · Σ exp(iS/ħ)

    此处仅计算相对值, 不做严格归一化。
    """
    phases = np.exp(1j * actions / hbar)
    return np.mean(phases)


# ================================================================
# 主程序
# ================================================================

def main():
    print("=" * 62)
    print("  费曼路径积分可视化 — 一维谐振子")
    print("  对所有路径求和 → 经典路径 ± 量子涨落")
    print("=" * 62)

    # ── 时间网格 ──
    t_grid = np.linspace(TI, TF, N_T)
    dt = t_grid[1] - t_grid[0]
    T = TF - TI

    print(f"\n  参数: m={MASS}, ω={OMEGA}, ħ={HBAR}")
    print(f"  端点: x_i={XI}, x_f={XF}")
    print(f"  时间: t_i={TI}, t_f={TF:.4f} (T={T:.4f})")
    print(f"  路径数: {N_PATHS}, 时间步: {N_T}, Fourier 模: {N_MODES}")

    # ── 经典路径 ──
    x_cl = classical_path(t_grid, XI, XF, TI, TF, OMEGA)
    S_cl = classical_action(XI, XF, TI, TF, MASS, OMEGA)
    phase_cl = np.exp(1j * S_cl / HBAR)

    print(f"\n  == 经典路径 ==")
    print(f"  S_cl = {S_cl:.6f}")
    print(f"  相位 exp(iS_cl/ħ) = {phase_cl.real:.4f} + {phase_cl.imag:.4f} i")
    print(f"  |exp(iS_cl/ħ)| = {abs(phase_cl):.6f}")

    # ── 生成所有随机路径 ──
    print(f"\n  == 生成 {N_PATHS} 条随机路径 == ")
    all_paths = []
    all_actions = []
    all_x_cl = []

    for i in range(N_PATHS):
        path, x_cl_i = generate_random_path(t_grid, XI, XF, TI, TF,
                                            N_MODES, FOURIER_SIGMA)
        all_paths.append(path)
        all_x_cl.append(x_cl_i)
        action = compute_action(path, t_grid, MASS, HBAR, OMEGA)
        all_actions.append(action)

    all_actions = np.array(all_actions)
    all_paths = np.array(all_paths)

    # ── 相位 ──
    all_phases = np.exp(1j * all_actions / HBAR)
    S_min = np.min(all_actions)
    S_max = np.max(all_actions)
    S_std = np.std(all_actions)

    print(f"  S 范围: [{S_min:.4f}, {S_max:.4f}]")
    print(f"  S_cl = {S_cl:.4f}")
    print(f"  S_std = {S_std:.4f}")
    print(f"  相位聚集度: N_cluster = "
          f"{np.sum(np.abs(all_actions - S_cl) < S_std)}/{N_PATHS}")
    print(f"  离散路径积分近似: |K| ≈ {abs(discretized_propagator(all_paths, all_actions, MASS, HBAR, dt)):.6f}")

    # ── 按 |S - S_cl| 排序路径 (用于颜色映射) ──
    dS = np.abs(all_actions - S_cl)
    sort_idx = np.argsort(dS)

    # ══════════════════════════════════════════════════════════════
    #  绘图
    # ══════════════════════════════════════════════════════════════

    # 中文字体
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC',
                                        'Source Han Sans CN',
                                        'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['mathtext.fontset'] = 'dejavusans'

    fig = plt.figure(figsize=(18, 8), facecolor='#0d1117')

    # ── 颜色方案 ──
    # 使用从蓝色 (靠近经典) → 品红 (远离经典) 的渐变
    # 对每条路径根据其 |S - S_cl| 着色
    dS_norm = dS / (dS.max() + 1e-12)

    def path_color(i):
        """根据路径与经典作用量的距离返回颜色 (色相从青到红)"""
        dn = dS_norm[i]
        # 色相: 0.55 (青/靠近经典) → 0.0 (红/远离经典)
        hue = 0.55 * (1.0 - dn)
        return hsv_to_rgb((hue, 0.85, 0.9))

    # ══════════════════════════════════════════════════════════════
    # Panel 1: 样本路径 + 经典路径
    # ══════════════════════════════════════════════════════════════
    ax1 = fig.add_subplot(1, 2, 1, facecolor='#0d1117')

    # 画 N_SHOW 条路径 (按与经典距离排序, 均匀采样)
    show_indices = np.linspace(0, N_PATHS - 1, N_SHOW, dtype=int)
    show_indices = sort_idx[show_indices]  # 按距离采样的路径索引

    for idx in show_indices:
        ax1.plot(t_grid, all_paths[idx], '-', color=path_color(idx),
                 linewidth=0.6, alpha=0.7, zorder=2)

    # 经典路径 (粗白线, 虚线描边)
    ax1.plot(t_grid, x_cl, '-', color='#ffffff', linewidth=3.0,
             zorder=5, label='经典路径 δS=0')
    ax1.plot(t_grid, x_cl, '-', color='#3fb950', linewidth=1.2,
             alpha=0.6, zorder=4)

    # 端点标记
    ax1.scatter([TI], [XI], c='#58a6ff', s=80, zorder=6,
                edgecolors='white', linewidths=1.5)
    ax1.scatter([TF], [XF], c='#ff7b72', s=80, zorder=6,
                edgecolors='white', linewidths=1.5)
    ax1.annotate(f'({TI},{XI})', xy=(TI, XI),
                 xytext=(TI + 0.15, XI + 0.2),
                 color='#58a6ff', fontsize=10, fontweight='bold')
    ax1.annotate(f'({TF:.1f},{XF})', xy=(TF, XF),
                 xytext=(TF - 0.6, XF + 0.25),
                 color='#ff7b72', fontsize=10, fontweight='bold')

    # 势阱背景提示
    V_max = 0.5 * MASS * OMEGA**2 * max(abs(XI), abs(XF), 1.2)**2
    ax1.fill_between([TI, TF], -1.0, 5.0, alpha=0.04, color='#58a6ff')
    ax1.text(0.5 * (TI + TF), 1.3,
             r'V(x) = ½mω²x²',
             color='#8b949e', fontsize=10, alpha=0.5, ha='center')

    ax1.set_xlabel('时间 t', color='#e6edf3', fontsize=13)
    ax1.set_ylabel('位置 x(t)', color='#e6edf3', fontsize=13)
    ax1.set_title(f'样本路径 (共{N_PATHS}条, 显示{N_SHOW}条)  +  经典路径',
                  color='#e6edf3', fontsize=13, fontweight='bold')
    ax1.legend(loc='upper right', facecolor='#161b22',
               edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
    ax1.tick_params(colors='#e6edf3')
    ax1.grid(True, alpha=0.10, color='#30363d')
    ax1.set_xlim(TI - 0.1, TF + 0.1)

    # ══════════════════════════════════════════════════════════════
    # Panel 2: 稳定相位聚类 — 复平面单位圆
    # ══════════════════════════════════════════════════════════════
    ax2 = fig.add_subplot(1, 2, 2, facecolor='#0d1117')

    # 单位圆
    theta_circle = np.linspace(0, 2 * np.pi, 300)
    ax2.plot(np.cos(theta_circle), np.sin(theta_circle),
             '-', color='#30363d', linewidth=1.0, alpha=0.6, zorder=1)

    # 画所有路径的相位点
    # 靠近经典的 (dS 小) → 大点, 亮色
    # 远离经典的 (dS 大) → 小点, 暗色
    for i in range(N_PATHS):
        phase = all_phases[i]
        dn = dS_norm[i]
        # 大小: 从 60 (靠近经典) 到 8 (远离经典)
        size = 60.0 * np.exp(-3.0 * dn) + 6.0
        # 透明度: 从 0.9 到 0.15
        alpha = 0.9 * np.exp(-2.5 * dn) + 0.1
        # 颜色: 用色相编码 |S - S_cl|
        hue = 0.15 * dn  # 0 (绿/靠近) → 0.15 (黄/远)
        color = hsv_to_rgb((hue, 0.7, 0.95))

        ax2.scatter(phase.real, phase.imag, s=size,
                    c=[color], alpha=alpha, zorder=2,
                    edgecolors='none')

    # 经典路径相位 (特殊标记)
    ax2.scatter([phase_cl.real], [phase_cl.imag], s=250,
                c='#ffffff', zorder=6, edgecolors='#3fb950',
                linewidths=2.5, marker='*')
    ax2.annotate('S=S_cl', xy=(phase_cl.real, phase_cl.imag),
                 xytext=(phase_cl.real + 0.25, phase_cl.imag + 0.2),
                 color='#3fb950', fontsize=11, fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#3fb950', lw=1.5))

    # 相位均值 (传播子方向)
    phase_mean = np.mean(all_phases)
    ax2.scatter([phase_mean.real], [phase_mean.imag], s=180,
                c='#d2a8ff', zorder=5, edgecolors='white',
                linewidths=1.5, marker='D')
    ax2.annotate('⟨exp(iS/ħ)⟩', xy=(phase_mean.real, phase_mean.imag),
                 xytext=(phase_mean.real + 0.15, phase_mean.imag - 0.3),
                 color='#d2a8ff', fontsize=10, fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#d2a8ff', lw=1.2))

    # 轴与标记
    ax2.axhline(0, color='#30363d', linewidth=0.5, alpha=0.4)
    ax2.axvline(0, color='#30363d', linewidth=0.5, alpha=0.4)
    ax2.set_xlabel('Re[exp(iS/ħ)]', color='#e6edf3', fontsize=13)
    ax2.set_ylabel('Im[exp(iS/ħ)]', color='#e6edf3', fontsize=13)
    ax2.set_title('稳定相位聚类: exp(iS/ħ) 在复平面单位圆上',
                  color='#e6edf3', fontsize=13, fontweight='bold')
    ax2.set_aspect('equal')
    ax2.set_xlim(-1.3, 1.3)
    ax2.set_ylim(-1.3, 1.3)
    ax2.tick_params(colors='#e6edf3')

    # 说明文本框
    n_cluster = np.sum(dS < S_std)
    cluster_text = (
        f"稳定相位原理:\n"
        f"● 靠近 S_cl 的路径 (约 {n_cluster} 条):  相位相近, 加强干涉\n"
        f"● 远离 S_cl 的路径:  相位随机, 相互抵消\n"
        f"● S_cl = {S_cl:.3f},  σ_S = {S_std:.3f}\n"
        f"● |⟨exp(iS/ħ)⟩| = {abs(phase_mean):.3f}"
    )
    ax2.text(0.02, 0.98, cluster_text,
             transform=ax2.transAxes,
             verticalalignment='top',
             fontsize=8.5, color='#8b949e',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#161b22',
                       edgecolor='#30363d', alpha=0.85))

    # ── 全局样式 ──
    for ax in [ax1, ax2]:
        for spine in ax.spines.values():
            spine.set_color('#30363d')

    fig.suptitle(
        f'费曼路径积分 — 对所有路径求和  |  '
        f'谐振子 V(x)=½mω²x²  |  '
        f'm={MASS}, ω={OMEGA}, ħ={HBAR}  |  '
        f'x({TI})={XI}, x({TF:.1f})={XF}  |  '
        f'{N_PATHS} 条路径, {N_T} 时间步',
        color='#e6edf3', fontsize=11, fontweight='bold', y=1.01)

    plt.tight_layout(rect=(0, 0, 1, 0.95))

    # ── 保存 ──
    save_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, 'path_integral_1d.png')
    fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    plt.close(fig)

    print(f"\n  [OK] 图像已保存: {save_path}")

    # ── 物理总结 ──
    print(f"\n{'='*62}")
    print("  物理结论")
    print(f"{'='*62}")
    print(f"""
  费曼路径积分的核心洞见:

  (1) 量子传播子 = 对所有可能路径求和:
      K(x_f, x_i; T) = ∫ D[x(t)] exp(iS[x]/ħ)

  (2) 稳定相位近似 (ħ→0):
      在经典极限下, 只有 δS=0 的路径 (经典路径) 对传播子
      有净贡献 — 其它路径的相位 exp(iS/ħ) 快速振荡, 相互抵消。

  (3) 数值验证 ({N_PATHS} 条路径):
      ● 经典作用量 S_cl = {S_cl:+.4f}
      ● 作用量标准差 σ_S = {S_std:.4f}
      ● 约 {n_cluster}/{N_PATHS} 条路径在 S_cl ± σ_S 范围内
      ● 这些路径的相位聚集于单位圆上的小弧段
      ● 远离经典的路径相位随机散布, 贡献趋近于零

  (4) 这解释了经典物理如何从量子物理「涌现」:
      ● 宏观世界 → 作用量巨大 (S ≫ ħ) → 只有 δS=0 的路径幸存
      ● 微观世界 → 作用量 ∼ ħ → 量子涨落显著
      ● 路径积分统一了经典力学 (极值原理) 与量子力学 (叠加原理)

  (5) 谐振子的特殊性:
      ● 只有经典路径对传播子有精确贡献 (Van Vleck 行列式为常数)
      ● 所有量子涨落精确抵消 → 半经典近似对谐振子是精确的!
      ● 传播子: K = (mω/2πiħ sin ωT)^{1/2} exp(iS_cl/ħ)

  对应书籍: §7.1 传播子概念 + §7.2 谐振子路径积分
""")

    return {
        'S_cl': S_cl,
        'actions': all_actions,
        'phases': all_phases,
        'phase_mean': phase_mean,
        'n_cluster': n_cluster,
    }


if __name__ == '__main__':
    result = main()

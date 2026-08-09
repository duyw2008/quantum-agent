#!/usr/bin/env python3
"""
Dyson 级数渐近收敛性 — φ⁴ 理论微扰展开的渐近行为演示

物理核心:
    量子场论的微扰展开是渐近级数 (asymptotic series):
    - 低阶项逐步趋近精确结果
    - 超过「最优截断」阶数后发散 (∼n! 增长)
    - 这是 φ⁴ 理论的基本特征, 不是数值误差

演示设计:
    Panel 1: Dyson 部分和 |Σ A_k|² 收敛到精确值 + 单阶 |A_n|² 下降趋势
    Panel 2: 相对误差 vs 阶数 (对数) + ∼n! 理论参考线
    在截断 Fock 空间 (dim=27) 中, n! 发散被部分压制;
    ∼n! 参考线展示完整理论中的预期发散行为。

对应书籍 §4.2 + §5.1 预备
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import math
import numpy as np
from src.qft import LatticePhi4


# ================================================================
# Dyson 级数高效计算 (相互作用绘景递归, 返回复振幅)
# ================================================================

def dyson_complex_amplitudes(H0: np.ndarray, V: np.ndarray,
                              psi_i: np.ndarray, psi_f: np.ndarray,
                              t: float, max_order: int,
                              n_steps: int = 100) -> list:
    """使用相互作用绘景递归计算各阶复振幅 <f|U^{(n)}(t)|i>

    递归:  |psi_I^{(0)}(t)> = |i>
           |psi_I^{(n)}(t)> = -i ∫_0^t ds V_I(s) |psi_I^{(n-1)}(s)>

    返回:  [A_0, A_1, ..., A_max]  各阶复振幅
    """
    dim = H0.shape[0]
    dt = t / n_steps

    eigvals, eigvecs = np.linalg.eigh(H0)

    def U0(tau):
        return eigvecs @ np.diag(np.exp(-1j * eigvals * tau)) @ eigvecs.conj().T

    # 预计算 V_I(t_j)
    VI = [U0(j * dt).conj().T @ V @ U0(j * dt) for j in range(n_steps + 1)]

    psi_I = [np.tile(psi_i, (n_steps + 1, 1))]
    complex_amps = [complex(np.dot(psi_f.conj(), U0(t) @ psi_i))]

    for n in range(1, max_order + 1):
        psi_curr = np.zeros((n_steps + 1, dim), dtype=complex)
        acc = np.zeros(dim, dtype=complex)
        for j in range(1, n_steps + 1):
            contrib = 0.5 * (VI[j-1] @ psi_I[n-1][j-1] +
                             VI[j]   @ psi_I[n-1][j]) * dt
            acc += contrib
            psi_curr[j] = -1j * acc
        psi_I.append(psi_curr)
        complex_amps.append(
            complex(np.dot(psi_f.conj(), U0(t) @ psi_curr[n_steps])))

    return complex_amps


def exact_probability(H_full: np.ndarray, psi_i: np.ndarray,
                      psi_f: np.ndarray, t: float) -> float:
    """精确对角化: |<f|e^{-iHt}|i>|^2"""
    eigvals, eigvecs = np.linalg.eigh(H_full)
    U_exact = eigvecs @ np.diag(np.exp(-1j * eigvals * t)) @ eigvecs.conj().T
    return float(abs(np.dot(psi_f.conj(), U_exact @ psi_i))**2)


# ================================================================
# 主程序
# ================================================================

def main():
    print("=" * 62)
    print("  Dyson 级数渐近行为 — φ^4 理论微扰展开")
    print("=" * 62)

    # ── 模型参数 ──
    N_sites = 3
    mass = 0.5
    coupling = 1.0
    N_fock = 3
    t = 1.0
    max_order = 10
    n_steps = 100

    dim = N_fock ** N_sites
    print(f"\n  格点: N={N_sites}, N_fock={N_fock} -> Hilbert dim = {dim}")
    print(f"  m={mass},  lambda={coupling},  t={t}")
    print(f"  最大阶数: {max_order},  时间步数: {n_steps}")

    # ── 构建模型 ──
    lpt = LatticePhi4(N_sites, mass=mass, coupling=0.0, N_fock=N_fock)
    print(f"  {lpt.summary()}")

    H0 = lpt.hamiltonian(0.0)
    H_full = lpt.hamiltonian(coupling)
    V = H_full - H0

    V_norm = np.linalg.norm(V)
    H0_norm = np.linalg.norm(H0)
    ratio = V_norm / H0_norm
    print(f"  ||V|| / ||H0|| = {V_norm:.4f} / {H0_norm:.4f} = {ratio:.4f}")

    # ── 初态 |i> = |1,0,0>,  末态 |f> = |0,1,0> ──
    # |n0,n1,n2>  idx = n0*9 + n1*3 + n2
    i_idx = 1 * 9 + 0 * 3 + 0   # |1,0,0>
    f_idx = 0 * 9 + 1 * 3 + 0   # |0,1,0>

    psi_i = np.zeros(dim, dtype=complex); psi_i[i_idx] = 1.0 + 0j
    psi_f = np.zeros(dim, dtype=complex); psi_f[f_idx] = 1.0 + 0j

    print(f"  |i> = |1,0,0> (idx={i_idx}),  |f> = |0,1,0> (idx={f_idx})")

    # ── 精确对角化 ──
    print("\n  == 精确对角化 ==")
    exact_prob = exact_probability(H_full, psi_i, psi_f, t)
    print(f"  |<f|e^{{-iHt}}|i>|^2 = {exact_prob:.10e}")

    # ── Dyson 级数 ──
    print(f"\n  == Dyson 级数 (0 ~ {max_order} 阶) ==")
    complex_amps = dyson_complex_amplitudes(H0, V, psi_i, psi_f, t,
                                            max_order, n_steps)

    indiv_probs = np.array([abs(a)**2 for a in complex_amps])
    cumul_amps = np.cumsum(complex_amps)
    cumul_probs = np.abs(cumul_amps)**2

    print(f"  {'n':>3s}  {'|A_n|^2':>14s}  {'|Sum A_k|^2':>14s}"
          f"  {'error':>14s}")
    print(f"  {'-'*3}  {'-'*14}  {'-'*14}  {'-'*14}")
    for n in range(max_order + 1):
        err = cumul_probs[n] - exact_prob
        print(f"  {n:3d}  {indiv_probs[n]:14.6e}  {cumul_probs[n]:14.10f}"
              f"  {err:+.2e}")

    # ── 最优截断 ──
    errors = np.abs(cumul_probs - exact_prob)
    opt_order = int(np.argmin(errors))
    print(f"\n  最优截断阶数: N* = {opt_order}")
    print(f"  最小误差:     Delta_min = {errors[opt_order]:.2e}")

    # ── 寻找单阶贡献的最小项位置 ──
    if len(indiv_probs) > 1:
        indiv_min_order = int(np.argmin(indiv_probs[1:])) + 1
    else:
        indiv_min_order = 0

    # ================================================================
    # 绘图
    # ================================================================
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # 中文字体
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC',
                                        'Source Han Sans CN',
                                        'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    # 使用 mathtext 渲染数学符号, 避免 CJK 字体缺失 Unicode 数学字符
    plt.rcParams['mathtext.fontset'] = 'dejavusans'

    orders = np.arange(0, max_order + 1)

    fig = plt.figure(figsize=(16, 6.5), facecolor='#0d1117')

    # ══════════════════════════════════════════════════════════════
    # Panel 1: 部分和收敛 + 单阶贡献 |A_n|^2
    # ══════════════════════════════════════════════════════════════
    ax1 = fig.add_subplot(1, 2, 1, facecolor='#0d1117')
    ax1b = ax1.twinx()

    # 主 y 轴 (左): 部分和 |Σ A_k|^2  (线性)
    ax1.plot(orders, cumul_probs, 'o-', color='#58a6ff', linewidth=2.2,
             markersize=8, label='Dyson 部分和', zorder=4)
    ax1.axhline(exact_prob, color='#3fb950', linewidth=2.5, linestyle='-',
                alpha=0.9, label=f'精确值 = {exact_prob:.4e}', zorder=2)

    # 最优截断标记
    ax1.axvline(opt_order, color='#d2a8ff', linewidth=1.2,
                linestyle=':', alpha=0.7)
    ax1.plot(opt_order, cumul_probs[opt_order], 'D', color='#3fb950',
             markersize=12, markeredgecolor='white', markeredgewidth=1.5,
             zorder=6)

    # 右 y 轴: 单阶贡献 |A_n|^2 (对数)
    ax1b.plot(orders, np.maximum(indiv_probs, 1e-20), 's--',
              color='#f0883e', linewidth=1.5, markersize=7,
              alpha=0.8, label='单阶 |A_n|^2', zorder=3)
    ax1b.set_yscale('log')

    # 标注
    ax1.annotate(
        f'N*={opt_order}\n误差={errors[opt_order]:.1e}',
        xy=(opt_order, cumul_probs[opt_order]),
        xytext=(opt_order + 1.8, cumul_probs[opt_order] * 0.85),
        arrowprops=dict(arrowstyle='->', color='#d2a8ff', lw=1.5),
        color='#d2a8ff', fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#161b22',
                  edgecolor='#30363d', alpha=0.9))

    ax1.set_xlabel('微扰阶数 n', color='#e6edf3', fontsize=12)
    ax1.set_ylabel('部分和 |Sum A_k|^2', color='#58a6ff', fontsize=12)
    ax1b.set_ylabel('单阶贡献 |A_n|^2 (对数)', color='#f0883e', fontsize=11)
    ax1.set_title('Dyson 级数收敛到精确值',
                  color='#e6edf3', fontsize=13, fontweight='bold')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='center right', facecolor='#161b22',
               edgecolor='#30363d', labelcolor='#e6edf3', fontsize=8.5)

    ax1.tick_params(colors='#58a6ff')
    ax1b.tick_params(colors='#f0883e')
    ax1.grid(True, alpha=0.12, color='#30363d')
    ax1.set_xlim(-0.5, max_order + 0.5)

    # ══════════════════════════════════════════════════════════════
    # Panel 2: 相对误差 (对数) + ~n! 理论发散
    # ══════════════════════════════════════════════════════════════
    ax2 = fig.add_subplot(1, 2, 2, facecolor='#0d1117')

    rel_err = np.where(errors > 0, errors / exact_prob,
                       np.full_like(errors, 1e-16))

    ax2.plot(orders, rel_err, 'o-', color='#58a6ff', linewidth=2.2,
             markersize=8, label='数值相对误差', zorder=4)
    ax2.set_yscale('log')

    # 最优截断
    ax2.axvline(opt_order, color='#d2a8ff', linewidth=1.2,
                linestyle=':', alpha=0.7)
    ax2.plot(opt_order, max(rel_err[opt_order], 1e-16), 'D',
             color='#3fb950', markersize=14, markeredgecolor='white',
             markeredgewidth=1.5, zorder=6)

    # ∼n! 理论参考线
    factorial_vals = np.array([math.factorial(n) if n > 0 else 1.0
                                for n in orders], dtype=float)
    # 归一化: 在 lowest-individual-term 位置匹配 scale
    ref_order = min(indiv_min_order, opt_order) if opt_order > 0 else 1
    if ref_order > 0:
        nf_scaled = (factorial_vals / factorial_vals[ref_order]
                     * rel_err[ref_order])
    else:
        nf_scaled = factorial_vals / factorial_vals[1] * rel_err[1]

    ax2.plot(orders, nf_scaled, '--', color='#ff7b72', linewidth=1.8,
             alpha=0.7, label='~n! 理论增长 (完整理论)', zorder=3)
    ax2.fill_between(orders, nf_scaled * 0.1, nf_scaled * 10,
                     alpha=0.06, color='#ff7b72')

    # 区域标注
    mid_converge = max(1, opt_order // 2)
    ax2.annotate('收敛区', xy=(mid_converge, rel_err[mid_converge] * 0.2),
                 color='#3fb950', fontsize=11, fontweight='bold',
                 ha='center', alpha=0.85)
    if opt_order < max_order - 1:
        mid_div = (opt_order + max_order) / 2
        ax2.annotate('截断稳定区\n(完整理论中 ~n! 发散)',
                     xy=(mid_div, rel_err[int(mid_div)] * 0.3),
                     color='#8b949e', fontsize=9, ha='center', alpha=0.75)

    ax2.set_xlabel('微扰阶数 n', color='#e6edf3', fontsize=12)
    ax2.set_ylabel('相对误差 (对数)', color='#e6edf3', fontsize=12)
    ax2.set_title('误差: 快速收敛 + ~n! 理论预期发散',
                  color='#e6edf3', fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right', facecolor='#161b22',
               edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
    ax2.tick_params(colors='#e6edf3')
    ax2.grid(True, alpha=0.12, color='#30363d', which='both')
    ax2.set_xlim(-0.5, max_order + 0.5)

    # ── 全局样式 ──
    for ax in [ax1, ax2]:
        for spine in ax.spines.values():
            spine.set_color('#30363d')

    fig.suptitle(
        f'φ^4 理论 Dyson 级数渐近收敛性  |  '
        f'N={N_sites}, m={mass}, λ={coupling}, t={t}  |  '
        f'Hilbert dim={dim}  |  最优截断 N*={opt_order}',
        color='#e6edf3', fontsize=11, fontweight='bold', y=1.02)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # ── 保存 ──
    save_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, 'dyson_convergence.png')
    fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    plt.close(fig)

    print(f"\n  [OK] 图像已保存: {save_path}")

    # ── 物理总结 ──
    print(f"\n{'='*62}")
    print("  物理结论")
    print(f"{'='*62}")
    print(f"""
  φ^4 理论的微扰展开是渐近级数 (asymptotic series):

  (1) 低阶 (n <= {opt_order}):   逐阶快速逼近精确值
      相对误差从 {rel_err[0]:.2e} 降至 {rel_err[opt_order]:.2e}

  (2) 最优截断 N* = {opt_order}:   误差最小 = {errors[opt_order]:.2e}
      超过此阶后, 更多微扰项不再改善精度

  (3) 理论上 (完整 Fock 空间):  更高阶项 ∼n! 增长, 级数发散
      本截断空间 (dim={dim}) 压制了发散, 但 ∼n! 参考线
      展示了完整理论中的预期行为

  根本原因:
    * 大阶数 Feynman 图数量 ∼ n!  (组合爆炸)
    * Dyson 论证:  微扰级数收敛半径为 0
    * 物理含义:  不能通过取更多微扰项逼近真实物理
    * 解决途径:  格点理论、非微扰方法 (Borel 求和等)

  对应书籍:  §4.2 Dyson 级数 + §5.1 格点 φ^4 理论
""")

    return {
        'exact_prob': exact_prob,
        'complex_amps': complex_amps,
        'indiv_probs': indiv_probs,
        'cumul_probs': cumul_probs,
        'errors': errors,
        'opt_order': opt_order,
    }


if __name__ == '__main__':
    result = main()

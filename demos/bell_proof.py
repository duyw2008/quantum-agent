#!/usr/bin/env python3
"""
Bell's Theorem: Proof of Nonlocality
=====================================
分步证明贝尔不等式违反意味着定域隐变量不存在。

五步动画：
  1. 假设：存在定域隐变量 λ → A(a,λ)=±1, B(b,λ)=±1
  2. 推导：CHSH 不等式 |S| ≤ 2
  3. 量子：Ψ⁻ 态 → E(a,b) = −cos(a−b)
  4. 代入：最优角 → S = 2√2 ≈ 2.828
  5. 矛盾：2.828 > 2 → 定域隐变量不存在

输出: demos/output/bell_proof.mp4
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.animation import FuncAnimation
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# CJK font
_cn_font_path = '/usr/share/fonts/adobe-source-han-sans/SourceHanSansCN-Regular.otf'
if os.path.exists(_cn_font_path):
    fm.fontManager.addfont(_cn_font_path)
    _cn_prop = fm.FontProperties(fname=_cn_font_path)
    plt.rcParams['font.family'] = _cn_prop.get_name()
    print(f"  Using CJK font: {_cn_prop.get_name()}")
else:
    _cn_prop = None
    print("  Warning: CJK font not found, Chinese may not render")


def dark_axes(ax):
    ax.set_facecolor('#0d1117')
    ax.tick_params(colors='#8b949e', labelsize=9)
    for spine in ax.spines.values():
        spine.set_color('#30363d')
    ax.grid(True, alpha=0.1, color='#30363d')


def step1_setup(fig):
    """假设定域隐变量"""
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2])
    ax_l = fig.add_subplot(gs[0]); dark_axes(ax_l)
    ax_r = fig.add_subplot(gs[1]); dark_axes(ax_r)

    # Left: source + Alice + Bob diagram
    ax_l.set_xlim(-5, 5); ax_l.set_ylim(-3, 3)
    ax_l.set_title('EPR 实验设置', color='#c9d1d9', fontsize=13, fontweight='bold')
    # Source
    ax_l.plot(0, 0, 'o', color='#f0883e', markersize=18, zorder=5)
    ax_l.text(0, -0.5, '源\nΨ⁻', color='#f0883e', fontsize=10, ha='center', fontweight='bold')
    # Particles
    ax_l.annotate('', xy=(-3, 0.3), xytext=(-1.5, 0.1),
                  arrowprops=dict(arrowstyle='->', color='#79c0ff', lw=2))
    ax_l.annotate('', xy=(3, 0.3), xytext=(1.5, 0.1),
                  arrowprops=dict(arrowstyle='->', color='#ff7b72', lw=2))
    ax_l.text(-3.2, 0.6, '光子 A', color='#79c0ff', fontsize=11, fontweight='bold')
    ax_l.text(2.7, 0.6, '光子 B', color='#ff7b72', fontsize=11, fontweight='bold')
    # Alice
    ax_l.add_patch(plt.Rectangle((-4.5, -0.5), 1.5, 2, fill=False,
                                  edgecolor='#79c0ff', linewidth=2))
    ax_l.text(-3.75, 1.1, 'Alice\n测 a', color='#79c0ff', fontsize=11, ha='center', fontweight='bold')
    # Bob
    ax_l.add_patch(plt.Rectangle((3, -0.5), 1.5, 2, fill=False,
                                  edgecolor='#ff7b72', linewidth=2))
    ax_l.text(3.75, 1.1, 'Bob\n测 b', color='#ff7b72', fontsize=11, ha='center', fontweight='bold')
    # Hidden variable
    ax_l.annotate('λ (隐变量)', xy=(0, -1.5), fontsize=12, color='#d2a8ff',
                  ha='center', fontweight='bold',
                  bbox=dict(boxstyle='round', facecolor='#d2a8ff20', edgecolor='#d2a8ff', alpha=0.8))
    ax_l.axis('off')

    # Right: assumptions
    ax_r.set_xlim(0, 10); ax_r.set_ylim(0, 10)
    ax_r.axis('off')
    ax_r.set_title('定域隐变量假设', color='#c9d1d9', fontsize=13, fontweight='bold')

    lines = [
        (1, 8.5, '定域性: Alice 的结果不依赖 Bob 的设置', '#79c0ff'),
        (1, 7.3, '  A = A(a, λ)    B = B(b, λ)', '#8b949e'),
        (1, 5.8, '实在性: 测量结果由 λ 预先决定', '#ff7b72'),
        (1, 4.6, '  A, B ∈ {+1, −1}', '#8b949e'),
        (1, 3.1, 'λ 的概率分布 ρ(λ) ≥ 0, ∫ρ(λ)dλ = 1', '#d2a8ff'),
        (1, 1.8, '关联: E(a,b) = ∫ A(a,λ)B(b,λ) ρ(λ) dλ', '#c9d1d9'),
    ]
    for x, y, text, color in lines:
        ax_r.text(x, y, text, color=color, fontsize=12, fontfamily='monospace',
                  va='center', fontweight='bold' if y in [8.5, 5.8, 3.1, 1.8] else 'normal')

    return fig


def step2_derive(fig):
    """推导 CHSH ≤ 2"""
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[0]); dark_axes(ax)

    ax.set_xlim(0, 14); ax.set_ylim(0, 12)
    ax.axis('off')

    ax.text(7, 11.5, 'Step 2: 推导 CHSH 不等式', color='#c9d1d9', fontsize=16,
            fontweight='bold', ha='center')

    derivation = [
        (1.0, 10.0, '定义 S = A(a₁)[B(b₁)+B(b₂)] + A(a₂)[B(b₁)−B(b₂)]', '#f0883e'),
        (1.5, 8.8, '因为 B=±1, 括号中一个为 0 (±2), 另一个为 ±2 (或 0)', '#8b949e'),
        (2.0, 7.6, '所以 S = A(a₁)·(0 or ±2) + A(a₂)·(±2 or 0) = ±2', '#8b949e'),
        (2.5, 6.4, '因为 A=±1, 两项的和只能是 ±2 或 0', '#8b949e'),
        (3.0, 5.2, '因此 |S| ≤ 2 对任何单次测量成立', '#79c0ff'),
        (3.5, 3.8, '求期望: |⟨S⟩| = |E(a₁,b₁)+E(a₁,b₂)+E(a₂,b₁)−E(a₂,b₂)| ≤ 2', '#3fb950'),
        (4.0, 2.2, '⇒ |S_CHSH| ≤ 2   Q.E.D.', '#3fb950'),
    ]
    for x, y, text, color in derivation:
        ax.text(x, y, text, color=color, fontsize=12, fontfamily='monospace', va='center')

    return fig


def step3_quantum(fig):
    """量子力学预测"""
    gs = fig.add_gridspec(1, 2)
    ax_l = fig.add_subplot(gs[0]); dark_axes(ax_l)
    ax_r = fig.add_subplot(gs[1]); dark_axes(ax_r)

    # Left: the quantum state
    ax_l.set_xlim(0, 10); ax_l.set_ylim(0, 10)
    ax_l.axis('off')
    ax_l.set_title('量子力学: Ψ⁻ 态', color='#c9d1d9', fontsize=14, fontweight='bold')

    qm_text = [
        (0.5, 9.0, '|Ψ⁻⟩ = (|01⟩ − |10⟩)/√2', '#79c0ff'),
        (0.5, 7.5, '测量 Aₐ = σ(a) = cos(a)σz + sin(a)σx', '#8b949e'),
        (0.5, 6.0, '测量 B_b = σ(b) = cos(b)σz + sin(b)σx', '#8b949e'),
        (0.5, 4.3, '关联:', '#c9d1d9'),
        (1.0, 3.0, 'E(a,b) = ⟨Ψ⁻|σ(a)⊗σ(b)|Ψ⁻⟩', '#f0883e'),
        (1.0, 1.7, 'E(a,b) = −cos(a − b)', '#3fb950'),
    ]
    for x, y, text, color in qm_text:
        fs = 15 if '|Ψ⁻⟩' in text else 13 if 'E(a,b)' in text else 11
        ax_l.text(x, y, text, color=color, fontsize=fs, fontfamily='monospace', va='center')

    # Right: correlation curve
    theta = np.linspace(0, np.pi, 200)
    eq = -np.cos(theta)
    e_classical = 1 - 2 * theta / np.pi  # optimal LHV

    ax_r.plot(theta * 180 / np.pi, eq, color='#79c0ff', linewidth=3, label='Quantum: −cos(θ)')
    ax_r.plot(theta * 180 / np.pi, e_classical, color='#ff7b72', linewidth=2,
              linestyle='--', label='Classical limit (LHV)')
    ax_r.set_xlabel('θ = |a−b| [°]', color='#c9d1d9', fontsize=11)
    ax_r.set_ylabel('E', color='#c9d1d9', fontsize=11)
    ax_r.set_title('关联函数 E(θ)', color='#c9d1d9', fontsize=13, fontweight='bold')
    ax_r.legend(loc='lower left', facecolor='#161b22', edgecolor='#30363d',
                labelcolor='#c9d1d9', fontsize=9)
    ax_r.set_ylim(-1.1, 1.1)
    ax_r.axhline(0, color='#30363d', linewidth=0.5)
    ax_r.axvline(45, color='#30363d', linewidth=0.5)
    ax_r.axvline(135, color='#30363d', linewidth=0.5)

    # Mark violation zone (where quantum deviates from classical)
    ax_r.fill_between(theta * 180 / np.pi, eq, e_classical, where=(eq < e_classical),
                       color='#f0883e20', alpha=0.4, label='Violation zone')
    ax_r.text(90, -0.3, '违规区', color='#f0883e', fontsize=10, ha='center', fontstyle='italic')

    return fig


def step4_plug_in(fig):
    """代入最优角度"""
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[0]); dark_axes(ax)
    ax.set_xlim(0, 14); ax.set_ylim(0, 12)
    ax.axis('off')

    ax.text(7, 11.5, 'Step 4: 代入最优测量角', color='#c9d1d9', fontsize=16,
            fontweight='bold', ha='center')

    steps = [
        (1.0, 10.2, 'Alice: a₁=0°, a₂=90°', '#79c0ff'),
        (1.0, 9.0, 'Bob:   b₁=45°, b₂=135°', '#ff7b72'),
        (1.5, 7.6, 'E(0°,45°)   = −cos(−45°)  = −1/√2', '#8b949e'),
        (1.5, 6.4, 'E(0°,135°)  = −cos(−135°) = +1/√2', '#8b949e'),
        (1.5, 5.2, 'E(90°,45°)  = −cos(45°)   = −1/√2', '#8b949e'),
        (1.5, 4.0, 'E(90°,135°) = −cos(−45°)  = −1/√2', '#ff7b72'),
        (2.0, 2.5, 'S = (−1/√2) + (1/√2) + (−1/√2) − (−1/√2)', '#f0883e'),
        (2.5, 1.0, '= −2/√2 = −√2 → |S| = 2√2 ≈ 2.828', '#3fb950'),
    ]
    for x, y, text, color in steps:
        fs = 14 if '2.828' in text else 12
        ax.text(x, y, text, color=color, fontsize=fs, fontfamily='monospace', va='center')

    return fig


def step5_conclusion(fig):
    """结论"""
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.2])
    ax_top = fig.add_subplot(gs[0]); dark_axes(ax_top)
    ax_bot = fig.add_subplot(gs[1]); dark_axes(ax_bot)

    # Top: S value comparison
    ax_top.set_xlim(-0.5, 1.5); ax_top.set_ylim(0, 3.5)
    labels = ['Classical\nBound', 'Quantum\n(Ψ⁻)']
    values = [2, 2.828]
    colors = ['#30363d', '#f0883e']
    bars = ax_top.bar(labels, values, color=colors, width=0.5)
    ax_top.bar_label(bars, fmt='%.3f', color='#c9d1d9', fontsize=16, fontweight='bold')
    ax_top.set_ylabel('|S|', color='#c9d1d9', fontsize=13)
    ax_top.set_title('Step 5: 矛盾！', color='#ff7b72', fontsize=16, fontweight='bold')
    ax_top.axhline(2, color='#ff7b72', linewidth=2, linestyle='--', alpha=0.7)
    ax_top.fill_between([-0.5, 1.5], 2, 3.3, color='#ff7b7220', alpha=0.3)
    ax_top.text(0.5, 2.88, '2.828 > 2\n违反贝尔不等式!', color='#f0883e',
                fontsize=13, ha='center', fontweight='bold')
    ax_top.set_ylim(0, 3.4)

    # Bottom: conclusion text
    ax_bot.set_xlim(0, 14); ax_bot.set_ylim(0, 6)
    ax_bot.axis('off')

    ax_bot.add_patch(FancyBboxPatch((0.3, 0.3), 13.4, 5.2,
                                     boxstyle='round,pad=0.3',
                                     facecolor='#f0883e15', edgecolor='#f0883e',
                                     linewidth=2, alpha=0.8))

    conclusion = [
        (0.8, 4.8, '结论:', '#c9d1d9'),
        (1.5, 3.5, '量子力学预测 |S| = 2√2 ≈ 2.828 > 2', '#f0883e'),
        (1.5, 2.2, '定域隐变量要求 |S| ≤ 2', '#ff7b72'),
        (2.2, 0.8, '⇒  定域隐变量理论不能复现量子力学预测', '#3fb950'),
    ]
    for x, y, text, color in conclusion:
        fs = 14 if '结论' in text else 13
        ax_bot.text(x, y, text, color=color, fontsize=fs, fontfamily='monospace',
                    va='center', fontweight='bold' if '结论' in text else 'normal')

    return fig


def main():
    print("Generating Bell Proof animation...")

    steps = [step1_setup, step2_derive, step3_quantum, step4_plug_in, step5_conclusion]
    step_names = [
        '1. 假设定域隐变量',
        '2. 推导 CHSH ≤ 2',
        '3. 量子力学预测',
        '4. 代入最优角度',
        '5. 矛盾！定域性不成立',
    ]

    # Pre-render all step figures
    step_figs = []
    for i, step_fn in enumerate(steps):
        fig = plt.figure(figsize=(14, 7), facecolor='#0d1117')
        step_fn(fig)
        fig.suptitle(f'Bell-CHSH Proof: {step_names[i]}',
                     color='#c9d1d9', fontsize=15, fontweight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.94])
        step_figs.append(fig)

    # Create animation: hold each step for some frames
    hold_frames = 60  # frames per step
    n_frames = len(steps) * hold_frames

    # Main animation figure
    anim_fig = plt.figure(figsize=(14, 7), facecolor='#0d1117')

    def draw_step_on_anim_fig(step_idx):
        anim_fig.clear()
        step_fn = steps[step_idx]
        step_fn(anim_fig)
        anim_fig.suptitle(f'Bell-CHSH Proof: {step_names[step_idx]}',
                          color='#c9d1d9', fontsize=15, fontweight='bold', y=0.98)

    def update(frame):
        step_idx = min(frame // hold_frames, len(steps) - 1)
        anim_fig.clear()
        steps[step_idx](anim_fig)
        anim_fig.suptitle(f'Bell-CHSH Proof: {step_names[step_idx]}',
                          color='#c9d1d9', fontsize=15, fontweight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.94])
        return []

    ani = FuncAnimation(anim_fig, update, frames=n_frames, interval=60, blit=False)

    outdir = os.path.join(os.path.dirname(__file__) or '.', 'output')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'bell_proof.mp4')
    try:
        ani.save(path, writer='ffmpeg', fps=15, dpi=120, bitrate=2500)
    except Exception:
        ani.save(path, writer='pillow', fps=8, dpi=100)
    plt.close('all')
    print(f"  Saved: {path}")


if __name__ == '__main__':
    main()

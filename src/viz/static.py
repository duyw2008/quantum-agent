"""静态可视化 — 能谱、相空间、矩阵元素"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import Optional, Tuple, List

from ..matrix.operators import NumericOperatorSystem


def plot_energy_levels(energies: np.ndarray, labels: List[str] = None,
                       dark: bool = True, title: str = "Energy Spectrum",
                       save_path: str = None):
    """绘制能级图"""
    from .animate import QM_DARK_THEME, QM_LIGHT_THEME
    theme = QM_DARK_THEME if dark else QM_LIGHT_THEME

    fig, ax = plt.subplots(figsize=(8, 6), facecolor=theme['bg'])
    ax.set_facecolor(theme['bg'])

    n = len(energies)
    for i, E in enumerate(energies):
        ax.hlines(y=E, xmin=0.1, xmax=0.9, color=theme['energy'],
                  linewidth=max(2, 8 - i * 0.5))
        label = labels[i] if labels else f'n={i}'
        ax.text(0.95, E, f'  {label}: {E:.4f}', color=theme['fg'],
                va='center', fontsize=11, fontfamily='monospace')

    ax.set_xlim(0, 1.5)
    ax.set_ylabel('Energy', color=theme['fg'], fontsize=12)
    ax.set_title(title, color=theme['fg'], fontsize=14)
    ax.set_xticks([])
    ax.tick_params(colors=theme['fg'])
    ax.grid(True, axis='y', alpha=0.2, color=theme['grid'])

    for spine in ax.spines.values():
        spine.set_color(theme['grid'])

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor=theme['bg'])
        plt.close(fig)
    return fig, ax


def plot_matrix_element(matrix: np.ndarray, title: str = "Matrix Element",
                        dark: bool = True, save_path: str = None,
                        abs_only: bool = False):
    """绘制矩阵元素的绝对值热图

    参数:
        matrix:    N×N 矩阵
        title:    图标题
        dark:     暗色主题
        save_path: 保存路径
        abs_only: 仅显示幅值 (不显示实部/虚部)
    """
    from .animate import QM_DARK_THEME, QM_LIGHT_THEME
    theme = QM_DARK_THEME if dark else QM_LIGHT_THEME

    n_plots = 1 if abs_only else 3
    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5),
                              facecolor=theme['bg'])

    if n_plots == 1:
        axes = [axes]

    # 绝对值
    im0 = axes[0].imshow(np.abs(matrix), cmap='plasma', aspect='equal')
    axes[0].set_title(f'|{title}|', color=theme['fg'], fontsize=12)
    axes[0].set_xlabel('j', color=theme['fg'])
    axes[0].set_ylabel('i', color=theme['fg'])
    axes[0].tick_params(colors=theme['fg'])
    plt.colorbar(im0, ax=axes[0])

    if not abs_only:
        # 实部
        vmax = max(abs(matrix.real).max(), abs(matrix.imag).max(), 1e-10)
        im1 = axes[1].imshow(matrix.real, cmap='RdBu_r', aspect='equal',
                              vmin=-vmax, vmax=vmax)
        axes[1].set_title(f'Re[{title}]', color=theme['fg'], fontsize=12)
        axes[1].set_xlabel('j', color=theme['fg'])
        axes[1].set_ylabel('i', color=theme['fg'])
        axes[1].tick_params(colors=theme['fg'])
        plt.colorbar(im1, ax=axes[1])

        # 虚部
        im2 = axes[2].imshow(matrix.imag, cmap='RdBu_r', aspect='equal',
                              vmin=-vmax, vmax=vmax)
        axes[2].set_title(f'Im[{title}]', color=theme['fg'], fontsize=12)
        axes[2].set_xlabel('j', color=theme['fg'])
        axes[2].set_ylabel('i', color=theme['fg'])
        axes[2].tick_params(colors=theme['fg'])
        plt.colorbar(im2, ax=axes[2])

    for ax in axes:
        ax.set_facecolor(theme['bg'])
        for spine in ax.spines.values():
            spine.set_color(theme['grid'])

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor=theme['bg'])
        plt.close(fig)
    return fig


def plot_phase_space(wf: 'WaveFunction', dark: bool = True,
                     save_path: str = None):
    """绘制 Wigner 准概率分布的近似 (动量 vs 位置散点)

    简化版: 显示 Δx, Δp 的不确定度椭圆。
    """
    from .animate import QM_DARK_THEME, QM_LIGHT_THEME
    theme = QM_DARK_THEME if dark else QM_LIGHT_THEME

    x_mean = wf.expectation_x()
    p_mean = wf.expectation_p()
    dx = wf.uncertainty_x()
    dp = wf.uncertainty_p()

    fig, ax = plt.subplots(figsize=(8, 8), facecolor=theme['bg'])
    ax.set_facecolor(theme['bg'])

    # 不确定度椭圆
    from matplotlib.patches import Ellipse
    ellipse = Ellipse((x_mean, p_mean), width=2 * dx, height=2 * dp,
                      edgecolor=theme['potential'], facecolor='none',
                      linewidth=2, linestyle='--')
    ax.add_patch(ellipse)

    # 中心点
    ax.plot(x_mean, p_mean, 'o', color=theme['prob_density'], markersize=8)

    ax.set_xlabel('x', color=theme['fg'], fontsize=12)
    ax.set_ylabel('p', color=theme['fg'], fontsize=12)
    ax.set_title(f'Phase Space | Δx·Δp = {dx * dp:.4f} ≥ 0.5',
                 color=theme['fg'], fontsize=14)
    ax.tick_params(colors=theme['fg'])
    ax.grid(True, alpha=0.2, color=theme['grid'])
    ax.axhline(y=0, color=theme['grid'], linewidth=0.5)
    ax.axvline(x=0, color=theme['grid'], linewidth=0.5)
    ax.set_aspect('equal', adjustable='datalim')

    for spine in ax.spines.values():
        spine.set_color(theme['grid'])

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor=theme['bg'])
        plt.close(fig)
    return fig, ax

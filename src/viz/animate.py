"""可视化模块 — 波函数动画与静态图

提供两类可视化:
    1. animate_wavefunction() — 生成波函数时间演化的 MP4/GIF 动画
    2. 静态图函数 — 势函数、本征态、能谱、相空间等

依赖: matplotlib, Pillow (用于 GIF)
"""

import numpy as np
from typing import Optional, Tuple, List
import os
import sys

import matplotlib
matplotlib.use('Agg')  # 非交互后端
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

from ..core.wave_function import WaveFunction, Grid
from ..core.potentials import Potential
from ..core.schrodinger import EvolutionResult


# ============================================================
# 配色方案
# ============================================================

# 量子力学主题配色 (暗色)
QM_DARK_THEME = {
    'bg': '#0d1117',
    'fg': '#e6edf3',
    'grid': '#30363d',
    'potential': '#58a6ff',
    'prob_density': '#79c0ff',
    'real_part': '#ff7b72',
    'imag_part': '#7ee787',
    'phase': '#d2a8ff',
    'energy': '#f0883e',
    'ket': '#ffa657',
}

# 量子力学亮色主题
QM_LIGHT_THEME = {
    'bg': '#ffffff',
    'fg': '#1f2328',
    'grid': '#d0d7de',
    'potential': '#0550ae',
    'prob_density': '#0969da',
    'real_part': '#cf222e',
    'imag_part': '#1a7f37',
    'phase': '#8250df',
    'energy': '#bc4c00',
    'ket': '#bf8700',
}


def get_theme(dark: bool = True) -> dict:
    """获取配色主题"""
    return QM_DARK_THEME if dark else QM_LIGHT_THEME


# ============================================================
# 静态图
# ============================================================

def plot_potential(V: Potential, x_range: Tuple[float, float] = (-5, 5),
                   n_points: int = 1000, ax=None, dark: bool = True,
                   title: str = None, save_path: str = None):
    """绘制势函数 V(x)

    参数:
        V:         势函数对象
        x_range:   x 轴范围
        n_points:  采样点数
        ax:        可选的 matplotlib axes
        dark:      使用暗色主题
        title:     图表标题
        save_path: 保存路径
    """
    theme = get_theme(dark)
    x = np.linspace(*x_range, n_points)
    Vx = V(x)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5), facecolor=theme['bg'])
    else:
        fig = ax.figure

    ax.set_facecolor(theme['bg'])
    ax.fill_between(x, Vx, Vx.min() - 0.1, alpha=0.3, color=theme['potential'])
    ax.plot(x, Vx, color=theme['potential'], linewidth=2)
    ax.axhline(y=0, color=theme['grid'], linestyle='--', alpha=0.5)
    ax.set_xlabel('x', color=theme['fg'], fontsize=12)
    ax.set_ylabel('V(x)', color=theme['fg'], fontsize=12)
    ax.set_title(title or f'Potential: {V.label}', color=theme['fg'], fontsize=14)
    ax.tick_params(colors=theme['fg'])
    ax.grid(True, alpha=0.2, color=theme['grid'])

    for spine in ax.spines.values():
        spine.set_color(theme['grid'])

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor=theme['bg'])
        plt.close(fig)

    return fig, ax


def plot_wavefunction(psi: WaveFunction, V: Potential = None,
                      show_potential: bool = True,
                      show_phase: bool = True,
                      dark: bool = True,
                      title: str = None, save_path: str = None):
    """绘制波函数快照

    显示:
        - 概率密度 |ψ(x)|² (填充)
        - 实部 Re[ψ(x)] (实线)
        - 虚部 Im[ψ(x)] (虚线)
        - 可选: 势函数 V(x) (背景)
        - 可选: 相位色条
    """
    theme = get_theme(dark)
    x = psi.grid.x

    fig, axes = plt.subplots(2 if show_phase else 1, 1,
                              figsize=(12, 8), facecolor=theme['bg'],
                              sharex=True,
                              gridspec_kw={'height_ratios': [2, 1]} if show_phase else None)

    if show_phase:
        ax_main, ax_phase = axes
    else:
        ax_main = axes

    ax_main.set_facecolor(theme['bg'])

    # 势函数
    if show_potential and V is not None:
        Vx = V(x)
        # 缩放到合理的显示范围
        v_max = max(np.abs(Vx).max(), 1.0)
        Vx_norm = Vx / v_max
        ax_main.fill_between(x, Vx_norm * 0.3 * psi.probability_density.max(),
                             alpha=0.15, color=theme['potential'], label='V(x)')

    # 概率密度
    prob = psi.probability_density
    ax_main.fill_between(x, prob, alpha=0.4, color=theme['prob_density'],
                         label=r'$|\psi(x)|^2$')

    # 实部和虚部
    psi_real = psi.psi.real
    psi_imag = psi.psi.imag
    scale = np.max(prob) ** 0.5
    ax_main.plot(x, psi_real * scale, color=theme['real_part'],
                 linewidth=1.5, label=r'$\mathrm{Re}[\psi(x)]$')
    ax_main.plot(x, psi_imag * scale, color=theme['imag_part'],
                 linewidth=1.5, linestyle='--', label=r'$\mathrm{Im}[\psi(x)]$')

    ax_main.axhline(y=0, color=theme['grid'], linewidth=0.5)
    ax_main.set_ylabel(r'$|\psi|^2$, $\psi$', color=theme['fg'], fontsize=12)
    ax_main.legend(loc='upper right', facecolor=theme['bg'],
                   edgecolor=theme['grid'], labelcolor=theme['fg'])
    ax_main.tick_params(colors=theme['fg'])
    ax_main.grid(True, alpha=0.2, color=theme['grid'])

    # 标题含统计信息
    title_text = title or f'Wavefunction at t = {psi.t:.3f}'
    title_text += f' | ⟨x⟩ = {psi.expectation_x():.3f}, Δx = {psi.uncertainty_x():.3f}'
    ax_main.set_title(title_text, color=theme['fg'], fontsize=14)

    # 相位
    if show_phase:
        ax_phase.set_facecolor(theme['bg'])
        phase = psi.phase
        # 只显示概率密度显著区域的相位
        mask = prob > 0.01 * prob.max()
        ax_phase.plot(x, phase, color=theme['phase'], linewidth=1)
        ax_phase.set_xlabel('x', color=theme['fg'], fontsize=12)
        ax_phase.set_ylabel(r'$\phi(x)$', color=theme['fg'], fontsize=12)
        ax_phase.set_ylim(-np.pi, np.pi)
        ax_phase.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        ax_phase.set_yticklabels([r'$-\pi$', r'$-\pi/2$', '0',
                                   r'$\pi/2$', r'$\pi$'])
        ax_phase.tick_params(colors=theme['fg'])
        ax_phase.grid(True, alpha=0.2, color=theme['grid'])

    for ax in (fig.axes if hasattr(fig, 'axes') else [fig.gca()]):
        for spine in ax.spines.values():
            spine.set_color(theme['grid'])

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor=theme['bg'])
        plt.close(fig)

    return fig


def plot_eigenstates(grid: Grid, V: Potential, n_states: int = 5,
                     dark: bool = True, save_path: str = None):
    """绘制前 n 个本征态的能量和波函数

    使用数值对角化计算本征态。
    """
    from scipy.linalg import eigh_tridiagonal

    theme = get_theme(dark)
    N = grid.n_points
    dx = grid.dx
    x = grid.x
    hbar = 1.0
    mass = 1.0

    # 构建哈密度量
    ke_diag = np.ones(N) * hbar**2 / (mass * dx**2)
    ke_off = np.ones(N - 1) * (-hbar**2 / (2 * mass * dx**2))
    Vx = V(x)
    diag = ke_diag + Vx

    eigvals, eigvecs = eigh_tridiagonal(diag, ke_off,
                                         select='i', select_range=(0, n_states - 1))

    fig, (ax_en, ax_wf) = plt.subplots(1, 2, figsize=(14, 8),
                                        facecolor=theme['bg'],
                                        gridspec_kw={'width_ratios': [1, 3]})

    # 左: 能级图
    ax_en.set_facecolor(theme['bg'])
    for n, E in enumerate(eigvals):
        ax_en.axhline(y=E, color=theme['energy'], linewidth=2)
        ax_en.text(n_states - 0.5, E, f'  E{n}={E:.4f}',
                   color=theme['fg'], va='center', fontsize=10,
                   fontfamily='monospace')
    ax_en.set_xlim(-0.5, n_states - 0.5)
    ax_en.set_ylabel('Energy', color=theme['fg'], fontsize=12)
    ax_en.set_title('Energy Levels', color=theme['fg'], fontsize=14)
    ax_en.tick_params(colors=theme['fg'])
    ax_en.grid(True, alpha=0.2, color=theme['grid'])

    # 右: 波函数
    ax_wf.set_facecolor(theme['bg'])
    offset_scale = (eigvals[1] - eigvals[0]) * 0.8 if n_states > 1 else 0.1
    for n in range(n_states):
        psi_n = eigvecs[:, n]
        prob_n = psi_n**2
        offset = eigvals[n]
        ax_wf.fill_between(x, offset, offset + prob_n / prob_n.max() * offset_scale,
                           alpha=0.4, color=theme['prob_density'],
                           label=f'|{n}⟩')
        ax_wf.plot(x, offset + psi_n / np.abs(psi_n).max() * offset_scale,
                   color=theme['real_part'], linewidth=1)
    ax_wf.set_xlabel('x', color=theme['fg'], fontsize=12)
    ax_wf.set_ylabel('Energy / Wavefunction', color=theme['fg'], fontsize=12)
    ax_wf.set_title('Eigenstate Wavefunctions', color=theme['fg'], fontsize=14)
    ax_wf.tick_params(colors=theme['fg'])
    ax_wf.grid(True, alpha=0.2, color=theme['grid'])
    ax_wf.legend(loc='upper right', facecolor=theme['bg'],
                 edgecolor=theme['grid'], labelcolor=theme['fg'])

    for ax in [ax_en, ax_wf]:
        for spine in ax.spines.values():
            spine.set_color(theme['grid'])

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor=theme['bg'])
        plt.close(fig)

    return fig, eigvals, eigvecs


# ============================================================
# 动画
# ============================================================

def animate_evolution(result: EvolutionResult, V: Potential,
                      save_path: str = None,
                      fps: int = 30, duration: float = None,
                      dpi: int = 150, dark: bool = True,
                      show_momentum: bool = False,
                      progress_callback=None) -> animation.FuncAnimation:
    """为 EvolutionResult 创建动画

    生成 MP4/GIF 动画，显示:
        - 上方: 概率密度 |ψ(x, t)|² 的演化 (heatmap 风格)
        - 下方: 势函数 + 当前概率密度 + ⟨x⟩ 轨迹
        - 可选: 动量空间分布

    参数:
        result:             演化结果
        V:                  势函数
        save_path:          保存路径 (.mp4 或 .gif)
        fps:                帧率
        duration:           动画时长 (秒), 默认匹配原始演化时间
        dpi:                分辨率
        dark:               暗色主题
        show_momentum:      同时显示动量空间
        progress_callback:  进度回调 (用于 agent 交互)
    """
    theme = get_theme(dark)
    x = result.grid.x
    Vx = V(x)
    prob = result.prob_density_snapshots
    times = result.times

    n_snapshots = prob.shape[0]

    # 预计算全局范围
    prob_max = prob.max()
    v_min, v_max = Vx.min(), Vx.max()
    v_scale = max(abs(v_min), abs(v_max), 1.0)

    # 帧间隔
    interval_ms = int(1000 / fps)

    # 图形尺寸
    n_rows = 3 if show_momentum else 2
    fig, axes = plt.subplots(n_rows, 1, figsize=(12, 4 * n_rows),
                              facecolor=theme['bg'],
                              gridspec_kw={'height_ratios': [1, 2, 2] if show_momentum else [1, 3]})

    if show_momentum:
        ax_obs, ax_main, ax_k = axes
    else:
        ax_obs, ax_main = axes

    k = result.grid.k
    dk = result.grid.dk

    # --- 上: 可观测量的时间演化 ---
    ax_obs.set_facecolor(theme['bg'])
    ax_obs.plot(times, result.expectation_x, color=theme['real_part'],
                linewidth=1.5, label=r'$\langle x \rangle$')
    ax_obs.plot(times, result.expectation_p, color=theme['imag_part'],
                linewidth=1.5, linestyle='--', label=r'$\langle p \rangle$')
    ax_obs.axhline(y=0, color=theme['grid'], linewidth=0.5)
    ax_obs.set_ylabel('Expectation', color=theme['fg'], fontsize=11)
    ax_obs.legend(loc='upper right', facecolor=theme['bg'],
                  edgecolor=theme['grid'], labelcolor=theme['fg'], fontsize=9)
    ax_obs.tick_params(colors=theme['fg'])
    ax_obs.grid(True, alpha=0.2, color=theme['grid'])
    # 添加时间指示线
    time_line = ax_obs.axvline(x=0, color=theme['potential'], linewidth=2, alpha=0.7)

    # --- 中: 概率密度 + 势 ---
    ax_main.set_facecolor(theme['bg'])
    ax2 = ax_main.twinx()

    # 势函数 (背景)
    Vx_norm = Vx / v_scale
    ax_main.fill_between(x, Vx_norm * prob_max * 0.4, alpha=0.2,
                         color=theme['potential'])
    ax_main.plot(x, Vx_norm * prob_max * 0.4, color=theme['potential'],
                 linewidth=1, alpha=0.5)
    ax_main.set_ylabel(r'$|\psi(x,t)|^2$', color=theme['fg'], fontsize=12)
    ax_main.set_xlabel('x', color=theme['fg'], fontsize=12)
    ax_main.tick_params(colors=theme['fg'])
    ax_main.grid(True, alpha=0.2, color=theme['grid'])

    # 期望值标记
    exp_x_line = ax_main.axvline(x=result.expectation_x[0],
                                  color=theme['potential'], linewidth=2,
                                  linestyle=':', alpha=0.8)
    prob_fill = ax_main.fill_between(x, 0, prob[0], alpha=0.5,
                                      color=theme['prob_density'])
    prob_line, = ax_main.plot(x, prob[0], color=theme['prob_density'],
                               linewidth=1.5)

    ax2.set_ylabel('V(x)', color=theme['potential'], fontsize=12)
    ax2.tick_params(colors=theme['potential'])

    # --- 下 (可选): 动量空间 ---
    if show_momentum:
        ax_k.set_facecolor(theme['bg'])
        k_fill = ax_k.fill_between(k, 0, np.zeros_like(k), alpha=0.5,
                                    color=theme['phase'])
        k_line, = ax_k.plot(k, np.zeros_like(k), color=theme['phase'],
                            linewidth=1.5)
        ax_k.set_xlabel('k', color=theme['fg'], fontsize=12)
        ax_k.set_ylabel(r'$|\tilde{\psi}(k,t)|^2$', color=theme['fg'], fontsize=12)
        ax_k.tick_params(colors=theme['fg'])
        ax_k.grid(True, alpha=0.2, color=theme['grid'])

    # 标题
    title = ax_main.set_title('', color=theme['fg'], fontsize=14)

    for ax in axes:
        for spine in ax.spines.values():
            spine.set_color(theme['grid'])

    plt.tight_layout()

    # 使用 mutable 容器追踪 fill 对象 (避免 remove 后引用失效)
    prob_fill_ref = [prob_fill]
    k_fill_ref = [k_fill] if show_momentum else [None]

    # 更新函数
    def update(frame_idx):
        # 移除上一帧的 fill
        prob_fill_ref[0].remove()
        prob_fill_ref[0] = ax_main.fill_between(x, 0, prob[frame_idx], alpha=0.5,
                                                  color=theme['prob_density'])
        prob_line.set_ydata(prob[frame_idx])
        time_line.set_xdata([times[frame_idx], times[frame_idx]])
        exp_x_line.set_xdata([result.expectation_x[frame_idx],
                               result.expectation_x[frame_idx]])

        # 标题
        t = times[frame_idx]
        energy = result.energy[frame_idx]
        title.set_text(f't = {t:.3f} | E = {energy:.4f} | '
                       f'⟨x⟩ = {result.expectation_x[frame_idx]:.3f} | '
                       f'⟨p⟩ = {result.expectation_p[frame_idx]:.3f}')

        if show_momentum:
            # 计算动量空间分布
            psi_k_data = result.psi_snapshots[frame_idx]
            psi_k_fft = np.fft.fft(psi_k_data)
            prob_k = np.abs(psi_k_fft)**2
            k_fill_ref[0].remove()
            k_fill_ref[0] = ax_k.fill_between(k, 0, prob_k, alpha=0.5,
                                               color=theme['phase'])
            k_line.set_ydata(prob_k)

        if progress_callback and frame_idx % 10 == 0:
            progress_callback(frame_idx, n_snapshots)

        artists = [prob_line, time_line, exp_x_line, title]
        if show_momentum:
            artists.extend([k_line])
        return artists

    # 创建动画
    total_frames = n_snapshots
    ani = animation.FuncAnimation(fig, update, frames=total_frames,
                                   interval=interval_ms, blit=False,
                                   repeat=True)

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)

        if save_path.endswith('.gif'):
            writer = animation.PillowWriter(fps=fps)
            ani.save(save_path, writer=writer, dpi=dpi)
        else:
            # MP4
            try:
                writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
                ani.save(save_path, writer=writer, dpi=dpi)
            except (FileNotFoundError, RuntimeError):
                # 回退到 GIF
                gif_path = save_path.rsplit('.', 1)[0] + '.gif'
                writer = animation.PillowWriter(fps=fps)
                ani.save(gif_path, writer=writer, dpi=dpi)
                print(f"FFmpeg not found, saved as GIF: {gif_path}")
                save_path = gif_path

    return ani


def animate_probability_density(result: EvolutionResult,
                                 save_path: str = None,
                                 fps: int = 30, dpi: int = 150,
                                 dark: bool = True):
    """概率密度热图动画 (紧凑版)

    仅显示 |ψ(x,t)|² 的瀑布/热图演化。
    """
    theme = get_theme(dark)
    prob = result.prob_density_snapshots
    x = result.grid.x
    times = result.times
    n_snapshots = len(times)

    fig, ax = plt.subplots(figsize=(10, 6), facecolor=theme['bg'])
    ax.set_facecolor(theme['bg'])

    # 初始帧
    im = ax.imshow(prob[:1].T, aspect='auto', origin='lower',
                    extent=[times[0], times[-1], x[0], x[-1]],
                    cmap='plasma', vmin=0, vmax=prob.max())
    ax.set_xlabel('t', color=theme['fg'], fontsize=12)
    ax.set_ylabel('x', color=theme['fg'], fontsize=12)
    ax.set_title(r'$|\psi(x,t)|^2$', color=theme['fg'], fontsize=14)
    ax.tick_params(colors=theme['fg'])
    plt.colorbar(im, ax=ax, label=r'$|\psi|^2$').set_label(
        r'$|\psi|^2$', color=theme['fg'])

    def update(frame):
        # 显示到当前帧的所有数据
        im.set_data(prob[:frame + 1].T)

    interval = int(1000 / fps)
    ani = animation.FuncAnimation(fig, update, frames=n_snapshots,
                                   interval=interval)

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)
        ani.save(save_path, writer='pillow', fps=fps, dpi=dpi)

    return ani

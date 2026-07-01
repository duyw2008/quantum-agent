#!/usr/bin/env python3
"""
Bell-CHSH Angle Scan Animation
===============================
展示 |S(a₁,a₂)| 曲面如何突破经典上界 |S| ≤ 2。

固定 Bob 的最优角度 b₁=45°, b₂=135°,
扫描 Alice 的角度 a₁,a₂ ∈ [0°, 360°],
展示 3D 曲面 + 2D 热力图, 违规区高亮。

输出: demos/output/bell_scan.mp4
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


def chsh_S(a1_deg, a2_deg, b1_deg=45, b2_deg=135):
    """CHSH |S| for given angles (degrees), using Ψ⁻ correlation E=-cos(a-b)"""
    a1, a2 = np.radians(a1_deg), np.radians(a2_deg)
    b1, b2 = np.radians(b1_deg), np.radians(b2_deg)
    e11 = -np.cos(a1 - b1)
    e12 = -np.cos(a1 - b2)
    e21 = -np.cos(a2 - b1)
    e22 = -np.cos(a2 - b2)
    return np.abs(e11 + e12 + e21 - e22)


def main():
    # Grid
    n = 120
    a1_vals = np.linspace(0, 360, n)
    a2_vals = np.linspace(0, 360, n)
    A1, A2 = np.meshgrid(a1_vals, a2_vals)
    S_map = chsh_S(A1, A2)

    S_max_theory = 2 * np.sqrt(2)  # ~2.828

    fig = plt.figure(figsize=(16, 7), facecolor='#0d1117')

    # ── Left: 3D surface ──
    ax3d = fig.add_subplot(121, projection='3d', facecolor='#0d1117')
    surf = ax3d.plot_surface(A1, A2, S_map, cmap='plasma', alpha=0.85,
                             linewidth=0, antialiased=True)

    # Classical bound plane at z=2
    bound_z = np.full_like(A1, 2.0)
    ax3d.plot_surface(A1, A2, bound_z, color='#ff7b72', alpha=0.25,
                      linewidth=0, antialiased=True)

    # Max point
    max_idx = np.unravel_index(np.argmax(S_map), S_map.shape)
    ax3d.scatter([A1[max_idx]], [A2[max_idx]], [S_map[max_idx]],
                 color='#ffffff', s=80, edgecolor='#f0883e', linewidth=2, zorder=10)

    ax3d.set_xlabel('a₁ [°]', color='#79c0ff', fontsize=11)
    ax3d.set_ylabel('a₂ [°]', color='#79c0ff', fontsize=11)
    ax3d.set_zlabel('|S|', color='#f0883e', fontsize=11)
    ax3d.set_title('|S(a₁,a₂)| — 3D Surface', color='#c9d1d9', fontsize=14,
                   fontweight='bold', pad=12)
    ax3d.tick_params(colors='#8b949e', labelsize=8)
    ax3d.xaxis.pane.fill = False; ax3d.yaxis.pane.fill = False
    ax3d.zaxis.pane.fill = False
    ax3d.xaxis.pane.set_edgecolor('#0d1117')
    ax3d.yaxis.pane.set_edgecolor('#0d1117')
    ax3d.zaxis.pane.set_edgecolor('#0d1117')
    ax3d.grid(True, alpha=0.15, color='#30363d')

    # Colorbar for 3D
    cbar3d = fig.colorbar(surf, ax=ax3d, shrink=0.6, pad=0.08)
    cbar3d.ax.yaxis.set_tick_params(color='#8b949e', labelsize=8)
    cbar3d.set_label('|S|', color='#c9d1d9', fontsize=10)
    plt.setp(cbar3d.ax.yaxis.get_ticklabels(), color='#8b949e')

    # ── Right: 2D heatmap ──
    ax2d = fig.add_subplot(122, facecolor='#0d1117')
    im = ax2d.pcolormesh(A1, A2, S_map, cmap='plasma', shading='auto')

    # Violation contour
    ax2d.contour(A1, A2, S_map, levels=[2.0], colors='#ff7b72',
                 linewidths=2, linestyles='--')
    ax2d.contourf(A1, A2, S_map, levels=[2.0, S_max_theory + 0.1],
                  colors='#ff7b7220', alpha=0.3)

    # Max point
    ax2d.scatter([A1[max_idx]], [A2[max_idx]], color='#ffffff',
                 s=60, edgecolor='#f0883e', linewidth=2, zorder=10)
    ax2d.annotate(f'  Max: {S_map[max_idx]:.3f}',
                  (A1[max_idx], A2[max_idx]),
                  color='#ffffff', fontsize=10, fontweight='bold',
                  xytext=(15, -10), textcoords='offset points')

    ax2d.set_xlabel('a₁ [°]', color='#79c0ff', fontsize=11)
    ax2d.set_ylabel('a₂ [°]', color='#79c0ff', fontsize=11)
    ax2d.set_title('|S(a₁,a₂)| — Heatmap', color='#c9d1d9', fontsize=14,
                   fontweight='bold', pad=10)
    ax2d.tick_params(colors='#8b949e', labelsize=8)
    for spine in ax2d.spines.values():
        spine.set_color('#30363d')
    ax2d.grid(True, alpha=0.1, color='#30363d')

    # Colorbar for 2D
    cbar2d = fig.colorbar(im, ax=ax2d, shrink=0.85, pad=0.03)
    cbar2d.ax.yaxis.set_tick_params(color='#8b949e', labelsize=8)
    cbar2d.set_label('|S|', color='#c9d1d9', fontsize=10)
    plt.setp(cbar2d.ax.yaxis.get_ticklabels(), color='#8b949e')

    # Annotations
    fig.text(0.5, 0.96, f'Bell-CHSH Angle Scan  |  b₁=45°, b₂=135°  |  max |S| = {S_max_theory:.3f}',
             color='#c9d1d9', fontsize=13, fontweight='bold',
             ha='center', transform=fig.transFigure)

    # Legend text
    ax2d.text(0.02, 0.97,
              '--- Violation boundary (|S|=2)\n'
              '    Violation zone (|S|>2)',
              transform=ax2d.transAxes, fontsize=9, fontfamily='monospace',
              color='#8b949e', va='top')

    # ── Animate: rotate 3D view ──
    n_frames = 120
    def update(frame):
        # Rotate 3D view
        elev = 25 + 10 * np.sin(frame * 2 * np.pi / n_frames)
        azim = (frame * 360 / n_frames) % 360
        ax3d.view_init(elev=elev, azim=azim)
        # Blink the classical bound annotation
        return [surf]

    ani = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)

    # Save
    outdir = os.path.join(os.path.dirname(__file__) or '.', 'output')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'bell_scan.mp4')
    try:
        ani.save(path, writer='ffmpeg', fps=20, dpi=120, bitrate=2500)
    except Exception:
        ani.save(path, writer='pillow', fps=10, dpi=100)
    plt.close(fig)
    print(f"  Saved: {path}")


if __name__ == '__main__':
    main()

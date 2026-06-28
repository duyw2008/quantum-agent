#!/usr/bin/env python3
"""
Bloch 球演化动画 — Rabi振荡 / 退相干 / Larmor进动
===================================================
生成三种场景的3D动画, 展示量子态在Bloch球上的实时演化。

输出: demos/output/bloch_*.mp4

用法:
    python bloch_evolution.py            # 生成所有三个动画
    python bloch_evolution.py rabi       # 只生成Rabi振荡
    python bloch_evolution.py decoherence # 只生成退相干
    python bloch_evolution.py larmor     # 只生成Larmor进动
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# ── Bloch球绘制工具 ──────────────────────────────────────────

def draw_bloch_sphere(ax, alpha=0.08):
    """绘制半透明Bloch球框架"""
    u = np.linspace(0, 2*np.pi, 48)
    v = np.linspace(0, np.pi, 24)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, rstride=4, cstride=4,
                      color='#404060', alpha=0.3, linewidth=0.3)
    # 坐标轴
    for (s, e, c, lbl) in [((-1.3,0,0),(1.3,0,0),'#ff444488','x'),
                            ((0,-1.3,0),(0,1.3,0),'#44ff4488','y'),
                            ((0,0,-1.3),(0,0,1.3),'#4488ff88','z')]:
        ax.plot([s[0],e[0]],[s[1],e[1]],[s[2],e[2]], color=c, linewidth=0.8)
        ax.text(e[0]*1.15, e[1]*1.15, e[2]*1.15, lbl,
                color=c[:7], fontsize=10, fontweight='bold')
    # 赤道
    th = np.linspace(0, 2*np.pi, 100)
    ax.plot(np.cos(th), np.sin(th), np.zeros_like(th),
            color='#606080', linewidth=0.5, alpha=0.5)


def setup_axes(title_text):
    """设置暗色主题3D坐标轴"""
    fig = plt.figure(figsize=(9, 8), facecolor='#0d1117')
    ax = fig.add_subplot(111, projection='3d', facecolor='#0d1117')
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_zlim(-1.3, 1.3)
    ax.set_title(title_text, color='#c9d1d9', fontsize=14, pad=10)
    ax.set_xlabel('x', color='#ff6666'); ax.set_ylabel('y', color='#66ff66')
    ax.set_zlabel('z', color='#6688ff')
    ax.grid(False)
    ax.xaxis.pane.fill = False; ax.yaxis.pane.fill = False; ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('#0d1117')
    ax.yaxis.pane.set_edgecolor('#0d1117')
    ax.zaxis.pane.set_edgecolor('#0d1117')
    # 隐藏刻度数字
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    return fig, ax


def save_animation(fig, ani, filename):
    """保存动画为mp4"""
    outdir = os.path.join(os.path.dirname(__file__) or '.', 'output')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, filename)
    try:
        ani.save(path, writer='ffmpeg', fps=30, dpi=120, bitrate=2000)
    except Exception:
        ani.save(path, writer='pillow', fps=15, dpi=100)
    plt.close(fig)
    print(f"  Saved: {path}")
    return path


# ── 场景1: Rabi振荡 ─────────────────────────────────────────

def rabi_oscillation():
    """
    共振驱动下的Rabi振荡: |0⟩ ↔ |1⟩ 来回翻转。
    等效于Bloch矢量绕x轴旋转。
    """
    print("\n  [Rabi] Driving |0⟩ with resonant field...")
    fig, ax = setup_axes("Rabi Oscillation — |0⟩ ↔ |1⟩ flip under resonant drive")

    draw_bloch_sphere(ax)

    # 状态: 绕x轴旋转, Ω_R = 2π (1周期=1秒)
    n_frames = 120
    t = np.linspace(0, 4, n_frames)  # 4个周期
    omega_rabi = 2 * np.pi  # Rabi频率

    states = []
    for ti in t:
        angle = omega_rabi * ti
        # 初始 |0⟩, 绕x轴旋转
        ry, rz = -np.sin(angle), np.cos(angle)
        states.append((0, ry, rz))

    # 轨迹尾迹
    trail_n = 40
    trail_x, trail_y, trail_z = [], [], []

    arrow_line, = ax.plot([], [], [], color='#ffaa00', linewidth=2.5)
    arrow_head, = ax.plot([], [], [], 'o', color='#ffcc44', markersize=8)
    trail_line, = ax.plot([], [], [], color='#ffaa0040', linewidth=0.8)

    # 状态标签
    state_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes,
                           color='#c9d1d9', fontsize=11, fontfamily='monospace')

    def init():
        arrow_line.set_data([], []); arrow_line.set_3d_properties([])
        arrow_head.set_data([], []); arrow_head.set_3d_properties([])
        trail_line.set_data([], []); trail_line.set_3d_properties([])
        return arrow_line, arrow_head, trail_line, state_text

    def update(i):
        sx, sy, sz = states[i]
        arrow_line.set_data([0, sx], [0, sy])
        arrow_line.set_3d_properties([0, sz])
        arrow_head.set_data([sx], [sy])
        arrow_head.set_3d_properties([sz])

        trail_x.append(sx); trail_y.append(sy); trail_z.append(sz)
        if len(trail_x) > trail_n:
            trail_x.pop(0); trail_y.pop(0); trail_z.pop(0)
        trail_line.set_data(trail_x, trail_y)
        trail_line.set_3d_properties(trail_z)

        # 文本: 当前状态
        prob0 = (1 + sz) / 2  # |⟨0|ψ⟩|²
        prob1 = 1 - prob0
        theta = np.arccos(sz) if abs(sz) < 1 else 0
        phi = np.arctan2(sy, sx)
        state_text.set_text(
            f"t={t[i]:.2f}s  |0⟩:{prob0:.3f}  |1⟩:{prob1:.3f}\n"
            f"θ={np.degrees(theta):.0f}°  φ={np.degrees(phi):.0f}°"
        )
        return arrow_line, arrow_head, trail_line, state_text

    ani = FuncAnimation(fig, update, frames=n_frames, init_func=init,
                        interval=40, blit=True)
    save_animation(fig, ani, 'bloch_rabi.mp4')


# ── 场景2: 退相干 ───────────────────────────────────────────

def decoherence():
    """
    退相干: Bloch矢量从球面向球心收缩。
    展示纯态→混合态的过程 (T1 + T2 衰减)。
    """
    print("\n  [Decoherence] Superposition decaying to mixture...")
    fig, ax = setup_axes("Decoherence — Pure → Mixed (T₁ + T₂ decay)")

    draw_bloch_sphere(ax)

    n_frames = 150
    t = np.linspace(0, 3, n_frames)
    gamma = 1.0  # 衰减率

    # 初始: |+⟩ = (|0⟩+|1⟩)/√2, Bloch矢量 (1,0,0)
    states = []
    for ti in t:
        r = np.exp(-gamma * ti)  # Bloch半径
        # 同时T2衰减 (xy平面收缩) 比T1快
        r_xy = r * np.exp(-0.3 * gamma * ti)  # T2稍快
        rz = 0  # |+⟩的z分量为0
        states.append((r_xy, 0, rz))

    trail_x, trail_y, trail_z = [], [], []
    arrow_line, = ax.plot([], [], [], color='#00ccff', linewidth=2.5)
    arrow_head, = ax.plot([], [], [], 'o', color='#44ddff', markersize=8)
    trail_line, = ax.plot([], [], [], color='#00ccff30', linewidth=0.6)

    # 画纯度指示圈
    for r_ring in [0.75, 0.5, 0.25]:
        th = np.linspace(0, 2*np.pi, 60)
        ax.plot(r_ring*np.cos(th), r_ring*np.sin(th), np.zeros_like(th),
                color='#ffffff15', linewidth=0.3)

    state_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes,
                           color='#c9d1d9', fontsize=11, fontfamily='monospace')

    def init():
        arrow_line.set_data([], []); arrow_line.set_3d_properties([])
        arrow_head.set_data([], []); arrow_head.set_3d_properties([])
        trail_line.set_data([], []); trail_line.set_3d_properties([])
        return arrow_line, arrow_head, trail_line, state_text

    def update(i):
        sx, sy, sz = states[i]
        arrow_line.set_data([0, sx], [0, sy])
        arrow_line.set_3d_properties([0, sz])
        arrow_head.set_data([sx], [sy])
        arrow_head.set_3d_properties([sz])

        trail_x.append(sx); trail_y.append(sy); trail_z.append(sz)
        if len(trail_x) > 50:
            trail_x.pop(0); trail_y.pop(0); trail_z.pop(0)
        trail_line.set_data(trail_x, trail_y)
        trail_line.set_3d_properties(trail_z)

        r_vec = np.sqrt(sx**2 + sy**2 + sz**2)
        purity = (1 + r_vec**2) / 2
        state_text.set_text(
            f"t={t[i]:.2f}s  r={r_vec:.3f}  purity={purity:.3f}\n"
            f"T₂ decay: e^(-γ·t)  Bloch→0"
        )
        return arrow_line, arrow_head, trail_line, state_text

    ani = FuncAnimation(fig, update, frames=n_frames, init_func=init,
                        interval=30, blit=True)
    save_animation(fig, ani, 'bloch_decoherence.mp4')


# ── 场景3: Larmor进动 ──────────────────────────────────────

def larmor_precession():
    """
    Larmor进动: 静磁场B∥z下自旋绕z轴匀速旋转。
    H = ω σz/2 → U(t) = exp(-iωt σz/2)
    """
    print("\n  [Larmor] Spin precessing around B-field (z-axis)...")
    fig, ax = setup_axes("Larmor Precession — Spin in static B-field ∥ z")

    draw_bloch_sphere(ax)

    # 画B场方向指示
    ax.quiver(0, 0, 0, 0, 0, 1.5, color='#4488ff', linewidth=2,
              arrow_length_ratio=0.15, alpha=0.7)
    ax.text(0, 0.1, 1.55, 'B', color='#6688ff', fontsize=12)

    n_frames = 180
    t = np.linspace(0, 3, n_frames)
    omega_l = 2 * np.pi * 1.5  # Larmor频率

    # 初始态: θ=60° 偏离z轴
    theta0 = np.radians(60)
    states = []
    for ti in t:
        phi = omega_l * ti
        sx = np.sin(theta0) * np.cos(phi)
        sy = np.sin(theta0) * np.sin(phi)
        sz = np.cos(theta0)
        states.append((sx, sy, sz))

    trail_x, trail_y, trail_z = [], [], []
    arrow_line, = ax.plot([], [], [], color='#ff44ff', linewidth=2.5)
    arrow_head, = ax.plot([], [], [], 'o', color='#ff88ff', markersize=8)
    trail_line, = ax.plot([], [], [], color='#ff44ff30', linewidth=0.6)

    # 画进动锥面
    th_cone = np.linspace(0, 2*np.pi, 80)
    r_cone = np.sin(theta0)
    z_cone = np.cos(theta0)
    ax.plot(r_cone*np.cos(th_cone), r_cone*np.sin(th_cone),
            np.full_like(th_cone, z_cone),
            color='#ff44ff20', linewidth=0.5)

    state_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes,
                           color='#c9d1d9', fontsize=11, fontfamily='monospace')

    def init():
        arrow_line.set_data([], []); arrow_line.set_3d_properties([])
        arrow_head.set_data([], []); arrow_head.set_3d_properties([])
        trail_line.set_data([], []); trail_line.set_3d_properties([])
        return arrow_line, arrow_head, trail_line, state_text

    def update(i):
        sx, sy, sz = states[i]
        arrow_line.set_data([0, sx], [0, sy])
        arrow_line.set_3d_properties([0, sz])
        arrow_head.set_data([sx], [sy])
        arrow_head.set_3d_properties([sz])

        trail_x.append(sx); trail_y.append(sy); trail_z.append(sz)
        if len(trail_x) > 60:
            trail_x.pop(0); trail_y.pop(0); trail_z.pop(0)
        trail_line.set_data(trail_x, trail_y)
        trail_line.set_3d_properties(trail_z)

        phi_deg = np.degrees(np.arctan2(sy, sx)) % 360
        state_text.set_text(
            f"t={t[i]:.2f}s  ωₗ={omega_l:.1f} rad/s\n"
            f"θ={np.degrees(theta0):.0f}° (const)  φ={phi_deg:.0f}°"
        )
        return arrow_line, arrow_head, trail_line, state_text

    ani = FuncAnimation(fig, update, frames=n_frames, init_func=init,
                        interval=30, blit=True)
    save_animation(fig, ani, 'bloch_larmor.mp4')


# ── 入口 ────────────────────────────────────────────────────

def main():
    scenes = {'rabi', 'decoherence', 'larmor'}
    if len(sys.argv) > 1:
        target = sys.argv[1].lower()
        if target not in scenes:
            print(f"Unknown scene: {target}. Choose: {', '.join(sorted(scenes))}")
            sys.exit(1)
        scenes = {target}

    if 'rabi' in scenes:
        rabi_oscillation()
    if 'decoherence' in scenes:
        decoherence()
    if 'larmor' in scenes:
        larmor_precession()

    print(f"\n  All done! Files in demos/output/")


if __name__ == '__main__':
    main()

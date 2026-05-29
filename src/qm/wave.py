"""波函数动力学 — 一维含时薛定谔方程

Split-Step Fourier Method (SSFM) 求解:
    iℏ ∂ψ/∂t = -ℏ²/2m ∂²ψ/∂x² + V(x)ψ

支持自由粒子、任意势函数的波包演化。
"""

import numpy as np

# numpy 兼容
if not hasattr(np, 'trapezoid'):
    np.trapezoid = np.trapz  # type: ignore


class WaveGrid:
    """一维空间网格"""
    def __init__(self, x_min=-10, x_max=10, N=1024):
        self.x = np.linspace(x_min, x_max, N)
        self.dx = self.x[1] - self.x[0]
        self.N = N
        # 动量空间
        self.k = 2 * np.pi * np.fft.fftfreq(N, self.dx)


def gaussian_wavepacket(grid: WaveGrid, x0=0.0, p0=5.0, sigma=1.0,
                         hbar=1.0) -> np.ndarray:
    """高斯波包

    ψ(x) = (πσ²)^{-1/4} exp[-(x-x0)²/(2σ²) + i p0 x / ℏ]
    """
    x = grid.x
    norm = (np.pi * sigma**2)**(-0.25)
    psi = norm * np.exp(-(x - x0)**2 / (2 * sigma**2) + 1j * p0 * x / hbar)
    return psi


def evolve_ssfm(psi0: np.ndarray, grid: WaveGrid, V_func=None,
                dt=0.005, t_max=5.0, hbar=1.0, mass=1.0,
                snapshots=200):
    """Split-Step Fourier 时间演化

    数值方法: 对称 Strang 拆分
        e^{-iHΔt/ħ} ≈ e^{-iVΔt/2ħ} · e^{-iTΔt/ħ} · e^{-iVΔt/2ħ}
    每步: V/2 → FFT → T → IFFT → V/2, 误差 O(Δt³)
    详见 docs/NUMERICAL_METHODS.md §1

    返回:
        result = {'times': [...], 'psi': [...], 'prob': [...], 'x_exp': [...], ...}
    """
    x = grid.x
    k = grid.k
    Vx = np.zeros(grid.N) if V_func is None else V_func(x)

    psi = psi0.copy()
    n_steps = int(t_max / dt)
    interval = max(1, n_steps // snapshots)

    result = {
        'times': [], 'psi': [], 'prob': [],
        'x_exp': [], 'p_exp': [], 'dx': [], 'energy': [],
        'grid': grid,
    }

    def save(t):
        result['times'].append(t)
        result['psi'].append(psi.copy())
        prob = np.abs(psi)**2
        result['prob'].append(prob)
        result['x_exp'].append(np.trapezoid(x * prob, x))
        # momentum expectation
        psi_k = np.fft.fft(psi)
        prob_k = np.abs(psi_k)**2
        result['p_exp'].append(np.trapezoid(k * prob_k, k) / np.trapezoid(prob_k, k))
        result['dx'].append(np.sqrt(np.trapezoid(x**2 * prob, x) - result['x_exp'][-1]**2))
        # energy
        T = np.trapezoid(0.5 * hbar**2 * k**2 / mass * prob_k, k) / np.trapezoid(prob_k, k)
        V_exp = np.trapezoid(Vx * prob, x)
        result['energy'].append(T + V_exp)

    save(0.0)

    # 预计算 SSFM 相位因子
    pe_half = np.exp(-0.5j * Vx * dt / hbar)
    ke_full = np.exp(-1.0j * hbar * k**2 * dt / (2 * mass))

    for step in range(1, n_steps + 1):
        # Step 1: 半步势能
        psi = pe_half * psi
        # Step 2: FFT → 动量空间
        psi_k = np.fft.fft(psi)
        # Step 3: 动能演化 (全步)
        psi_k *= ke_full
        # Step 4: 逆 FFT
        psi = np.fft.ifft(psi_k)
        # Step 5: 半步势能
        psi = pe_half * psi

        if step % interval == 0:
            save(step * dt)

    return result


def animate_wave(result, save_path='wave_evolution.mp4', fps=20):
    """生成波函数演化动画"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    x = result['grid'].x
    prob = result['prob']
    times = result['times']
    x_exp = np.array(result['x_exp'])
    dx = np.array(result['dx'])

    fig, (ax_top, ax_main) = plt.subplots(2, 1, figsize=(12, 7),
                                           gridspec_kw={'height_ratios': [1, 3]},
                                           facecolor='#0d1117')

    # Top: ⟨x⟩ vs t
    ax_top.set_facecolor('#0d1117')
    ax_top.plot(times, x_exp, color='#79c0ff', linewidth=2)
    ax_top.fill_between(times, x_exp - dx, x_exp + dx, alpha=0.2, color='#79c0ff')
    t_line = ax_top.axvline(0, color='#ff7b72', linewidth=1.5, alpha=0.7)
    ax_top.set_ylabel(r'$\langle x \rangle$', color='#e6edf3', fontsize=11)
    ax_top.tick_params(colors='#e6edf3')
    ax_top.grid(True, alpha=0.15, color='#30363d')

    # Main: |ψ|²
    ax_main.set_facecolor('#0d1117')
    fill = ax_main.fill_between(x, 0, prob[0], alpha=0.5, color='#79c0ff')
    line, = ax_main.plot(x, prob[0], color='#79c0ff', linewidth=1.5)
    ax_main.set_xlabel('x', color='#e6edf3', fontsize=12)
    ax_main.set_ylabel(r'$|\psi(x,t)|^2$', color='#e6edf3', fontsize=12)
    ax_main.tick_params(colors='#e6edf3')
    ax_main.grid(True, alpha=0.15, color='#30363d')
    title = ax_main.set_title('', color='#e6edf3', fontsize=14)

    ax_main.set_xlim(x[0], x[-1])
    ax_main.set_ylim(0, max(p.max() for p in prob) * 1.1)

    for ax in [ax_top, ax_main]:
        for spine in ax.spines.values():
            spine.set_color('#30363d')

    fill_ref = [fill]  # mutable container for update

    def update(i):
        fill_ref[0].remove()
        fill_ref[0] = ax_main.fill_between(x, 0, prob[i], alpha=0.5, color='#79c0ff')
        line.set_ydata(prob[i])
        t_line.set_xdata([times[i], times[i]])
        title.set_text(f't = {times[i]:.2f}  |  '
                       f'⟨x⟩ = {x_exp[i]:.3f}  |  '
                       f'Δx = {dx[i]:.3f}  |  '
                       f'E = {result["energy"][i]:.4f}')
        return [line, t_line, title]

    ani = animation.FuncAnimation(fig, update, frames=len(times),
                                   interval=1000 // fps, blit=False)

    import os
    os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)
    try:
        writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
        ani.save(save_path, writer=writer, dpi=150)
    except (FileNotFoundError, RuntimeError):
        gif_path = save_path.rsplit('.', 1)[0] + '.gif'
        writer = animation.PillowWriter(fps=fps)
        ani.save(gif_path, writer=writer, dpi=150)
        print(f"Saved as GIF: {gif_path}")
        save_path = gif_path

    plt.close(fig)
    print(f"Animation saved: {save_path}")
    return save_path


# ================================================================
# 特色势函数工厂
# ================================================================

def double_well(grid: WaveGrid, a=3.0, depth=10.0, separation=6.0):
    """双阱势: V(x) = depth * [(x/separation)² - 1]² / 2
    
    参数:
        a:        阱宽 (控制抛物线形状)
        depth:    阱深
        separation: 两阱间距
    """
    x = grid.x
    s = separation
    return depth * ((x/s)**2 - 1)**2 / 2


def periodic_potential(grid: WaveGrid, amplitude=2.0, period=4.0):
    """周期势 (余弦光晶格): V(x) = amplitude * cos(2π x / period)
    
    参数:
        amplitude: 势的幅度
        period:    空间周期
    """
    x = grid.x
    return amplitude * np.cos(2 * np.pi * x / period)


def delta_barrier(grid: WaveGrid, x0=0.0, strength=10.0):
    """δ 势垒的近似: V(x) = strength / (π σ²) * exp(-(x-x0)²/σ²)
    
    用窄高斯近似 δ 函数。σ 自动取为网格间距的 5 倍。
    
    参数:
        x0:       势垒位置
        strength: 势垒强度 (积分 ∫ V dx ≈ strength)
    """
    x = grid.x
    sigma = grid.dx * 5
    return strength / (np.sqrt(np.pi) * sigma) * np.exp(-(x - x0)**2 / sigma**2)


def finite_well(grid: WaveGrid, x0=0.0, width=4.0, depth=5.0):
    """有限深方势阱: V(x) = -depth (|x-x0| < width/2), 0 otherwise
    
    参数:
        x0:     阱中心
        width:  阱宽度
        depth:  阱深 (正值)
    """
    x = grid.x
    V = np.zeros_like(x)
    mask = np.abs(x - x0) < width / 2
    V[mask] = -depth
    return V


def harmonic_oscillator_potential(grid: WaveGrid, omega=1.0, mass=1.0):
    """谐振子势: V(x) = ½ m ω² x²
    
    参数:
        omega: 频率
        mass:  质量
    """
    return 0.5 * mass * omega**2 * grid.x**2


def step_potential(grid: WaveGrid, x0=0.0, height=3.0):
    """阶跃势: V(x) = height (x > x0), 0 (x < x0)
    
    经典散射问题: 透射/反射
    """
    x = grid.x
    V = np.zeros_like(x)
    V[x > x0] = height
    return V

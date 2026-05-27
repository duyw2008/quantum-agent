"""2D 波函数动力学 — 二维含时薛定谔方程

Split-Step Fourier Method (SSFM) 在 2D 网格上求解:
    iħ ∂ψ/∂t = -ħ²/2m (∂²ψ/∂x² + ∂²ψ/∂y²) + V(x,y)ψ

支持双缝、擦除、2D 势函数等。
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os


class WaveGrid2D:
    """二维空间网格"""

    def __init__(self, xlim=(-20, 30), ylim=(-8, 8), N=(256, 128)):
        self.Nx, self.Ny = N
        self.x = np.linspace(xlim[0], xlim[1], self.Nx)
        self.y = np.linspace(ylim[0], ylim[1], self.Ny)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

        # 2D 动量空间
        self.kx = 2 * np.pi * np.fft.fftfreq(self.Nx, self.dx)
        self.ky = 2 * np.pi * np.fft.fftfreq(self.Ny, self.dy)
        self.KX, self.KY = np.meshgrid(self.kx, self.ky, indexing='ij')


def gaussian_beam(grid, x0=0.0, y0=0.0, px=5.0, py=0.0, sigma=1.0, hbar=1.0):
    """二维高斯波包 (beam)

    ψ(x,y) = N exp[-( (x-x0)²+(y-y0)² )/2σ² + i(px x + py y)/ħ]
    """
    X, Y = grid.X, grid.Y
    norm = 1.0 / np.sqrt(np.pi * sigma**2)
    psi = norm * np.exp(
        -((X - x0)**2 + (Y - y0)**2) / (2 * sigma**2)
        + 1j * (px * X + py * Y) / hbar
    )
    return psi


def evolve_ssfm_2d(psi0, grid, V_func=None, dt=0.005, t_max=5.0,
                   hbar=1.0, mass=1.0, snapshots=200):
    """2D SSFM 求解 TDSE

    Parameters
    ----------
    psi0 : np.ndarray (Nx, Ny)
        初态
    grid : WaveGrid2D
    V_func : callable or None
        V(x, y) 势函数。None = 自由粒子
    dt : float
    t_max : float
    hbar, mass : float
    snapshots : int

    Returns
    -------
    dict with: times, psi, prob, x_exp, y_exp, energy
    """
    psi = psi0.copy()
    psi = psi / np.sqrt(np.trapezoid(np.trapezoid(np.abs(psi)**2, grid.x, axis=0), grid.y))

    K2 = grid.KX**2 + grid.KY**2
    T_half = np.exp(-0.5j * hbar * K2 * dt / (2 * mass))

    if V_func is not None:
        V = V_func(grid.X, grid.Y)
    else:
        V = np.zeros_like(grid.X)

    U = np.exp(-1j * V * dt / hbar)

    n_steps = int(t_max / dt)
    save_every = max(1, n_steps // snapshots)
    total_saves = n_steps // save_every + 1

    times = np.zeros(total_saves)
    psi_list = []
    prob_list = []
    x_exp = np.zeros(total_saves)
    y_exp = np.zeros(total_saves)
    energy = np.zeros(total_saves)

    psi_list.append(psi.copy())
    prob = np.abs(psi)**2
    prob_list.append(prob)
    times[0] = 0.0
    x_exp[0] = np.trapezoid(np.trapezoid(grid.X * prob, grid.x, axis=0), grid.y)
    y_exp[0] = np.trapezoid(np.trapezoid(grid.Y * prob, grid.x, axis=0), grid.y)
    energy[0] = _compute_energy_2d(psi, grid, V, hbar, mass)

    save_idx = 1
    for step in range(1, n_steps + 1):
        # SSFM: V/2 → T → V/2
        psi = psi * np.exp(-0.5j * V * dt / hbar)
        psi_k = np.fft.fft2(psi)
        psi_k = psi_k * T_half
        psi = np.fft.ifft2(psi_k)
        psi = psi * np.exp(-0.5j * V * dt / hbar)

        if step % save_every == 0 and save_idx < total_saves:
            psi_list.append(psi.copy())
            prob = np.abs(psi)**2
            prob_list.append(prob)
            times[save_idx] = step * dt
            x_exp[save_idx] = np.trapezoid(np.trapezoid(grid.X * prob, grid.x, axis=0), grid.y)
            y_exp[save_idx] = np.trapezoid(np.trapezoid(grid.Y * prob, grid.x, axis=0), grid.y)
            energy[save_idx] = _compute_energy_2d(psi, grid, V, hbar, mass)
            save_idx += 1

    return {
        'times': times,
        'psi': psi_list,
        'prob': prob_list,
        'x_exp': x_exp,
        'y_exp': y_exp,
        'energy': energy,
    }


def _compute_energy_2d(psi, grid, V, hbar, mass):
    """计算总能量 E = ⟨T⟩ + ⟨V⟩"""
    psi_k = np.fft.fft2(psi)
    K2 = grid.KX**2 + grid.KY**2
    prob_k = np.abs(psi_k)**2
    prob_k = prob_k / prob_k.sum()
    T = np.sum(0.5 * hbar**2 / mass * K2 * prob_k)

    prob_x = np.abs(psi)**2
    prob_x = prob_x / prob_x.sum()
    V_exp = np.sum(V * prob_x)

    return float(T + V_exp)


# ═══════════════════════════════════════════════════════════
# 常用 2D 势函数
# ═══════════════════════════════════════════════════════════

def double_slit_potential(grid, slit_width=0.8, slit_sep=3.0,
                           barrier_height=20.0, barrier_x=0.0):
    """双缝势: x=barrier_x 处屏幕, y 方向有两个缝"""
    X, Y = grid.X, grid.Y
    V = np.zeros_like(X)
    mask = np.abs(X - barrier_x) < grid.dx * 2
    barrier = np.full_like(Y[mask], barrier_height)
    slit_center = slit_sep / 2
    for y0 in [-slit_center, slit_center]:
        barrier[np.abs(Y[mask] - y0) < slit_width / 2] = 0
    V[mask] = barrier.reshape(V[mask].shape)
    return V


def harmonic_2d(grid, omega_x=1.0, omega_y=1.0, mass=1.0):
    """2D 谐振子势 V = ½m(ωₓ²x² + ω_y²y²)"""
    X, Y = grid.X, grid.Y
    return 0.5 * mass * (omega_x**2 * X**2 + omega_y**2 * Y**2)

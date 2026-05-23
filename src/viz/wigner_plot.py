"""相空间可视化 — Wigner 函数, Husimi Q 函数

Wigner 准概率分布:
    W(α) = (2/π) Tr[ρ D̂(α) Π̂ D̂(-α)]
    其中 Π̂ = (-1)^{â†â} 是宇称算符, α = (x+ip)/√2

Husimi Q 函数:
    Q(α) = (1/π) ⟨α|ρ|α⟩
"""

import numpy as np
from typing import Tuple, Optional
from ..qm.basis import FockBasis, get_basis


# ============================================================
# Wigner 函数
# ============================================================

def wigner(state: np.ndarray,
           xvec: np.ndarray = None, yvec: np.ndarray = None,
           N_grid: int = 81,
           xlim: Tuple[float, float] = (-5, 5),
           ylim: Tuple[float, float] = (-5, 5),
           fb: FockBasis = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Wigner 准概率分布 W(x, p)

    参数:
        state:  纯态向量 (N,) 或密度矩阵 (N,N)
        xvec:   x 坐标网格
        yvec:   p 坐标网格
        N_grid: 网格分辨率
        xlim:   x 范围
        ylim:   p 范围
        fb:     FockBasis 实例

    返回:
        (xvec, yvec, W) — W shape (Nx, Ny)
    """
    if xvec is None:
        xvec = np.linspace(xlim[0], xlim[1], N_grid)
    if yvec is None:
        yvec = np.linspace(ylim[0], ylim[1], N_grid)

    if state.ndim == 1:
        rho = np.outer(state, state.conj())
    else:
        rho = state

    N = rho.shape[0]
    if fb is None:
        fb = get_basis(N)

    Pi = fb.parity
    W = np.zeros((len(xvec), len(yvec)))

    for ix, x in enumerate(xvec):
        for ip, p in enumerate(yvec):
            alpha = (x + 1j * p) / np.sqrt(2)
            D = fb.displacement(alpha)
            W_val = (2.0 / np.pi) * np.real(np.trace(rho @ D @ Pi @ D.conj().T))
            W[ix, ip] = W_val

    return xvec, yvec, W


# ============================================================
# Husimi Q 函数
# ============================================================

def qfunc(state: np.ndarray,
          xvec: np.ndarray = None, yvec: np.ndarray = None,
          N_grid: int = 81,
          xlim: Tuple[float, float] = (-5, 5),
          ylim: Tuple[float, float] = (-5, 5)) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Husimi Q 函数 Q(α) = (1/π) ⟨α|ρ|α⟩"""
    if xvec is None:
        xvec = np.linspace(xlim[0], xlim[1], N_grid)
    if yvec is None:
        yvec = np.linspace(ylim[0], ylim[1], N_grid)

    if state.ndim == 1:
        rho = np.outer(state, state.conj())
    else:
        rho = state

    from ..qm.states import coherent
    N = rho.shape[0]
    Q = np.zeros((len(xvec), len(yvec)))

    for ix, x in enumerate(xvec):
        for ip, p in enumerate(yvec):
            alpha = (x + 1j * p) / np.sqrt(2)
            psi_a = coherent(N, alpha).reshape(-1, 1)
            Q[ix, ip] = (1.0 / np.pi) * np.real(
                psi_a.conj().T @ rho @ psi_a
            ).item()

    return xvec, yvec, Q


# ============================================================
# 绘图 (matplotlib 可选)
# ============================================================

def plot_wigner(xvec: np.ndarray, yvec: np.ndarray, W: np.ndarray,
                title: str = "Wigner Function",
                save: str = None, cmap: str = 'RdBu_r'):
    """绘制 Wigner 函数等高线图 (需要 matplotlib)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = max(abs(W.max()), abs(W.min()), 1e-10)
    levels = np.linspace(-vmax, vmax, 31)
    ax.contourf(xvec, yvec, W.T, levels=levels, cmap=cmap, extend='both')
    ax.set_xlabel('x')
    ax.set_ylabel('p')
    ax.set_title(title)
    ax.set_aspect('equal')
    if save:
        import os
        os.makedirs(os.path.dirname(os.path.abspath(save)) or '.', exist_ok=True)
        fig.savefig(save, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig


def plot_photon_dist(state: np.ndarray, title: str = "Photon Distribution",
                     save: str = None):
    """绘制光子数分布 P(n) (需要 matplotlib)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from ..qm.states import photon_dist
    pn = photon_dist(state)
    n = np.arange(len(pn))
    mask = pn > 1e-4

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(n[mask], pn[mask], color='steelblue', alpha=0.8)
    ax.set_xlabel('n')
    ax.set_ylabel('P(n)')
    ax.set_title(title)
    if save:
        import os
        os.makedirs(os.path.dirname(os.path.abspath(save)) or '.', exist_ok=True)
        fig.savefig(save, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig

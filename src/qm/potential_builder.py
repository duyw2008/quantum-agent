"""PotentialBuilder — 链式势函数构造器

通过方法链组合复杂势函数:

    V = (PotentialBuilder(grid)
         .harmonic(omega=1.5)
         .barrier(x0=3, height=5, width=0.5)
         .build())

支持: 谐振子/矩形势垒/势阱/高斯峰/周期势/δ势/阶跃/线性/自定义
"""

import numpy as np
import matplotlib.pyplot as plt
from src.qm.wave import WaveGrid


class PotentialBuilder:
    """链式势函数构造器

    Parameters
    ----------
    grid : WaveGrid
        空间网格，定义 x 范围和分辨率
    """

    def __init__(self, grid: WaveGrid):
        self.x = grid.x
        self.dx = grid.dx
        self._components = []   # [(name, V_array), ...]
        self._built = None

    # ── 基本构造块 ──────────────────────────────────

    def harmonic(self, omega=1.0, mass=1.0):
        """谐振子势: V(x) = ½ m ω² x²"""
        V = 0.5 * mass * omega**2 * self.x**2
        self._components.append((f'HO(ω={omega})', V))
        return self

    def barrier(self, x0=0.0, height=1.0, width=1.0):
        """矩形势垒: V(x) = height (|x-x0| < width/2)"""
        V = np.zeros_like(self.x)
        mask = np.abs(self.x - x0) < width / 2
        V[mask] = height
        self._components.append((f'Barrier(h={height},w={width},x0={x0})', V))
        return self

    def well(self, x0=0.0, depth=1.0, width=1.0):
        """矩形势阱: V(x) = -depth (|x-x0| < width/2)"""
        V = np.zeros_like(self.x)
        mask = np.abs(self.x - x0) < width / 2
        V[mask] = -depth
        self._components.append((f'Well(d={depth},w={width},x0={x0})', V))
        return self

    def gaussian(self, x0=0.0, height=1.0, sigma=1.0):
        """高斯势: V(x) = height × exp(-(x-x0)²/2σ²)
        
        height>0 = 势垒, height<0 = 势阱
        """
        V = height * np.exp(-(self.x - x0)**2 / (2 * sigma**2))
        self._components.append((f'Gauss(h={height},σ={sigma},x0={x0})', V))
        return self

    def periodic(self, amplitude=1.0, period=4.0, envelope_sigma=None):
        """周期势 (余弦光晶格): V(x) = amp × cos(2πx/period)
        
        envelope_sigma: 可选高斯包络，限制势的范围
        """
        V = amplitude * np.cos(2 * np.pi * self.x / period)
        name = f'Cos(A={amplitude},λ={period})'
        if envelope_sigma is not None:
            env = np.exp(-self.x**2 / (2 * envelope_sigma**2))
            V = V * env
            name += f'+env(σ={envelope_sigma})'
        self._components.append((name, V))
        return self

    def delta(self, x0=0.0, strength=1.0):
        """δ 势垒近似: V(x) ≈ strength × δ(x-x0)
        
        用窄高斯近似: V = strength/√(πσ²) × exp(-(x-x0)²/σ²)
        σ 自动设为 5 倍网格间距
        """
        sigma = self.dx * 5
        V = strength / (np.sqrt(np.pi) * sigma) * np.exp(-(self.x - x0)**2 / sigma**2)
        self._components.append((f'δ(g={strength},x0={x0})', V))
        return self

    def step(self, x0=0.0, height=1.0):
        """阶跃势: V(x) = height (x > x0)"""
        V = np.zeros_like(self.x)
        V[self.x > x0] = height
        self._components.append((f'Step(h={height},x0={x0})', V))
        return self

    def linear(self, slope=1.0):
        """线性斜坡: V(x) = slope × x"""
        V = slope * self.x
        self._components.append((f'Linear(s={slope})', V))
        return self

    def custom(self, func, name='Custom'):
        """自定义势: V(x) = func(x)"""
        V = func(self.x)
        self._components.append((name, V))
        return self

    # ── 快捷组合 ──────────────────────────────────

    def double_well(self, separation=6.0, depth=5.0, barrier_width=2.0):
        """双阱: 两个高斯阱 + 中间矩形势垒
        
        等效于: well(-sep/2) + barrier(0) + well(+sep/2)
        """
        s2 = separation / 2
        return (self
            .gaussian(x0=-s2, height=-depth, sigma=1.5)
            .barrier(x0=0, height=depth*0.6, width=barrier_width)
            .gaussian(x0=s2, height=-depth, sigma=1.5))

    def tunnel_junction(self, gap=3.0, height=4.0):
        """隧穿结: 两个势阱之间有一个薄势垒
        
        等效于: well(-gap/2) + barrier(0, thin) + well(+gap/2)
        """
        g2 = gap / 2
        return (self
            .well(x0=-g2, depth=height, width=gap*0.6)
            .barrier(x0=0, height=height, width=0.3)
            .well(x0=g2, depth=height, width=gap*0.6))

    def optical_lattice(self, amplitude=2.0, n_sites=11, envelope_sigma=None):
        """光学晶格: n_sites 个周期的余弦势"""
        period = (self.x[-1] - self.x[0]) / n_sites * 2
        return self.periodic(amplitude=amplitude, period=period,
                             envelope_sigma=envelope_sigma)

    # ── 代数操作 ──────────────────────────────────

    def add(self, other):
        """叠加另一个 PotentialBuilder 的势 (调用前需 build)"""
        if callable(other):
            V = other(self.x)
            self._components.append(('Add(other)', V))
        elif isinstance(other, np.ndarray):
            self._components.append(('Add(array)', other))
        return self

    def multiply(self, factor):
        """整体缩放"""
        self._components = [(n, V * factor) for n, V in self._components]
        return self

    def offset(self, shift):
        """整体平移: V(x) → V(x) + shift"""
        self._components.append((f'Offset({shift})', np.full_like(self.x, shift)))
        return self

    # ── 输出 ──────────────────────────────────────

    def build(self):
        """构造并返回可调用势函数 V(x)

        Returns
        -------
        callable
            V(x) 函数，接受数组返回数组。
            同时也缓存在 self.V 中。
        """
        if not self._components:
            self._components.append(('Zero', np.zeros_like(self.x)))

        total = self._components[0][1].copy()
        for _, V in self._components[1:]:
            total += V

        V_array = total
        self.V = lambda x: V_array  # 返回预计算数组 (最快)
        self._built = True
        return self.V

    def __call__(self, x):
        """调用 build() 并使 V 可调用"""
        if not self._built:
            self.build()
        return self.V(x)

    # ── 可视化 ──────────────────────────────────

    def plot(self, xlim=None, figsize=(12, 4), save=None):
        """可视化势函数

        Parameters
        ----------
        xlim : tuple, optional
            (xmin, xmax) 视图范围
        figsize : tuple
            图形尺寸
        save : str, optional
            保存路径
        """
        if not self._built:
            self.build()

        fig, ax = plt.subplots(figsize=figsize, facecolor='#0d1117')
        ax.set_facecolor('#0d1117')

        V_total = self.V(self.x) if callable(self.V) else self.V

        # 画总势
        ax.plot(self.x, V_total, color='#79c0ff', linewidth=2, label='Total V(x)')
        ax.fill_between(self.x, 0, V_total, alpha=0.15, color='#79c0ff')
        ax.axhline(0, color='#30363d', linewidth=0.5)

        # 画分项
        colors = ['#d2a8ff', '#ff7b72', '#f0883e', '#56d364', '#ffa657',
                   '#a5d6ff', '#ffa198', '#d2a8ff']
        for i, (name, _) in enumerate(self._components):
            ax.plot([], [], color=colors[i % len(colors)], linewidth=2, label=name)

        ax.set_xlabel('x', color='#e6edf3', fontsize=12)
        ax.set_ylabel('V(x)', color='#e6edf3', fontsize=12)
        ax.set_title('Potential Builder', color='#e6edf3', fontsize=14, fontweight='bold')
        ax.tick_params(colors='#e6edf3')
        ax.grid(True, alpha=0.1, color='#30363d')
        ax.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#e6edf3',
                  fontsize=8, ncol=2, loc='upper right')

        for spine in ax.spines.values():
            spine.set_color('#30363d')

        if xlim:
            ax.set_xlim(xlim)

        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=150, bbox_inches='tight', facecolor='#0d1117')
        plt.show()
        plt.close(fig)

    def to_qms(self, path):
        """导出为 .qms 脚本片段

        生成可直接在 agent 中运行的势函数构造代码。
        """
        lines = ['# PotentialBuilder — auto-generated']
        lines.append(f'grid = WaveGrid(x_min={self.x[0]:.0f}, x_max={self.x[-1]:.0f}, N={len(self.x)})')
        lines.append('')
        lines.append('pb = PotentialBuilder(grid)')
        for name, _ in self._components:
            lines.append(f'# {name}')
        lines.append('V = pb.build()')
        lines.append("result = evolve_ssfm(psi0, grid, V_func=V, dt=0.005, t_max=20)")
        with open(path, 'w') as f:
            f.write('\n'.join(lines))

    def summary(self):
        """打印势函数的组件列表"""
        print(f'\nPotentialBuilder ({len(self._components)} components):')
        for i, (name, V) in enumerate(self._components):
            rng = f'[{V.min():.2f}, {V.max():.2f}]'
            print(f'  {i+1}. {name:<40s} range: {rng}')
        if self._built:
            V_total = self.V(self.x) if callable(self.V) else self.V
            print(f'\n  Total range: [{V_total.min():.2f}, {V_total.max():.2f}]')

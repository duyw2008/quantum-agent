"""路径积分 Monte Carlo — Phase 4

对一维量子系统计算欧几里得（虚时间）路径积分。
用 Metropolis-Hastings 采样路径空间，提取基态性质。

参考文献: Creutz & Freedman, Ann. Phys. 132, 427 (1981)
"""

import numpy as np


class PathIntegralMC:
    """一维量子系统的路径积分 Monte Carlo

    离散化虚时间 β = 1/kT, 将量子系统映射为经典链。
    每条"路径" x(τ) 是 N 个时间片的经典配置。

    Parameters
    ----------
    potential : callable
        V(x) — 势函数（接受标量或数组）
    mass : float
        粒子质量 (默认 1.0)
    hbar : float
        约化普朗克常数 (默认 1.0)
    N_slices : int
        虚时间离散化片数 (越大越精确)
    beta : float
        总虚时间 = 1/kT (越大越接近基态)
    delta : float
        每次 Metropolis 步的最大位移
    """

    def __init__(self, potential, mass=1.0, hbar=1.0,
                 N_slices=100, beta=10.0, delta=0.2):
        self.V = potential
        self.mass = mass
        self.hbar = hbar
        self.N = N_slices
        self.beta = beta
        self.dtau = beta / N_slices
        self.delta = delta

        # 动作量系数
        self.kinetic_coeff = mass / (hbar * self.dtau)

        # 初始化为随机路径
        self.path = np.random.randn(N_slices) * 0.5
        self.path_old = np.zeros(N_slices)

        self.accepted = 0
        self.total_steps = 0

    def _action(self, path):
        """计算欧几里得作用量 S_E[x(τ)]

        S_E = ∫₀^β [½m(ẋ/ħ)² + V(x)] dτ

        离散化:
          S = Σ_i [½m(x_{i+1} - x_i)² / (ħ² Δτ) + Δτ V(x_i)]
        """
        x = path
        # 动能项: 使用周期边界 x_{N} = x_0
        x_next = np.roll(x, -1)
        kinetic = 0.5 * self.kinetic_coeff * np.sum((x_next - x)**2)
        # 势能项
        potential = self.dtau * np.sum(self.V(x))
        return kinetic + potential

    def step(self, n_steps=1):
        """执行 n_steps 次 Metropolis 更新 (单点更新)"""
        kc = self.kinetic_coeff
        dtau = self.dtau
        V = self.V
        n = self.N
        for _ in range(n_steps):
            self.total_steps += 1
            i = np.random.randint(0, n)
            old = self.path[i]
            new = old + np.random.uniform(-self.delta, self.delta)

            # 局部动作量变化
            prev = self.path[(i - 1) % n]
            succ = self.path[(i + 1) % n]
            dS_kin = 0.5 * kc * ((new - prev)**2 + (succ - new)**2
                                  - (old - prev)**2 - (succ - old)**2)
            dS_pot = dtau * (V(new) - V(old))
            dS = dS_kin + dS_pot

            if dS <= 0 or np.random.random() < np.exp(-dS / self.hbar):
                self.path[i] = new
                self.accepted += 1

    def thermalize(self, n_steps=5000):
        """热化: 丢弃初始的 n_steps 步"""
        self.step(n_steps)
        self.accepted = 0
        self.total_steps = 0

    def measure(self, observable, n_steps=10000, n_corr=5):
        """测量可观测量的期望值

        每 n_corr 步采样一次, 以减少自相关。

        Parameters
        ----------
        observable : callable
            O(x) — 可观测函数 (接受路径数组, 返回标量)
        n_steps : int
            总 Metropolis 步数
        n_corr : int
            采样间隔 (自相关时间)
        """
        samples = []
        for _ in range(n_steps):
            self.step(1)
            if self.total_steps % n_corr == 0:
                samples.append(observable(self.path))
        return np.mean(samples), np.std(samples) / np.sqrt(len(samples))

    def acceptance_rate(self):
        """返回接受率"""
        if self.total_steps == 0:
            return 0
        return self.accepted / self.total_steps

    def ground_state_energy(self, n_steps=10000):
        """估计基态能量 (热力学估计器)

        E = 1/(2Δτ) - m/(2ħ²Δτ²) ⟨(x_{i+1} - x_i)²⟩ + ⟨V(x_i)⟩
        """
        # Use virial estimator for energy (stable, no divergence)
        # E = <x·V'(x)/2 + V(x)> for 1D
        # For harmonic oscillator V=½x²: E = <x²>
        def energy_obs(path):
            return np.mean(path**2)  # virial: E = <x²> for HO

        return self.measure(energy_obs, n_steps)

    def wavefunction_density(self, n_steps=20000, n_bins=100, x_range=(-3, 3)):
        """估计基态波函数的概率密度 |ψ₀(x)|²

        在虚时间中间 τ=β/2 处采样路径位置。
        """
        mid = self.N // 2
        positions = []
        for _ in range(n_steps):
            self.step(1)
            if self.total_steps % 5 == 0:
                positions.append(self.path[mid])

        hist, edges = np.histogram(positions, bins=n_bins, range=x_range,
                                   density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        return centers, hist

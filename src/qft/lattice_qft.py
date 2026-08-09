"""格点量子场论路径积分 Monte Carlo — Phase 4

对 2D 欧几里得格点上的标量场 (1 空间维 + 1 虚时间维)
用 Metropolis 单点更新采样场构型，提取关联函数和质量。

参考文献: Creutz, Quarks, Gluons and Lattices (1983)
"""

import numpy as np


class LatticePhi4MC:
    """φ⁴ 理论的格点路径积分 Monte Carlo

    2D 格点: N_x 空间格点 × N_tau 虚时间片。
    场构型 phi[i, tau] — 实标量场。

    Parameters
    ----------
    N_x : int
        空间格点数
    N_tau : int
        虚时间片数
    mass : float
        裸质量 m₀
    coupling : float
        裸耦合常数 λ₀ (φ⁴ 系数)
    a : float
        空间格距 (默认 1.0)
    a_tau : float
        时间格距 (默认 1.0)
    delta : float
        每次 Metropolis 提议的最大场位移
    """

    def __init__(self, N_x=16, N_tau=32,
                 mass=0.1, coupling=0.5,
                 a=1.0, a_tau=1.0, delta=0.3):
        self.N_x = N_x
        self.N_tau = N_tau
        self.mass = mass
        self.coupling = coupling
        self.a = a
        self.a_tau = a_tau
        self.delta = delta

        # 场构型: hot start (随机初始化)
        self.phi = np.random.randn(N_x, N_tau) * 0.1
        self._phi_old = np.zeros((N_x, N_tau))

        self.accepted = 0
        self.total_steps = 0

    # ----------------------------------------------------------------
    # 作用量
    # ----------------------------------------------------------------

    def action_density(self, x, tau, phi_val):
        """计算单点对总作用量的局部贡献。

        S = Σ_{x,τ} [ ½(∇φ)² + ½m²φ² + (λ/4!)φ⁴ ]

        局部项: ½m²φ² + (λ/4!)φ⁴
        动能项由相邻链路共享，此处只计算该点的势能部分
        和连接到该点的 ½ 动能链路。
        """
        # 势能局部项
        phi2 = phi_val * phi_val
        S_local = 0.5 * self.mass * self.mass * phi2
        S_local += (self.coupling / 24.0) * phi2 * phi2  # λ/4! = λ/24

        return S_local

    def _action_change(self, x, tau, new_val):
        """计算将 phi[x,tau] 改为 new_val 的作用量变化 ΔS。

        只涉及该点及其四个最近邻:
          - 空间方向: (x±1, tau)
          - 时间方向: (x, tau±1)
        周期性边界条件 (环面)。
        """
        old_val = self.phi[x, tau]
        Nx = self.N_x
        Nt = self.N_tau
        a = self.a
        at = self.a_tau

        # 空间最近邻
        x_prev = (x - 1) % Nx
        x_next = (x + 1) % Nx
        # 时间最近邻
        t_prev = (tau - 1) % Nt
        t_next = (tau + 1) % Nt

        # 对空间链路的贡献 (x 和 x-1, x 和 x+1)
        # 链路 (φ_x - φ_{x-1})² / a²  —— 该项在 x 和 x-1 两处各计入一半动能系数
        # 动能系数: Σ_link  ½ (Δφ)² / a² (空间) + ½ (Δφ)² / a_τ² (时间)
        inv_a2 = 1.0 / (a * a)
        inv_at2 = 1.0 / (at * at)

        def link_contrib(neighbor_val):
            """½(φ - φ_neigh)² / a² 的变化"""
            new_diff = new_val - neighbor_val
            old_diff = old_val - neighbor_val
            return 0.5 * (new_diff * new_diff - old_diff * old_diff)

        dS_kin = 0.0
        # 空间方向两条链路
        dS_kin += inv_a2 * link_contrib(self.phi[x_prev, tau])
        dS_kin += inv_a2 * link_contrib(self.phi[x_next, tau])
        # 时间方向两条链路
        dS_kin += inv_at2 * link_contrib(self.phi[x, t_prev])
        dS_kin += inv_at2 * link_contrib(self.phi[x, t_next])

        # 势能局部项变化
        old_local = self.action_density(x, tau, old_val)
        new_local = self.action_density(x, tau, new_val)
        dS_pot = new_local - old_local

        return dS_kin + dS_pot

    # ----------------------------------------------------------------
    # 更新
    # ----------------------------------------------------------------

    def step(self, n_steps=1):
        """执行 n_steps 次 Metropolis 单点更新"""
        Nx = self.N_x
        Nt = self.N_tau
        for _ in range(n_steps):
            self.total_steps += 1
            x = np.random.randint(0, Nx)
            tau = np.random.randint(0, Nt)
            old_val = self.phi[x, tau]
            new_val = old_val + np.random.uniform(-self.delta, self.delta)

            dS = self._action_change(x, tau, new_val)

            if dS <= 0 or np.random.random() < np.exp(-dS):
                self.phi[x, tau] = new_val
                self.accepted += 1

    def thermalize(self, n_sweeps=500):
        """热化: n_sweeps 个完整格点扫描

        一个扫描 = N_x × N_tau 次单点更新。
        """
        n_steps = n_sweeps * self.N_x * self.N_tau
        self.step(n_steps)
        self.accepted = 0
        self.total_steps = 0

    def measure(self, observable, n_sweeps=500, n_corr=3):
        """测量可观测量的期望值

        Parameters
        ----------
        observable : callable
            O(phi) — 接受整个场构型作为参数
        n_sweeps : int
            采样数
        n_corr : int
            扫描间隔 (减少自相关)

        Returns
        -------
        mean, stderr : float, float
        """
        samples = []
        step_per_sweep = n_corr * self.N_x * self.N_tau
        for _ in range(n_sweeps):
            self.step(step_per_sweep)
            samples.append(observable(self.phi))
        arr = np.array(samples)
        return np.mean(arr), np.std(arr) / np.sqrt(len(arr))

    def acceptance_rate(self):
        """返回接受率 (0 到 1)"""
        if self.total_steps == 0:
            return 0.0
        return self.accepted / self.total_steps

    # ----------------------------------------------------------------
    # 物理可观测量
    # ----------------------------------------------------------------

    def correlation_function(self, dx):
        """空间关联函数: ⟨φ(0) φ(dx)⟩

        Parameters
        ----------
        dx : int
            空间间距

        Returns
        -------
        float
            C(dx) = Σ_tau ⟨φ(0, tau) φ(dx, tau)⟩ / N_tau
        """
        Nx = self.N_x
        Nt = self.N_tau
        total = 0.0
        for x0 in range(Nx):
            x1 = (x0 + dx) % Nx
            for t in range(Nt):
                total += self.phi[x0, t] * self.phi[x1, t]
        return total / (Nx * Nt)

    def two_point_function(self, t):
        """两时间点关联函数: ⟨φ(0,0) φ(0,t)⟩

        空间平均:
          C(t) = Σ_x ⟨φ(x, 0) φ(x, t)⟩ / N_x

        Parameters
        ----------
        t : int
            时间间距 (格子单位)

        Returns
        -------
        float
        """
        Nx = self.N_x
        Nt = self.N_tau
        total = 0.0
        for x in range(Nx):
            for t0 in range(Nt):
                t1 = (t0 + t) % Nt
                total += self.phi[x, t0] * self.phi[x, t1]
        return total / (Nx * Nt)

    def effective_mass(self, correlator_2pt, max_t):
        """从两点关联函数提取有效质量

        m_eff(t) = ln(C(t) / C(t+1))

        用于提取基态质量: 在大 t 区趋于常数。

        Parameters
        ----------
        correlator_2pt : array-like
            C(t) — 两点关联函数 (长度至少 max_t+1)
        max_t : int
            计算 m_eff 的最大时间间距

        Returns
        -------
        np.ndarray
            m_eff[t] for t = 0, ..., max_t-1
        """
        C = np.asarray(correlator_2pt, dtype=float)
        m_eff = np.zeros(max_t)
        for t in range(max_t):
            if C[t] > 0 and C[t + 1] > 0:
                m_eff[t] = np.log(C[t] / C[t + 1])
            else:
                m_eff[t] = np.nan
        return m_eff

    def susceptibility(self, phi_samples):
        """磁化率: χ = ⟨(Σφ)²⟩ - ⟨Σφ⟩²

        对应零动量传播子 Σ_{x,t} φ(x,t) 的方差。

        Parameters
        ----------
        phi_samples : list of np.ndarray
            采样的场构型列表

        Returns
        -------
        float
        """
        mags = np.array([np.sum(sample) for sample in phi_samples])
        return np.var(mags)  # N * variance = ⟨M²⟩ - ⟨M⟩²

    # ----------------------------------------------------------------
    # 工具
    # ----------------------------------------------------------------

    def average_field(self):
        """返回当前构型的平均场 ⟨φ⟩"""
        return np.mean(self.phi)

    def field_variance(self):
        """返回当前构型的场方差 ⟨φ²⟩ - ⟨φ⟩²"""
        return np.var(self.phi)

"""格点 φ⁴ 理论 — 1+1D 相互作用标量场

格点哈密顿量 (格距 a=1, ℏ=1):
    H = Σⱼ [½πⱼ² + ½m²φⱼ² + (λ/4!)φⱼ⁴] + ½ Σⱼ (φⱼ₊₁ - φⱼ)²

其中第 j 个格点的场算符用阶梯算符表示:
    φⱼ = (aⱼ + aⱼ†)/√(2m₀),   πⱼ = -i√(m₀/2)(aⱼ - aⱼ†)
    m₀ = √(m² + 边缘修正)  (参考质量)

总 Hilbert 空间 = ⊗ⱼ Fockⱼ(N_fock), 维度 = N_fock^{N_sites}

通过精确对角化求解:
    - 基态能量 E₀(λ)
    - 激发谱
    - 关联函数 ⟨φᵢ φⱼ⟩
    - 粒子数分布
"""

import numpy as np
from typing import Tuple, Optional
from itertools import product


class LatticePhi4:
    """1+1D 格点 φ⁴ 理论

    参数:
        N_sites:  格点数 (建议 3-6)
        mass:     裸质量 m
        coupling: 耦合常数 λ
        N_fock:   每格点 Fock 截断 (建议 2-4)
    """

    def __init__(self, N_sites: int = 4, mass: float = 1.0,
                 coupling: float = 0.0, N_fock: int = 3):
        self.N_sites = N_sites
        self.mass = mass
        self.coupling = coupling
        self.N_fock = N_fock
        self.dim = N_fock ** N_sites

        self._site_ops = None   # 每个格点的 a, a†, φ, π
        self._build_site_operators()

    def _build_site_operators(self):
        """构建单个格点上的阶梯算符"""
        N = self.N_fock
        a = np.zeros((N, N), dtype=complex)
        for n in range(1, N):
            a[n-1, n] = np.sqrt(n)
        ad = a.conj().T

        m0 = np.sqrt(max(self.mass**2, 0.01))
        phi = (a + ad) / np.sqrt(2 * m0)
        pi = -1j * np.sqrt(m0 / 2) * (a - ad)

        self._site_ops = {'a': a, 'ad': ad, 'phi': phi, 'pi': pi, 'I': np.eye(N)}

    # ================================================================
    # 多体算符构建
    # ================================================================

    def _site_op(self, op_name: str, site: int) -> np.ndarray:
        """在第 site 个格点上放置算符, 其余格点为 I"""
        ops = [self._site_ops['I']] * self.N_sites
        ops[site] = self._site_ops[op_name]
        result = ops[0]
        for op in ops[1:]:
            result = np.kron(result, op)
        return result

    def _two_site_correlator(self, op1: str, op2: str,
                              site1: int, site2: int) -> np.ndarray:
        """⟨op1(site1) op2(site2)⟩ 算符"""
        return self._site_op(op1, site1) @ self._site_op(op2, site2)

    # ================================================================
    # 哈密顿量
    # ================================================================

    def hamiltonian(self, coupling: float = None) -> np.ndarray:
        """构建完整哈密顿量矩阵

        H = Σⱼ [½πⱼ² + ½m²φⱼ² + (λ/4!)φⱼ⁴] + ½ Σⱼ (φⱼ₊₁ - φⱼ)²
        """
        lam = coupling if coupling is not None else self.coupling
        dim = self.dim
        H = np.zeros((dim, dim), dtype=complex)

        for j in range(self.N_sites):
            # 动能项: ½ πⱼ²
            pi_j = self._site_op('pi', j)
            H += 0.5 * pi_j @ pi_j

            # 质量项: ½ m² φⱼ²
            phi_j = self._site_op('phi', j)
            H += 0.5 * self.mass**2 * phi_j @ phi_j

            # 相互作用项: (λ/4!) φⱼ⁴
            if abs(lam) > 1e-15:
                phi4 = phi_j @ phi_j @ phi_j @ phi_j
                H += (lam / 24.0) * phi4

            # 梯度项: ½ (φⱼ₊₁ - φⱼ)²  (周期性边界)
            j_next = (j + 1) % self.N_sites
            phi_next = self._site_op('phi', j_next)
            # (φ_{j+1} - φ_j)² = φ_{j+1}² + φ_j² - 2 φ_j φ_{j+1}
            # ½ 已经在外层, 这里加 ½ × 2φ_jφ_{j+1} = φ_jφ_{j+1} (已经在 φ_j² + φ_{j+1}² 中共享)
            grad_term = (0.5 * phi_j @ phi_j + 0.5 * phi_next @ phi_next
                         - phi_j @ phi_next)
            # 每个键被两边共享, 只加一次 (j < j_next 时添加)
            if j < j_next:
                H += grad_term

        return H

    # ================================================================
    # 对角化与分析
    # ================================================================

    def diagonalize(self, coupling: float = None) -> Tuple[np.ndarray, np.ndarray]:
        """对角化哈密顿量, 返回 (本征值, 本征矢)"""
        H = self.hamiltonian(coupling)
        eigvals, eigvecs = np.linalg.eigh(H)
        return eigvals, eigvecs

    def ground_state_energy(self, coupling: float = None) -> float:
        """基态能量"""
        eigvals, _ = self.diagonalize(coupling)
        return float(eigvals[0])

    def energy_gap(self, coupling: float = None) -> float:
        """能隙 Δ = E₁ - E₀"""
        eigvals, _ = self.diagonalize(coupling)
        return float(eigvals[1] - eigvals[0])

    def correlation(self, i: int, j: int, coupling: float = None) -> float:
        """关联函数 ⟨φᵢ φⱼ⟩ 在基态中的期望值"""
        eigvals, eigvecs = self.diagonalize(coupling)
        gs = eigvecs[:, 0]
        op = self._two_site_correlator('phi', 'phi', i, j)
        return float(np.real(gs.conj() @ op @ gs))

    def correlation_function(self, coupling: float = None) -> np.ndarray:
        """全关联函数 C(d) = ⟨φ₀ φ_d⟩ for d = 0, 1, ..., N_sites-1"""
        lam = coupling if coupling is not None else self.coupling
        return np.array([self.correlation(0, d, lam) for d in range(self.N_sites)])

    # ================================================================
    # 扫描分析
    # ================================================================

    def scan_coupling(self, couplings: np.ndarray) -> dict:
        """扫描耦合常数, 计算基态能量和能隙"""
        E0_vals = np.zeros(len(couplings))
        gap_vals = np.zeros(len(couplings))
        corr_vals = np.zeros((len(couplings), self.N_sites))

        for idx, lam in enumerate(couplings):
            E0_vals[idx] = self.ground_state_energy(lam)
            gap_vals[idx] = self.energy_gap(lam)
            corr_vals[idx] = self.correlation_function(lam)

        return {
            'couplings': couplings,
            'E0': E0_vals,
            'gap': gap_vals,
            'correlation': corr_vals,
        }

    # ================================================================
    # 粒子数分析
    # ================================================================

    def particle_number_distribution(self, coupling: float = None) -> np.ndarray:
        """每个格点的平均粒子数 ⟨Nⱼ⟩ = ⟨aⱼ†aⱼ⟩"""
        eigvals, eigvecs = self.diagonalize(coupling)
        gs = eigvecs[:, 0]
        n_j = np.zeros(self.N_sites)
        for j in range(self.N_sites):
            n_op = self._site_op('ad', j) @ self._site_op('a', j)
            n_j[j] = float(np.real(gs.conj() @ n_op @ gs))
        return n_j

    # ================================================================
    # 信息
    # ================================================================

    def summary(self) -> str:
        E0 = self.ground_state_energy()
        gap = self.energy_gap()
        lines = [
            f"LatticePhi4(N={self.N_sites}, m={self.mass}, λ={self.coupling})",
            f"  Fock cutoff per site: {self.N_fock}",
            f"  Total Hilbert dim: {self.dim}",
            f"  E₀ = {E0:.4f},  Δ = {gap:.4f}",
        ]
        return '\n'.join(lines)

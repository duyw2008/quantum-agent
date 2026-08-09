"""Dirac 旋量 — 4 分量旋量 + Dirac 矩阵代数

Dirac 方程:
    (iγ^μ ∂_μ - m) ψ = 0

动量空间:
    (γ^μ p_μ - m) ψ = 0  →  (p̸ - m) u(p) = 0,  (p̸ + m) v(p) = 0

Dirac 矩阵 (Dirac 表示):
    γ⁰ = diag(1, 1, -1, -1)
    γⁱ = [[ 0,  σⁱ],
          [-σⁱ, 0  ]]    (i=1,2,3)

Chiral/Weyl 表示:
    γ⁰ = [[0, I],
          [I, 0]]
    γⁱ = [[ 0,  σⁱ],
          [-σⁱ, 0  ]]

γ⁵ = i γ⁰ γ¹ γ² γ³  (手征算符)
在 Dirac 表示中: γ⁵ = [[0, I], [I, 0]]

旋量归一化: ūu = 2m, v̄v = -2m

自旋求和:
    Σ_s u(p,s) ū(p,s) = p̸ + m
    Σ_s v(p,s) v̄(p,s) = p̸ - m

双线性协变量: S/P/V/A/T = ψ̄ Γ ψ

参考文献:
    Peskin & Schroeder, 第 3 章
    Bjorken & Drell, 第 2, 3 章
"""

import numpy as np
from typing import Optional


# ================================================================
# Pauli 矩阵
# ================================================================

SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA = [SIGMA_X, SIGMA_Y, SIGMA_Z]  # σ¹, σ², σ³

I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
I4 = np.eye(4, dtype=complex)
Z4 = np.zeros((4, 4), dtype=complex)


# ================================================================
# Gamma 矩阵类
# ================================================================

class GammaMatrices:
    """Dirac gamma 矩阵代数

    参数:
        representation: 'dirac' (标准) 或 'chiral' (Weyl)

    属性:
        gamma:  [g0, g1, g2, g3] — 4 个 4×4 gamma 矩阵
        g0, g1, g2, g3:  单独访问
        gamma5:  γ⁵ = i γ⁰ γ¹ γ² γ³
        sigma:   σ^{μν} = (i/2)[γ^μ, γ^ν]
    """

    def __init__(self, representation: str = 'dirac'):
        self.representation = representation
        self._build_gamma()
        self._build_gamma5()
        self._build_sigma()

    def _build_gamma(self):
        """构建 gamma 矩阵"""
        if self.representation == 'dirac':
            # Dirac 表示
            self.g0 = np.block([[I2, Z2], [Z2, -I2]])
            self.g1 = np.block([[Z2, SIGMA_X], [-SIGMA_X, Z2]])
            self.g2 = np.block([[Z2, SIGMA_Y], [-SIGMA_Y, Z2]])
            self.g3 = np.block([[Z2, SIGMA_Z], [-SIGMA_Z, Z2]])

        elif self.representation in ('chiral', 'weyl'):
            # Chiral/Weyl 表示
            self.g0 = np.block([[Z2, I2], [I2, Z2]])
            self.g1 = np.block([[Z2, SIGMA_X], [-SIGMA_X, Z2]])
            self.g2 = np.block([[Z2, SIGMA_Y], [-SIGMA_Y, Z2]])
            self.g3 = np.block([[Z2, SIGMA_Z], [-SIGMA_Z, Z2]])

        else:
            raise ValueError(f"未知表示: {self.representation} (可选 'dirac' 或 'chiral')")

        self.gamma = [self.g0, self.g1, self.g2, self.g3]

    def _build_gamma5(self):
        """γ⁵ = i γ⁰ γ¹ γ² γ³"""
        self.gamma5 = 1j * (self.g0 @ self.g1 @ self.g2 @ self.g3)
        # 确保是厄米的
        if not np.allclose(self.gamma5, self.gamma5.conj().T):
            self.gamma5 = (self.gamma5 + self.gamma5.conj().T) / 2

    def _build_sigma(self):
        """σ^{μν} = (i/2)[γ^μ, γ^ν]"""
        self.sigma = {}
        for mu in range(4):
            for nu in range(mu + 1, 4):
                self.sigma[(mu, nu)] = \
                    (1j / 2) * (self.gamma[mu] @ self.gamma[nu] -
                                self.gamma[nu] @ self.gamma[mu])

    def get_sigma_munu(self, mu: int, nu: int) -> np.ndarray:
        """获取 σ^{μν} = (i/2)[γ^μ, γ^ν]

        参数:
            mu, nu: 时空指标 (0..3)

        返回:
            4×4 复数矩阵
        """
        if mu == nu:
            return Z4.copy()
        key = (min(mu, nu), max(mu, nu))
        mat = self.sigma[key]
        if mu > nu:
            return -mat  # 反对称
        return mat

    # ================================================================
    # 性质检查
    # ================================================================

    def check_clifford_algebra(self) -> dict:
        """验证 Clifford 代数 {γ^μ, γ^ν} = 2 g^{μν} I₄

        返回:
            { (μ,ν): max_deviation } 字典
        """
        g = np.diag([1, -1, -1, -1])
        results = {}
        for mu in range(4):
            for nu in range(mu, 4):
                anticom = self.gamma[mu] @ self.gamma[nu] + \
                          self.gamma[nu] @ self.gamma[mu]
                expected = 2 * g[mu, nu] * I4
                dev = float(np.max(np.abs(anticom - expected)))
                results[(mu, nu)] = dev
        return results

    def check_gamma5_properties(self) -> dict:
        """验证 γ⁵ 性质: (γ⁵)² = I, {γ⁵, γ^μ} = 0, (γ⁵)† = γ⁵

        返回:
            {'square_is_I', 'anticommutes', 'hermitian'}
        """
        return {
            'square_is_I': np.allclose(self.gamma5 @ self.gamma5, I4),
            'anticommutes': all(
                np.allclose(self.gamma5 @ g + g @ self.gamma5, Z4)
                for g in self.gamma
            ),
            'hermitian': np.allclose(self.gamma5, self.gamma5.conj().T),
        }

    def summary(self) -> str:
        lines = [f"GammaMatrices(representation='{self.representation}')"]
        lines.append(f"  γ⁵ 迹 = {np.trace(self.gamma5):.2f}  (应为 0)")
        clifford = self.check_clifford_algebra()
        max_dev = max(clifford.values())
        lines.append(f"  Clifford 代数最大偏差: {max_dev:.2e}  (应为 ≈0)")
        return '\n'.join(lines)


# ================================================================
# Dirac 斜线算符
# ================================================================

def dirac_slash(p_4vec: np.ndarray,
                gm: GammaMatrices = None) -> np.ndarray:
    """p̸ = γ^μ p_μ = γ⁰ p₀ - γ¹ p₁ - γ² p₂ - γ³ p₃

    使用约定: p_μ = (E, px, py, pz), γ^μ p_μ

    参数:
        p_4vec: 四动量 (4,) 或 (N, 4)
        gm:     Gamma 矩阵对象 (默认: Dirac 表示)

    返回:
        4×4 矩阵 或 (N, 4, 4) 张量
    """
    if gm is None:
        gm = GammaMatrices()

    p = np.asarray(p_4vec, dtype=complex)

    if p.ndim == 1:
        # p̸ = γ^μ p_μ = γ⁰ p₀ - γ¹ p₁ - γ² p₂ - γ³ p₃
        # p_4vec 约定为 (E, px, py, pz) = p^μ (逆变), p_μ = (E, -px, -py, -pz)
        return p[0] * gm.g0 - p[1] * gm.g1 - p[2] * gm.g2 - p[3] * gm.g3
    else:
        # 批量: (N, 4, 4)
        result = np.zeros((len(p), 4, 4), dtype=complex)
        # p̸ = γ⁰ p₀ - γ¹ p₁ - γ² p₂ - γ³ p₃
        result += p[:, 0, None, None] * gm.gamma[0]
        for mu in range(1, 4):
            result -= p[:, mu, None, None] * gm.gamma[mu]
        return result


# ================================================================
# Dirac 旋量类
# ================================================================

class DiracSpinor:
    """Dirac 旋量 u(p,s) 和 v(p,s)

    参数:
        gm: Gamma 矩阵对象 (默认: Dirac 表示)
    """

    def __init__(self, gm: GammaMatrices = None):
        self.gm = gm if gm is not None else GammaMatrices()

    # ================================================================
    # 正频率旋量 u(p,s)
    # ================================================================

    def u_spinor(self, p_3vec: np.ndarray, mass: float,
                 spin: int = 1) -> np.ndarray:
        """正频率旋量 u(p,s)

        Dirac 表示中:
            u(p,s) = √(E+m) [         χ_s        ]
                              [ σ·p/(E+m) χ_s  ]

        其中 χ_{+1} = [1,0]^T, χ_{-1} = [0,1]^T

        归一化: ū u = 2m

        参数:
            p_3vec: 三维动量 (px, py, pz)
            mass:   费米子质量
            spin:   +1 (自旋向上) 或 -1 (自旋向下)

        返回:
            (4,) 复数旋量
        """
        p = np.asarray(p_3vec, dtype=complex)
        p2 = float(np.dot(p_3vec, p_3vec))
        E = np.sqrt(p2 + mass**2)

        # 自旋二分量
        if spin == 1:
            chi = np.array([1.0, 0.0], dtype=complex)
        elif spin == -1:
            chi = np.array([0.0, 1.0], dtype=complex)
        else:
            raise ValueError(f"spin 必须是 ±1, 得到 {spin}")

        sigma_dot_p = SIGMA_X * p[0] + SIGMA_Y * p[1] + SIGMA_Z * p[2]
        factor = 1.0 / (E + mass) if E + mass > 1e-15 else 0.0

        upper = np.sqrt(E + mass) * chi
        lower = np.sqrt(E + mass) * factor * (sigma_dot_p @ chi)

        return np.concatenate([upper, lower])

    # ================================================================
    # 负频率旋量 v(p,s)
    # ================================================================

    def v_spinor(self, p_3vec: np.ndarray, mass: float,
                 spin: int = 1) -> np.ndarray:
        """负频率旋量 v(p,s) (反粒子)

        Dirac 表示中:
            v(p,s) = √(E+m) [ σ·p/(E+m) ξ_s  ]
                              [      ξ_s       ]

        其中 ξ_{+1} = [0,1]^T, ξ_{-1} = -[1,0]^T

        归一化: v̄ v = -2m

        参数:
            p_3vec: 三维动量
            mass:   费米子质量
            spin:   +1 或 -1

        返回:
            (4,) 复数旋量
        """
        p = np.asarray(p_3vec, dtype=complex)
        p2 = float(np.dot(p_3vec, p_3vec))
        E = np.sqrt(p2 + mass**2)

        if spin == 1:
            xi = np.array([0.0, 1.0], dtype=complex)
        elif spin == -1:
            xi = np.array([-1.0, 0.0], dtype=complex)
        else:
            raise ValueError(f"spin 必须是 ±1, 得到 {spin}")

        sigma_dot_p = SIGMA_X * p[0] + SIGMA_Y * p[1] + SIGMA_Z * p[2]
        factor = 1.0 / (E + mass) if E + mass > 1e-15 else 0.0

        upper = np.sqrt(E + mass) * factor * (sigma_dot_p @ xi)
        lower = np.sqrt(E + mass) * xi

        return np.concatenate([upper, lower])

    # ================================================================
    # Dirac 共轭
    # ================================================================

    def adjoint(self, psi: np.ndarray) -> np.ndarray:
        """Dirac 共轭 ψ̄ = ψ† γ⁰

        参数:
            psi: 4 分量旋量

        返回:
            (4,) 复数行向量
        """
        return (psi.conj().T @ self.gm.g0).flatten()

    def bar(self, psi: np.ndarray) -> np.ndarray:
        """ψ̄ = ψ† γ⁰ (别名)"""
        return self.adjoint(psi)

    # ================================================================
    # 检验
    # ================================================================

    def check_normalization(self, p_3vec: np.ndarray, mass: float) -> dict:
        """验证旋量归一化 ūu = 2m, v̄v = -2m

        返回:
            {'u_bar_u': ..., 'v_bar_v': ..., 'u_pass': bool, 'v_pass': bool}
        """
        u = self.u_spinor(p_3vec, mass, 1)
        v = self.v_spinor(p_3vec, mass, 1)

        u_bar = self.adjoint(u)
        v_bar = self.adjoint(v)

        ubar_u = float(np.dot(u_bar, u))
        vbar_v = float(np.dot(v_bar, v))

        return {
            'ūu': ubar_u,
            'v̄v': vbar_v,
            'ūu ≈ 2m': abs(ubar_u - 2 * mass) < 1e-10,
            'v̄v ≈ -2m': abs(vbar_v + 2 * mass) < 1e-10,
        }


# ================================================================
# 自旋求和
# ================================================================

def spin_sum_u(p_4vec: np.ndarray, mass: float,
               gm: GammaMatrices = None) -> np.ndarray:
    """正频率旋量自旋求和: Σ_s u(p,s) ū(p,s) = p̸ + m

    参数:
        p_4vec: 四动量 (E, px, py, pz)
        mass:   费米子质量
        gm:     Gamma 矩阵

    返回:
        4×4 复数矩阵
    """
    return dirac_slash(p_4vec, gm) + mass * I4


def spin_sum_v(p_4vec: np.ndarray, mass: float,
               gm: GammaMatrices = None) -> np.ndarray:
    """负频率旋量自旋求和: Σ_s v(p,s) v̄(p,s) = p̸ - m

    参数:
        p_4vec: 四动量
        mass:   费米子质量
        gm:     Gamma 矩阵

    返回:
        4×4 复数矩阵
    """
    return dirac_slash(p_4vec, gm) - mass * I4


# ================================================================
# 验证自旋求和 (用旋量显式构造)
# ================================================================

def spin_sum_u_from_spinors(p_3vec: np.ndarray, mass: float,
                             gm: GammaMatrices = None) -> np.ndarray:
    """用旋量显式构造 Σ_s u(p,s) ū(p,s)

    参数:
        p_3vec: 三维动量
        mass:   费米子质量
        gm:     Gamma 矩阵

    返回:
        4×4 矩阵 (应与 p̸+m 一致)
    """
    ds = DiracSpinor(gm)
    total = np.zeros((4, 4), dtype=complex)
    for s in [1, -1]:
        u = ds.u_spinor(p_3vec, mass, s)
        u_bar = ds.adjoint(u)
        total += np.outer(u, u_bar)
    return total


def spin_sum_v_from_spinors(p_3vec: np.ndarray, mass: float,
                             gm: GammaMatrices = None) -> np.ndarray:
    """用反旋量显式构造 Σ_s v(p,s) v̄(p,s)"""
    ds = DiracSpinor(gm)
    total = np.zeros((4, 4), dtype=complex)
    for s in [1, -1]:
        v = ds.v_spinor(p_3vec, mass, s)
        v_bar = ds.adjoint(v)
        total += np.outer(v, v_bar)
    return total


# ================================================================
# 双线性协变量
# ================================================================

def bilinear(psi_bar: np.ndarray, Gamma: np.ndarray,
             psi: np.ndarray) -> complex:
    """双线性协变量 ψ̄ Γ ψ

    Γ 类型:
        标量 (S):  Γ = I₄
        赝标量 (P): Γ = γ⁵
        矢量 (V):  Γ = γ^μ
        轴矢量 (A): Γ = γ^μ γ⁵
        张量 (T):  Γ = σ^{μν}

    参数:
        psi_bar: Dirac 共轭旋量 ψ̄ (4,) 行向量
        Gamma:   4×4 Dirac 矩阵
        psi:     旋量 (4,)

    返回:
        复数标量 ψ̄ Γ ψ
    """
    return complex(psi_bar @ Gamma @ psi)


# 便捷函数: 主要双线性

def bilinear_scalar(psi_bar: np.ndarray, psi: np.ndarray) -> complex:
    """ψ̄ ψ (标量)"""
    return complex(psi_bar @ psi)


def bilinear_pseudoscalar(psi_bar: np.ndarray, psi: np.ndarray,
                           gm: GammaMatrices = None) -> complex:
    """ψ̄ γ⁵ ψ (赝标量)"""
    if gm is None:
        gm = GammaMatrices()
    return complex(psi_bar @ gm.gamma5 @ psi)


def bilinear_vector(psi_bar: np.ndarray, psi: np.ndarray,
                     mu: int, gm: GammaMatrices = None) -> complex:
    """ψ̄ γ^μ ψ (矢量流)"""
    if gm is None:
        gm = GammaMatrices()
    return complex(psi_bar @ gm.gamma[mu] @ psi)


def bilinear_axial(psi_bar: np.ndarray, psi: np.ndarray,
                    mu: int, gm: GammaMatrices = None) -> complex:
    """ψ̄ γ^μ γ⁵ ψ (轴矢量流)"""
    if gm is None:
        gm = GammaMatrices()
    return complex(psi_bar @ (gm.gamma[mu] @ gm.gamma5) @ psi)


# ================================================================
# Dirac 方程检验
# ================================================================

def dirac_equation_check(spinor: np.ndarray, p_4vec: np.ndarray,
                          mass: float, gm: GammaMatrices = None) -> dict:
    """检验 Dirac 方程 (p̸ - m) ψ ≈ 0

    对于 u 旋量: (p̸ - m) u = 0
    对于 v 旋量: (p̸ + m) v = 0

    参数:
        spinor: 4 分量旋量
        p_4vec: 四动量 (E, px, py, pz)
        mass:   费米子质量
        gm:     Gamma 矩阵

    返回:
        {'norm': 残差范数, 'passes': bool, 'is_u': ...}
    """
    p_slash = dirac_slash(p_4vec, gm)
    residual_u = (p_slash - mass * I4) @ spinor
    residual_v = (p_slash + mass * I4) @ spinor

    norm_u = float(np.linalg.norm(residual_u))
    norm_v = float(np.linalg.norm(residual_v))

    tol = 1e-10 * max(float(np.linalg.norm(spinor)), 1.0)
    is_u = norm_u < tol
    is_v = norm_v < tol

    return {
        'norm_u': norm_u,
        'norm_v': norm_v,
        'is_u_solution': is_u,
        'is_v_solution': is_v,
        'passes': is_u or is_v,
    }

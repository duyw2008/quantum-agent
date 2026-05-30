"""多体量子工具 — 张量积 / 偏迹 / 纠缠度量

支持复合量子系统的构造与分析。
"""

import numpy as np


# ═══════════════════════════════════════════════════════════
# 张量积 & 复合系统
# ═══════════════════════════════════════════════════════════

def tensor(*args):
    """多矩阵/多向量的张量积 (Kronecker 积)

    Parameters
    ----------
    *args : np.ndarray
        任意数量的矩阵或向量

    Returns
    -------
    np.ndarray
        张量积结果

    Examples
    --------
    >>> rho_AB = tensor(rho_A, rho_B)     # 复合密度矩阵
    >>> psi = tensor(psi_A, psi_B)         # 复合态向量
    >>> H = tensor(H_A, I_B) + tensor(I_A, H_B)  # 相互作用
    """
    result = args[0]
    for a in args[1:]:
        result = np.kron(result, a)
    return result


def partial_trace(rho, dims, keep=0):
    """偏迹 — 对部分子系统求迹

    Tr_keep[rho] — 保留 keep 子系统, 对其余求迹.

    Parameters
    ----------
    rho : np.ndarray
        复合密度矩阵 (d1*d2*... × d1*d2*...)
    dims : tuple
        各子系统维度 (d1, d2, ...)
    keep : int or tuple
        保留的子系统索引 (从 0 开始)
    """
    rho = np.asarray(rho)
    dims = list(dims)
    n = len(dims)
    if isinstance(keep, int):
        keep = (keep,)

    rho_t = rho.reshape(dims + dims)
    trace_out = sorted([i for i in range(n) if i not in keep], reverse=True)

    for i in trace_out:
        rho_t = np.trace(rho_t, axis1=i, axis2=i + len(dims) - len(trace_out))

    kd = int(np.prod([dims[i] for i in keep]))
    return rho_t.reshape(kd, kd)



def entropy_vn(rho):
    """冯·诺依曼熵: S(ρ) = -Tr[ρ log₂ ρ]

    两体纯态: S(ρ_A) = S(ρ_B) = 纠缠熵。

    Parameters
    ----------
    rho : np.ndarray
        密度矩阵

    Returns
    -------
    float
        熵值 (以 bit 为单位)
    """
    eigenvals = np.linalg.eigvalsh(rho)
    eigenvals = eigenvals[eigenvals > 1e-15]
    if len(eigenvals) == 0:
        return 0.0
    return float(-np.sum(eigenvals * np.log2(eigenvals)))


def concurrence(psi):
    psi = np.asarray(psi)
    """Wootters 并发度 — 两量子比特纠缠度量

    对两体纯态 |ψ⟩ = α|00⟩ + β|01⟩ + γ|10⟩ + δ|11⟩:
        C = 2|αδ - βγ|

    对混合态用 Wootters 公式 (concurrence of formation)。

    Parameters
    ----------
    psi : np.ndarray
        纯态向量 (shape (4,)) 或密度矩阵 (shape (4,4))

    Returns
    -------
    float
        并发度 C ∈ [0, 1]。0 = 可分离, 1 = 最大纠缠 (Bell 态)
    """
    if psi.ndim == 1:
        # 纯态: C = 2|αδ - βγ|
        if len(psi) != 4:
            raise ValueError("Concurrence for 2-qubit systems only (dim=4)")
        return float(2 * abs(psi[0] * psi[3] - psi[1] * psi[2]))

    # 混合态: Wootters 公式
    if psi.shape != (4, 4):
        raise ValueError("Concurrence for 2-qubit systems only")

    sigma_y = np.array([[0, -1j], [1j, 0]])
    R = psi @ tensor(sigma_y, sigma_y) @ psi.conj() @ tensor(sigma_y, sigma_y)
    eigenvals = np.sqrt(np.maximum(np.sort(np.real(np.linalg.eigvals(R)))[::-1], 0))
    return float(max(0, eigenvals[0] - eigenvals[1] - eigenvals[2] - eigenvals[3]))


def bell_states():
    """返回 4 个 Bell 态的字典

    |Φ⁺⟩, |Φ⁻⟩, |Ψ⁺⟩, |Ψ⁻⟩ — 每个 concurrence = 1.0
    """
    phi_plus = np.array([1, 0, 0, 1]) / np.sqrt(2)
    phi_minus = np.array([1, 0, 0, -1]) / np.sqrt(2)
    psi_plus = np.array([0, 1, 1, 0]) / np.sqrt(2)
    psi_minus = np.array([0, 1, -1, 0]) / np.sqrt(2)
    return {
        'Φ⁺': phi_plus,
        'Φ⁻': phi_minus,
        'Ψ⁺': psi_plus,
        'Ψ⁻': psi_minus,
    }

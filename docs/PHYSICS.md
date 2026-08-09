# 物理基础

> 深入理解请读 [KNOWLEDGE_HANDBOOK.md](KNOWLEDGE_HANDBOOK.md)（五卷，25KB）

---

## 量子力学核心

量子态 `|ψ⟩` 在 Fock 空间中演化。测量由 Born 规则 `P = |⟨k|ψ⟩|²` 给出。

| 概念 | 函数 | 公式 |
|------|------|------|
| 相干态 | `coherent(N, α)` | `|α⟩ = e^{-|α|²/2} Σ αⁿ/√(n!) |n⟩` |
| 光子统计 | `g2(state)` | g²(0) = ⟨a†a†aa⟩/⟨a†a⟩² |
| 对易子 | `commutator(A, B)` | [x̂, p̂] = iħ |
| Wigner | `wigner(state)` | W(x,p) < 0 ⇒ 非经典 |

---

## 波函数动力学

TDSE: `iħ ∂ψ/∂t = -ħ²/2m ∂²ψ/∂x² + V(x)ψ`

SSFM 分步傅里叶法求解。`PotentialBuilder` 构造任意 V(x)。

---

## 量子场论

### 自由标量场

1+1D 自由标量场的正则量子化：

    φ̂(x) = Σ_k (â_k e^{ikx} + â†_k e^{-ikx}) / √(2ω_k L)

其中 ω_k = √(k² + m²), [â_k, â†_q] = δ_{kq}。

关键可观测量：
- [φ̂(x), φ̂(y)] ≠ 0 在类空区域 — 微观因果性体现在对易子中
- ⟨0|φ̂²|0⟩ = Σ_k 1/(2ω_k L) — 真空涨落 (紫外截断依赖)
- D_F(x-y) = ⟨0|T{φ̂(x)φ̂(y)}|0⟩ — Feynman 传播子

### φ⁴ 相互作用 (微扰)

相互作用拉氏量: L_int = -(λ/4!) φ⁴

Feynman 规则:
- 传播子: i/(p² - m² + iε)
- 顶点: -iλ
- 对称因子: 1/4! → 顶点排列

Wick 定理将时序积 T{φ(x₁)...φ(x_n)} 展开为所有可能的收缩对的乘积。

### 重整化群

φ⁴ 理论在 1+1D 中是超可重整化的 (只发散到对数阶)。
单圈 β 函数: β(λ) = μ dλ/dμ = 3λ²/(16π²)

物理质量 m_phys² = m² + δm = m² + Π(0)
物理耦合 λ_phys = λ + δλ

### QED 散射

拉氏量: L_QED = ψ̄(iγ^μ D_μ - m)ψ - ¼F_{μν}F^{μν}

Feynman 规则:
- 电子传播子: i(p̸ + m)/(p² - m²)
- 光子传播子: -ig_{μν}/k²
- 顶点: -ieγ^μ

Klein-Nishina 公式描述 Compton 散射 γe⁻ → γe⁻ 的微分截面。

### 有效势与对称破缺

单圈有效势 (Coleman-Weinberg):
    V_eff(φ_c) = V_0(φ_c) + (ħ/64π²) Σ_i M_i⁴(φ_c)[ln(M_i²/μ²) - 3/2]

当 m² < 0 时，有效势在 φ ≠ 0 处取得最小值 → 自发对称破缺。
序参量 = ⟨φ⟩ ≠ 0, Goldstone 定理保证出现无质量玻色子。

| 模块 | 核心物理 |
|------|---------|
| `ScalarField` | 自由场对易子, 因果律, 传播子 |
| `LatticePhi4` | 格点精确对角化, 关联函数, 能隙 |
| `scattering` | Wick 定理, φ⁴ 树图振幅, Dyson 级数 |
| `renormalization` | 单圈自能, 顶点修正, counterterm, β 函数 |
| `gauge` | U(1) 规范场, 光子传播子, Ward 恒等式 |
| `dirac` | γ 矩阵, Dirac 旋量, 自旋求和 |
| `qed` | Klein-Nishina, 对湮灭, Møller 散射 |
| `effective_potential` | Coleman-Weinberg, SSB, Goldstone |
| `PathIntegralMC` | 1D 量子力学路径积分 |
| `LatticePhi4MC` | 2D φ⁴ 场构型 Metropolis 采样 |

---

## 知识手册篇章

[KNOWLEDGE_HANDBOOK.md](KNOWLEDGE_HANDBOOK.md):

- 第〇卷：为什么需要 Fock 空间
- 第一卷：闵氏空间与 Fock 空间的缝合
- 第二卷：因果关系在量子力学中的表现
- 第三卷：为什么需要 Wigner 函数
- 第四卷：退相干的机制

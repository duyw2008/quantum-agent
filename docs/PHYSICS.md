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

| 模块 | 说明 |
|------|------|
| `ScalarField` | 标量场对易子 [φ̂(x), π̂(y)] |
| `LatticePhi4` | 格点 φ⁴ 精确对角化 |
| `PathIntegralMC` | 路径积分 Monte Carlo |

---

## 知识手册篇章

[KNOWLEDGE_HANDBOOK.md](KNOWLEDGE_HANDBOOK.md):

- 第〇卷：为什么需要 Fock 空间
- 第一卷：闵氏空间与 Fock 空间的缝合
- 第二卷：因果关系在量子力学中的表现
- 第三卷：为什么需要 Wigner 函数
- 第四卷：退相干的机制

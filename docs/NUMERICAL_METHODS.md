# 核心数值方法

> Quantum Agent 使用的算法详解

---

## 目录

1. [Split-Step Fourier Method (SSFM)](#1-split-step-fourier-method-ssfm)
2. [精确对角化](#2-精确对角化)
3. [RK4 积分](#3-rk4-积分)
4. [Wigner 函数的 Fock 基计算](#4-wigner-函数的-fock-基计算)
5. [路径积分 Monte Carlo](#5-路径积分-monte-carlo)
6. [Metropolis-Hastings 采样](#6-metropolis-hastings-采样)

---

## 1. Split-Step Fourier Method (SSFM)

### 问题

求解一维含时薛定谔方程：

```
iħ ∂ψ/∂t = Ĥψ = (T̂ + V̂)ψ

T̂ = -ħ²/2m · ∂²/∂x²    (动能)
V̂ = V(x)               (势能)
```

### 核心思想

T̂ 和 V̂ 不对易，所以 `e^{-iĤΔt/ħ} ≠ e^{-iT̂Δt/ħ} · e^{-iV̂Δt/ħ}`。

但当 Δt 足够小时，可以对称拆分（Strang splitting）：

```
e^{-iĤΔt/ħ} ≈ e^{-iV̂Δt/2ħ} · e^{-iT̂Δt/ħ} · e^{-iV̂Δt/2ħ}
                ↑ 半步势能       ↑ 全步动能      ↑ 半步势能
```

**误差** `O(Δt³)`，比简单拆分 `O(Δt²)` 高一个量级。

### 为什么选"拆分"

| 算符 | 在对角表象 | 演化操作 |
|------|:---:|------|
| V̂ | 位置空间 | `ψ(x) *= exp(-iV(x)Δt/ħ)` — 逐点相乘 |
| T̂ | 动量空间 | `ψ̃(k) *= exp(-iħk²Δt/2m)` — 逐点相乘 |

在对角表象里，演化算符只是**逐点乘一个相位因子**——`O(N)` 操作。

### 每次时间步

```
ψ(x)  ──────────→  ψ(x) · e^{-iVΔt/2ħ}          ① 半步势能 (实空间)
  ↓ FFT
ψ̃(k)  ──────────→  ψ̃(k) · e^{-iħk²Δt/2m}         ② 全步动能 (动量空间)
  ↓ IFFT
ψ(x)  ──────────→  ψ(x) · e^{-iVΔt/2ħ}          ③ 半步势能 (实空间)
```

**复杂度**：每步 = 2 次 FFT + 3 次逐点乘 = `O(N log N)`。对 `N=1024` 约 0.5ms。

### 为什么不用差分法

| 方法 | 空间精度 | 稳定性 | 速度 |
|------|:---:|:---:|:---:|
| Crank-Nicolson | O(Δx²) | 无条件稳定 | 需解三对角 → 慢 |
| FFT 谱方法 | 谱精度 (指数收敛) | Δt 小即可 | O(N log N) |
| SSFM | 谱精度 | Δt < Δx² 即可 | 最快 |

SSFM 对光滑波函数给出**谱精度**——误差随 Δx 指数衰减，而非多项式。

### 在我们的代码中

```python
# src/qm/wave.py lines 78-96
pe_half = np.exp(-0.5j * Vx * dt / hbar)    # 半步 V̂ 因子
ke_full = np.exp(-1.0j * hbar * k**2 * dt / (2 * mass))  # 全步 T̂ 因子

for step in range(1, n_steps + 1):
    psi *= pe_half          # ① 半步 V̂ (位置空间)
    psi_k = np.fft.fft(psi) # ② FFT → 动量空间
    psi_k *= ke_full        # ③ 全步 T̂ (动量空间)
    psi = np.fft.ifft(psi_k)# ④ IFFT → 位置空间
    psi *= pe_half          # ⑤ 半步 V̂ (位置空间)
```

### 相关函数

| 函数 | 维度 | 文件 |
|------|:---:|------|
| `evolve_ssfm()` | 1D | `src/qm/wave.py` |
| `evolve_ssfm_2d()` | 2D | `src/qm/wave2d.py` |

2D 版本同理——用 `fft2/ifft2` + `KX²+KY²` 代替 `k²`。

---

## 2. 精确对角化

### 问题

求解 Schrödinger 方程 `iħ d|ψ⟩/dt = H|ψ⟩` 在 Fock 空间中。

### 方法

H 是 `N×N` 的 Hermitian 矩阵。对角化后：

```
H = U Λ U†     (Λ = diag(λ₁, ..., λ_N))

|ψ(t)⟩ = U e^{-iΛt/ħ} U† |ψ(0)⟩
         ↑ 酉变换    ↑ 对角矩阵乘 (精确)
```

**复杂度** `O(N³)` 的对角化只需做一次。之后每步 `O(N²)` 计算任意 t 的 |ψ(t)⟩。

### 守恒大检验

因为用的是精确对角化（非差分近似），能量守恒到机器精度 `~10⁻¹⁵`。

### 在我们的代码中

```python
# src/qm/dynamics.py — sesolve()
eigenvals, eigenvecs = np.linalg.eigh(H)
U = eigenvecs
U_dag = U.conj().T

def psi_at(t):
    phases = np.exp(-1j * eigenvals * t / hbar)
    return U @ (phases * (U_dag @ psi0))
```

**优势**：任意 t 无需逐时间步积分。在 agent 中常用 `sesolve` 验证理论预测。

**限制**：仅适用于**不含时**哈密顿量 `H ≠ H(t)`。

---

## 3. RK4 积分

### 问题

求解 Lindblad 主方程：

```
dρ/dt = -i[H, ρ] + Σ_k γ_k (L_k ρ L_k† - ½{L_k† L_k, ρ})
```

右边是非线性的（矩阵乘），不能直接对角化。

### 方法

经典 4 阶 Runge-Kutta：

```
k₁ = f(ρ)
k₂ = f(ρ + Δt·k₁/2)
k₃ = f(ρ + Δt·k₂/2)
k₄ = f(ρ + Δt·k₃)

ρ(t+Δt) = ρ + Δt/6 · (k₁ + 2k₂ + 2k₃ + k₄)
```

误差 `O(Δt⁴)`。

### 在我们的代码中

```python
# src/qm/dynamics.py — mesolve()
for t in tlist[1:]:
    rho = _rk4(lambda r: lindblad_rhs(H, r, c_ops, hbar), rho, dt)
```

**验证**：对单光子衰减，`⟨n⟩(t) = |α|² e^{-γt}` 精确吻合。

---

## 4. Wigner 函数的 Fock 基计算

### 公式

```
W(α) = (2/π) Tr[ρ D(α) Π D†(α)]

D(α) = exp(α a† - α* a)     (位移算符)
Π = (-1)^N                   (宇称算符)
α = (x + ip)/√2              (复相空间坐标)
```

### 计算

对每个网格点 `(x, p)` → 构造 `α` → 位移算符 → 求迹 → W 值。

**复杂度** `O(N_grid² × N³)`，N_grid=61 时约 1 秒。

### 边界效应

有限 Fock 截断在相空间边界处产生数值伪影——这是 Wigner 函数在网格边缘可能显示负值的原因。物理区域（相空间中心）精度完好。

---

## 5. 路径积分 Monte Carlo

### 欧几里得路径积分

虚时间 `τ = it` 下，传播子变为 Boltzmann 权重：

```
Z = ∫ Dx e^{-S_E[x]/ħ}

S_E[x] = ∫₀^β [½m(ẋ/ħ)² + V(x)] dτ
```

离散化：`β = N_slices × Δτ`，路径变为 `(x₁, x₂, ..., x_N)`

```
S_E ≈ Σ_i [½m(x_{i+1} - x_i)²/(ħ²Δτ) + Δτ V(x_i)]
```

### Metropolis 采样

1. 随机扰动路径上所有点
2. 计算 `ΔS = S_E[新] - S_E[旧]`
3. 若 `ΔS < 0` → 接受（新构型"更可能"）
4. 若 `ΔS > 0` → 以概率 `e^{-ΔS}` 接受

热化后，采样获得 Boltzmann 分布 → 可测量基态能量、|ψ₀|²。

### β → ∞ 极限

虚时间越大，系统越接近基态（高能态的 Boltzmann 权重指数衰减）。`β=10, N_slices=100` 对谐振子是良好的基态逼近（E₀ 误差 ~1%）。

### 在我们的代码中

```python
# src/qft/path_integral.py
pimc = PathIntegralMC(V, mass=1.0, hbar=1.0, N_slices=100, beta=10.0)
pimc.thermalize(5000)
E0, err = pimc.ground_state_energy(10000)
x, psi2 = pimc.wavefunction_density(20000)
```

**验证**：谐振子 `E₀ = 0.5`（理论），PIMC 得 `E₀ ≈ 0.49`。

---

## 6. Metropolis-Hastings 采样

通用的 Markov Chain Monte Carlo 方法，用于从目标分布 `P(x) ∝ e^{-S[x]}` 采样。

### 算法

```
1. 取当前构型 x
2. 提议 x' = x + 随机扰动
3. 计算接受比 r = P(x')/P(x) = e^{-(S[x'] - S[x])}
4. 以概率 min(1, r) 接受 x' → 新构型, 否则保留 x
5. 重复 1-4 直至收敛
```

**接受率调谐**：接受率 ≈ 30-50% 为最优。通过 `delta`（步长）控制。太小 → 慢收敛，太大 → 高拒绝率。

---

## 相关文档

- [PHYSICS.md](PHYSICS.md) — 物理基础
- [KNOWLEDGE_HANDBOOK.md](KNOWLEDGE_HANDBOOK.md) — 五卷深度解析
- [USER_GUIDE.md](USER_GUIDE.md) — 函数参考

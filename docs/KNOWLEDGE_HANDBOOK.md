# Quantum Agent 知识手册

> 量子力学与量子场论核心概念深度解析

---

## 目录

**第〇卷：量子力学基础**

0. [为什么需要 Fock 空间](#0-为什么需要-fock-空间) $\cdot$ [相干态](#07-相干态--最经典的量子态)

**第一卷：量子场论基础**

1. [闵可夫斯基空间与 Fock 空间的缝合](#1-闵可夫斯基空间与-fock-空间的缝合)

**第二卷：因果**

2. [因果关系在量子力学中的表现](#2-因果关系在量子力学中的表现)

**第三卷：量子几何**

3. [Bloch 球的物理意义](#3-bloch-球的物理意义)

**第四卷：量子相空间**

4. [为什么需要 Wigner 函数](#4-为什么需要-wigner-函数)

**第五卷：退相干**

5. [退相干的机制](#5-退相干的机制)

**第六卷：对称性与守恒律**

6. [Noether 定理](#6-noether-定理)

**第七卷：封闭系统的三种等价表述**

7. [Heisenberg $\cdot$ Schrödinger $\cdot$ Feynman — 三条路线](#7-封闭系统的三种等价表述)

**第八卷：开放量子系统**

8. [Lindblad 主方程 — 从幺正到非幺正](#8-lindblad-主方程--开放量子系统)

---

## 0. 为什么需要 Fock 空间

### 0.1 固定粒子数 Hilbert 空间的局限

单粒子量子力学的 Hilbert 空间 H_1 维度固定——一个谐振子、一个自旋、一个氢原子。但自然界充满了粒子数**变化的**过程：

| 现象 | 过程 |
|------|------|
| 光子发射/吸收 |原子 $|e\rangle \rightarrow$|
| 粒子-反粒子产生 |$\gamma \rightarrow e^+ + e^-$|
| Cooper 对形成 |$2e^- \rightarrow 一对$|
| 声子激发 | 晶格振动量子数的增减 |

这些过程的始态和终态粒子数不同——不能写在同一个固定维度的矩阵里。

**核心矛盾**：


固定 N 的 Hilbert 空间 H_N:
  |\psi\rangle = \sum cᵢ_1...ᵢ_n |i_1,...,i_N\rangle
  只能描述恰好 N 个粒子的态

真实情况:
  需要同时描述 N=0,1,2,... 的态
  并且算符可以在不同 N 之间跃迁


### 0.2 Fock 空间的构造

Fock 空间把不同粒子数的 Hilbert 空间全部直和在一起：

```
F = H₀ ⊕ H₁ ⊕ H₂ ⊕ H₃ ⊕ ... ⊕ Hₙ ⊕ ...
  真空   1粒子  2粒子  3粒子      n粒子
```

每个 H_n 是 n 个**不可区分**粒子的 Hilbert 空间：

- **玻色子**：对称子空间 $H_n^+$（波函数交换对称）
- **费米子**：反对称子空间 $H_n^-$（Slater 行列式）

Fock 空间的任意态可以写作：


|Ψ\rangle = c_0|0\rangle + \sum_k c_1(k)|k\rangle + \sumₖ_1,ₖ_2 c_2(k_1,k_2)|k_1,k_2\rangle + ...
      ↑真空      ↑1粒子              ↑2粒子


粒子数不再固定——它是可观测的动力学变量。

### 0.3 产生 / 湮灭算符

在 Fock 空间中，粒子数变化通过升降算符实现：

**玻色子**：

[\hat{a}, \hat{a}^\dagger] = 1

\hat{a}|n\rangle   = \sqrt n \cdot |n-1\rangle      湮灭
\hat{a}^\dagger|n\rangle  = \sqrt{n+1} \cdot |n+1\rangle  产生


**费米子**：

{ĉ, \hat{c}^\dagger} = 1,   ĉ^2 = 0

ĉ|0\rangle = 0,  ĉ|1\rangle = |0\rangle
\hat{c}^\dagger|0\rangle = |1\rangle,  \hat{c}^\dagger|1\rangle = 0  (Pauli 不相容)


产生和湮灭算符把不同 N 的子空间连接起来，是 Fock 空间的"桥梁算符"。

### 0.4 为什么恰好是谐振子代数

玻色子产生/湮灭算符的对易关系 [$\hat{a}, \hat{a}^\dagger] = 1$ 和量子谐振子的升降算符**一模一样**。这不是巧合——本质原因：


量子谐振子能级:       |0\rangle, |1\rangle, |2\rangle, ...
Fock 态的粒子数:      |0\rangle, |1\rangle, |2\rangle, ...
                          ↑
                  数学结构完全相同


但物理意义不同：
- 谐振子中，$|n\rangle$ 是第 n 个**能量本征态**，能量为 $\hbar\omega(n+$½)
- Fock 空间中，$|n\rangle$ 是**有 n 个量子激发**的态，能量为 n$\hbar\omega$（不计零点能零点重定义后）

> 同一个数学结构，承载了两种完全不同的物理

### 0.5 Fock 空间解决了什么

| 问题 | 固定 N 的 QM | Fock 空间 |
|------|:---:|:---:|
| 粒子数变化 | ✗ 不同 N 需要不同空间 | ✓ 所有 N 在一个空间 |
| 产生/湮灭 | ✗ 没有这俩算符 |$✓ \hat{a}^\dagger/\hat{a}$|
| 不可区分性 | 手动对称化 | 自动由交换关系保证 |
| 二次量子化 | 需要手动构造 | 自然框架 |
| 真空涨落 | 不存在 |$\langle0\$|$\hat{\phi}^2\$|$0\rangle \neq 0$|

### 0.6 在 Quantum Agent 中的验证

```python
from src.qm import FockBasis, fock, coherent

# 截断到 N=50 的 Fock 空间
fb = FockBasis(N=50)

# |n⟩ 表示 n 个量子激发的态
psi_3 = fock(50, 3)       # 3 量子 Fock 态

# 升降算符验证
import numpy as np
a_on_3 = fb.a @ psi_3     # â|3⟩ → √3|2⟩
a_dag_on_3 = fb.a_dag @ psi_3   # â†|3⟩ → √4|4⟩

# 相干态 = 所有 Fock 态的叠加 (粒子数不固定!)
psi_coh = coherent(50, alpha=2.0)   # |α⟩ = e^-|α|²/2 ∑ αⁿ/√(n!) |n⟩
mean_photon(psi_coh, fb)             # ⟨n⟩ = |α|² = 4.0
g2(psi_coh, fb)                      # g²(0) = 1.0 (Poisson)
```

### 0.7 相干态 — 最"经典"的量子态

相干态 $|\alpha\rangle$ 是量子光学中最重要的态之一，它是最接近经典电磁波的量子态。

**定义**

$$
|\alpha\rangle = e^{-|\alpha|^2/2} \sum_{n=0}^{\infty} \alpha^n/\sqrt{n!} |n\rangle
$$

$\alpha$ 是复数，编码振幅和相位：$\alpha = |\alpha| e^{i\theta}$。$|\alpha|^2$ 是平均光子数 $\langle n\rangle$。相角 $\theta$ 是光场的经典相位。

**核心性质**

| 性质 | 相干态 | 对比 |
|------|:------:|------|
| 光子数分布 | Poisson，均值/方差 = |$\alpha$|^2 | 热态：Bose-Einstein（超 Poisson） |
| g^2(0) | 1.0 | Fock 态：0 (反聚束); 热态：2 (聚束) |
| Mandel Q | 0 | Q<0 亚 Poisson; Q>0 超 Poisson |
|$\Delta x\cdot\Delta p$|$\hbar/2$| **最小不确定态** — 等于基态 |
| Wigner |$处处 \geq 0$| 纯态中唯一全正的 — 最"经典" |
|$\hat{a}$|$\alpha\rangle = \alpha$|$\alpha\rangle$| 湮灭算符本征态 | 拿走一个光子态不变 |
| 损耗后 | 仍是相干态 | 只是 |$\alpha$|^2 衰减 — 保持相干性 |
| 相位噪声 | 标准量子极限 (SQL) | 压缩态可突破 SQL |

**为什么湮灭掉一个光子态不变？**

这是相干态最具标志性的性质。对于任意光子态 $|\psi\rangle$，测到一个光子后态变为 $\hat{a}|\psi\rangle$（归一化后）。对相干态：


\hat{a}|\alpha\rangle = e^{-|\alpha|^2/2} \sum \alpha^n/\sqrt{n!} \hat{a}|n\rangle
     = e^{-|\alpha|^2/2} \sum \alpha^n/\sqrt{n!} \sqrt n |n-1\rangle
     = e^{-|\alpha|^2/2} \sum_{m=0}^{\infty} \alpha^{m+1}/\sqrt{m!} |m\rangle    (令 m=n-1)
     = \alpha \cdot e^{-|\alpha|^2/2} \sum \alpha^m/\sqrt{m!} |m\rangle
     = \alpha |\alpha\rangle


拿走一个光子 $\rightarrow$ 态除了乘以 $\alpha$ 外完全不变。光子数分布维持 Poisson。这意味着你无法通过"偷看"光子来改变相干态的统计性质——它是没有反作用的测量。

**光子统计详析**

相干态的光子数概率服从 Poisson：

$$
P(n) = |\langle n|\alpha\rangle|^2 = e^{-|\alpha|^2} |\alpha|^{2n} / n!
     = e^{-\langle n\rangle} \langle n\rangle^n / n!
$$

均值和方差相等：$\langle$ n$\rangle =$ Var(n) = $|\alpha|^2$。Mandel Q 参数定义：

$$
Q = Var(n)/\langle n\rangle - 1
$$

- Q = 0: Poisson（相干态）
- Q = -1: 无涨落（Fock 态 $|n\rangle$）
- Q > 0: 超 Poisson（热态/混沌光）

经典激光在远高于阈值时输出 Poisson 光子统计 $\rightarrow$ 相干态是激光的量子描述。

**相干态不正交**

不同 $\alpha$ 的相干态不正交：

$$
\langle\beta|\alpha\rangle = e^{-|\alpha|^2/2 - |\beta|^2/2 + \beta*\alpha}
|\langle\beta|\alpha\rangle|^2 = e^{-|\alpha-\beta|^2}
$$

它们形成一个**过完备基**（overcomplete basis）。任何态可以用相干态展开，但展开系数不唯一。关键关系：


(1/\pi) \int d^2\alpha |\alpha\rangle\langle\alpha| = \hat{I}    (单位算符，d^2\alpha = d(Re \alpha) d(Im \alpha))


这是 Glauber-Sudarshan P 表象的基础。

**Displacement 算符**

相干态可以从真空通过位移产生：


|\alpha\rangle = \hat{D}(\alpha) |0\rangle
\hat{D}(\alpha) = exp(\alpha \hat{a}^\dagger - \alpha* \hat{a})    (幺正算符)


$\hat{D}$ 将真空在相空间平移 $|\alpha|$ 并在相位 $\theta$ 方向旋转。这也是相干态是"位移后的真空"说法的来源。

**Wigner 函数**

相干态的 Wigner 函数是一个高斯：

$$
W_\alpha(x, p) = (1/\pi\hbar) exp[-(x-x_0)^2/2\sigma^2 - (p-p_0)^2/(2\sigma_p^2)]
x_0 = \sqrt{2\hbar} Re(\alpha), p_0 = \sqrt{2\hbar} Im(\alpha)
$$

纯量子态中唯一 Wigner 处处 $\geq 0$ 的态——所有其他纯态（Fock、压缩、猫态）都有负值区域。这赋予了相干态"最经典"的地位：Wigner 的非负性意味着它可以用经典概率分布完美描述（尽管量子干涉仍然存在）。

**在实验室生成**

- 高度稳定的激光远高于阈值 $\rightarrow$ 近似相干态
- 衰减经典电流驱动天线 $\rightarrow$ 相干态光子
- 完美的相干态是理想化概念，真实激光有微量相位扩散

**与压缩态的对比**

| | 相干态 | 压缩真空 | 压缩相干态 |
|---|:---:|:---:|:---:|
|$\Delta x$|$\sigma$|$\sigma e^{-r}$|$\sigma e^{-r}$|
|$\Delta p$|$\sigma$|$\sigma e^{+r}$|$\sigma e^{+r}$|
|$\Delta x\cdot\Delta p$|$\hbar/2$|$\hbar/2$|$\hbar/2$|
|$\langle n\rangle$| |$\alpha$|^2 | sinh^2 r | > |$\alpha$|^2 |
| g^2(0) | 1 | > 1 (偶光子) | 非平凡 |

压缩态仍然满足最小不确定度，但将噪声从一个正交分量"挤压"到另一个——用于增强引力波探测（LIGO）的灵敏度。

**代码实现**

```python
coherent(N=50, alpha=2.0+0.5j)

# N: Fock 截断维度 (≥ 2|α|² 保证收敛)
# alpha: 复数 = |α|e^{iθ}
# 
# 算法: 迭代乘积而非调用 factorial，避免溢出
#   norm = exp(-|α|²/2)
#   ap = α⁰, α¹, α², ...       (累积乘积)
#   fact = 0!, 1!, 2!, ...      (累积阶乘)
#   ket[n] = norm × αⁿ / √(n!)
```

**关键实验基准**

| 实验 | 式子 | 期待值 |
|------|------|:---:|
| 平均光子 | $\langle n\rangle = \langle\psi|\hat{a}^\dagger\hat{a}|\psi\rangle$ | $|\alpha|^2$ |
| 方差 | $\Delta n^2 = \langle n^2\rangle - \langle n\rangle^2$ | $|\alpha|^2$ |
| $g^2(0)$ | $\langle\hat{a}^{\dagger 2}\hat{a}^2\rangle / \langle\hat{a}^\dagger\hat{a}\rangle^2$ | 1.0 |
| 湮灭不变 | $\|\hat{a}|\psi\rangle\| - \||\psi\rangle\||$ | $\approx 0$ |
| Wigner 最小值 | $\min(W)$ | $> 0$（全正）|

---

### 1.1 核心问题

经典场论：场 $\phi(x,t)$ 是一个时空函数，值在每点确定。

量子场论：场变成算符 $\hat{\phi}(x)$，作用在 Fock 空间上。粒子是场的量子激发。

**关键缝合**：闵氏坐标 (x,t) 是场算符的参数，Fock 空间的态描述"这个模式有几个粒子"。

闵可夫斯基空间 M^4  ──\rightarrow  场算符 \hat{\phi}(x)  ──\rightarrow  Fock 空间 F
 (x, t)           (x是坐标参数)       {|n_1, n_2, ...\rangle}


### 1.2 第一步：闵氏空间 \rightarrow Fourier 模式分解

经典实标量场 $\phi(x,t)$ 满足 Klein-Gordon 方程 (□ + m^2)$\phi = 0$。Fourier 按平面波展开：


\phi(x,t) = \int d^3k/(2\pi)^3 \cdot 1/\sqrt{2\omega_k} \cdot [a(k) e^-ikx + a*(k) e^+ikx]

其中: kx = \omega_k t - k\cdot x,   \omega_k = \sqrt{k^2 + m^2}


**关键**：每个波矢 **k** 对应一个独立的振动模式——一个经典谐振子，频率 $\omega_k$。

> 一个场 = 无穷多个独立谐振子的集合

闵氏空间的连续性被离散化为可数的模式标签。动量 `k` 成为模式的"名字"。

### 1.3 第二步：每个模式 \rightarrow 量子升降算符

对每个 k 模式做正则量子化（[$\hat{\phi}, \pi$̂] = i$\hbar$），经典振幅提升为算符：

```
a(k)  →  â(k)     湮灭算符
a*(k) →  â†(k)    产生算符

对易关系: [â(k), â†(k')] = (2π)³ δ³(k - k')
```

每个 k 模式现在是一个量子谐振子，有自己的 Fock 空间：


|n_k\rangle  =  第 k 模式有 n_k 个量子激发

\hat{a}(k)|n_k\rangle   = \sqrt n_k \cdot |n_k - 1\rangle
\hat{a}^\dagger(k)|n_k\rangle  = \sqrt{n_k + 1} \cdot |n_k + 1\rangle


**物理意义**：

- $\hat{a}^\dagger(k)$ 产生一个动量为 `k`、能量为 $\omega_k$ 的粒子
- $\hat{a}(k)$ 湮灭一个这样的粒子
- 粒子是场的量子激发——不是独立的实体

### 1.4 第三步：直积 \rightarrow Fock 空间

所有 k 模式的 Fock 空间做张量积：


F = ⨂_k  F_k

基矢:  |n_1, n_2, n_3, ...\rangle = |n_1\rangle \otimes |n_2\rangle \otimes |n_3\rangle \otimes ...


不同 k 模式上的粒子数构成 Fock 空间的完整态。

| 态 | 物理意义 |
|---|---|
| |$0\rangle =$|$0, 0, 0, ...\rangle$| 真空——所有模式基态 |
|$\hat{a}^\dagger(k_1)$|$0\rangle =$|$1, 0, 0, ...\rangle$| 一个动量为 k_1 的粒子 |
|$\hat{a}^\dagger(k_1)\hat{a}^\dagger(k_2)$|$0\rangle$| 两个粒子，动量 k_1 和 k_2 |
|$(\hat{a}^\dagger(k))^n/\sqrt{n!}$|$0\rangle$| n 个动量相同的玻色子 |

**粒子数算符**：


N̂ = \int d^3k \cdot \hat{a}^\dagger(k)\hat{a}(k)      总粒子数（所有模式求和）
\hat{H} = \int d^3k \cdot \omega_k \cdot \hat{a}^\dagger(k)\hat{a}(k)  总能量


### 1.5 第四步：缝合——场算符 \hat{\phi}(x)

场算符将两者统一：

$$
\hat{\phi}(x) = \int d^3k/(2\pi)^3 \cdot 1/\sqrt{2\omega_k} \cdot [\hat{a}(k) e^-ikx + \hat{a}^\dagger(k) e^+ikx]
$$

- **x = (t, x)** 是闵氏坐标——算符的**参数**
- **$\hat{a}(k), \hat{a}^\dagger(k)**$ 作用在 Fock 空间——算符的**算符部分**

$\hat{\phi}(x)$ 就是缝合线。

**场算符的物理效果**：


\hat{\phi}(x)|0\rangle          粒子在时空点 x 产生（自真空涨落）
\langle0|\hat{\phi}(x)\hat{\phi}(y)|0\rangle    粒子从 x 传播到 y（Feynman 传播子）
\langle k_1k_2|\hat{\phi}(x)|0\rangle    在 x 点可以"找到"这两个粒子


闵氏距离 (x-y)^2 直接出现在传播子分母中——时空结构编码在 Fock 空间的矩阵元里。

### 1.6 物理后果

**粒子 = 场的量子激发**：不是"粒子穿行于时空"，而是"场在 Fock 空间中的激发在闵氏空间中表现为粒子"。

**真空不空**：真空 $|0\rangle$ 是 Fock 空间的一个态——所有模式的基态。由于零点能 $\hbar\omega/2$，真空有非零的场涨落：$\langle0|\hat{\phi}(x)^2|0\rangle \neq 0$。Casimir 效应、Lamb 移位都由此而来。

**因果关系在 Fock 空间中编码**：类空间隔 (x-y)^2 > 0 时，[$\hat{\phi}(x), \hat{\phi}(y)] = 0$——场算符的对易子在类空间隔上为零，保证了微观因果性。这是闵氏几何通过场算符对易子映射到 Fock 空间的结果。

### 1.7 图示总结


闵可夫斯基空间 M^4                       Fock 空间 F
                                              │
经典场 \phi(x,t)  ──Fourier──\rightarrow  {k 模式}  ──量子化──\rightarrow  \hat{a}(k), \hat{a}^\dagger(k)
                                              │
                                              ├──\rightarrow  |n_k\rangle  (单模式 Fock 态)
                                              │
                                              └──\rightarrow  |n_1, n_2, ...\rangle  (全 Fock 空间)
                                                    ↑
场算符 \hat{\phi}(x)  \leftarrow────────── 缝合 ──────────────┘
(参数: x∈M^4)      ã(k)e^-ikx + h.c.


| 概念 | 闵氏侧 | Fock 侧 |
|------|:------:|:------:|
| 基本对象 | 时空坐标 x = (t, x) |$粒子数态 \$|$n_1, n_2, ...\rangle$|
| 动力学量 | k (波矢/动量) | n_k (占据数) |
| 观测 |$\hat{\phi}(x)——在 x 处测量场$|$\langle N̂\rangle——统计粒子数$|
| 传播 |$传播子 \langle\hat{\phi}(x)\hat{\phi}(y)\rangle$|$产生\rightarrow湮灭振幅$|
| 因果性 |$(x-y)^2 > 0 \rightarrow 类空间隔$|$[\hat{\phi}(x), \hat{\phi}(y)] = 0$|
| 真空 | — |$\$|$0\rangle (基点)$|

> **闵氏空间告诉你"在哪里"，Fock 空间告诉你"有多少"。**
> **场算符 $\hat{\phi}(x)$ 同时理解两种语言。**

### 1.8 在 Quantum Agent 中的验证

**标量场对易子**：
```python
from src.qft import ScalarField
sf = ScalarField(mass=1.0)
sf.commutator(0, 0)      # 验证等时对易子 [\hat{\phi}(x), \hat{\pi}(y)] = i δ(x-y)
```

**格点 $\phi^4$ 理论**：
```python
from src.qft import LatticePhi4
lat = LatticePhi4(N_sites=10, mass=1.0)
E0, psi0 = lat.ground_state()   # 精确对角化找基态 (真空)
```

**传播子和 Feynman 振幅**：
```python
from src.qft import feynman_amplitude_phi4_2to2
amp = feynman_amplitude_phi4_2to2(s=10, t=-2, u=-8, coupling=1.0)
# 分母含 (s-m²)——闵氏不变量直接出现在振幅里
```

---

## 2. 因果关系在量子力学中的表现

量子力学没有摧毁因果律，而是分层次地重新定义了它。

### 2.1 第一层：Schrödinger 演化——比经典更因果

波函数 $\psi(x,t)$ 完全由初始条件 $\psi(x,0)$ 通过薛定谔方程决定：


给定 \psi(0) ⟶ Schrödinger 方程 i\hbar \partial\psi/\partial t = \hat{H}\psi ⟶ \psi(t) 唯一确定


这是拉普拉斯式的决定论。没有混沌、没有蝴蝶效应——波函数演化严格酉。在这个层面，量子力学比经典力学更"因果"。

### 2.2 第二层：测量——因果断裂

测量坍缩是非酉、不可逆、概率性的：

```
确定论演化 ──测量──→ 随机坍缩 ──确定论演化──→ ...
```

坍缩瞬间相同原因产生不同结果——这是"因果律"断裂的核心位置：


|\psi\rangle = \alpha|↑\rangle + \beta|↓\rangle

测量前: 因果演化确定
测量:   概率 |\alpha|^2 \rightarrow |↑\rangle, |\beta|^2 \rightarrow |↓\rangle
测量后: 因果演化继续


**退相干给出部分解答**：坍缩不是神秘事件，而是系统与环境纠缠 $\rightarrow$ 相位信息不可逆泄漏 $\rightarrow$ 叠加态表观上退化为经典概率混合。但退相干不能解释"为什么得到这一个结果"——只解释了为什么看不到叠加。

### 2.3 第三层：纠缠——相关 \neq 因果

爱因斯坦最著名的困惑：两个纠缠粒子，测量 A 瞬间"影响"B 的结果，无论距离多远：


|Φ^+\rangle = (|↑↑\rangle + |↓↓\rangle)/\sqrt2

测 A 得 ↑ ⟹ B 必定 ↑  (无论 B 在多远)
测 A 得 ↓ ⟹ B 必定 ↓


**这是相关性，不是因果性**。你不能用纠缠超光速传递信息——A 的测量结果是随机的，无法编码信息。量子力学的深刻结论：

> **相关可以不需要因果**

| 经典直觉 | 量子事实 |
|---|---|
| 完美相关 ⟹ 因果联系 | 完美相关 ⇏ 因果联系 |
| 非定域性 ⟹ 超光速因果 | 非定域性 ⇏ 可用的信号 |
| 因果是客观关系 | 因果可能是涌现性质 |

### 2.4 第四层：QFT——微观因果性作为定理

回到闵氏$\times$ Fock 的缝合——在量子场论层面，因果性是一个严格的定理：

```
类空间隔 (x-y)² > 0 时:  [\hat{\phi}(x), \hat{\phi}(y)] = 0
```

两个类空间隔上的场算符对易——在一点的操作不能在另一点被观测到。信息传递速度 $\leq$ c 被**严格保证**。

```
经典因果: A 在 B 的光锥内 → A 可以影响 B
量子因果: [Ô(A), Ô(B)] = 0  ⇔  A 和 B 类空间隔
```

对易子为零就是量子版本的"没有因果联系"——它是一个数学定理，不是假设。

**为什么对易子为零保证因果性？**
- 测量 A = 从 $|\psi\rangle$ 投影到 Ô(A) 的本征态
- 类空间隔上 [Ô(A), Ô(B)] = 0 ⟹ Ô(A) 和 Ô(B) 有共同本征基
- 测 A 后再测 B $\neq$ 测 B 后再测 A（时间序依赖）
- 但类空间隔上两者对易 $\rightarrow$ 测量序可交换 $\rightarrow$ A 不影响 B

### 2.5 第五层：量子因果推断——因果序本身可以叠加

近年最激进的发展：量子版本的 Pearl 因果图。经典系统中国果关系用有向无环图（DAG）表示：

```
经典: 因果图 G 是固定的 DAG 结构
      A → B → C  或  A ← B → C
      因果序唯一确定
```

量子系统中，因果关系本身可能处于叠加态：


量子: 因果序可以叠加 (indefinite causal order)
      |\psi\rangle = (|A then B\rangle + |B then A\rangle)/\sqrt2


**量子 switch 实验**已实现：操作 A 是否在 B 之前——这个"顺序"本身可以是叠加的。这颠覆了关于因果最基本的假设——因果序不可叠加。

### 2.6 因果的五个层次总结

```
第五层  量子因果推断    因果序可以叠加 ─────────→ 量子 switch
第四层  QFT 对易子     [\hat{\phi}(x), \hat{\phi}(y)] = 0  ────→ 微观因果定理
第三层  纠缠          相关 ≠ 因果 ──────────→ Bell 不等式
第二层  测量坍缩      因果断裂 ────────────→ 退相干
第一层  Schrödinger   严格决定论 ──────────→ 酉演化
```

| 层次 | 因果的表现 | 核心数学 |
|---|---|---|
| Schrödinger 演化 | 严格决定论，比经典更因果 |$i\hbar \partial\psi/\partial t = \hat{H}\psi$|
| 测量坍缩 | 因果断裂——相同因不同果 | Born 规则 P = |$\langle k$|$\psi\rangle$|^2 |
| 纠缠 |$相关 \neq 因果$| Bell 态 + CHSH 不等式 |
| QFT 对易子 | 微观因果性作为定理 |$[\hat{\phi}(x), \hat{\phi}(y)] = 0 on spacelike$|
| 量子因果推断 | 因果序可以叠加 | process matrix / quantum switch |

### 2.7 Bell 不等式与定域实在论

Bell 不等式把哲学问题变成可实验检验的数学不等式。

**两个假设**

| 假设 | 含义 |
|------|------|
| 定域性 | Alice 的测量不瞬间影响 Bob（类空间隔无因果） |
| 实在性 | 粒子在被测之前有确定的自旋值（不靠测量"创造"） |

任何同时满足这两条的理论，其预测必须遵守 Bell 不等式。

**CHSH 形式**

Alice 测 a1 或 a2，Bob 测 b1 或 b2。结果 A,B ∈ {+1,-1}。关联函数 E(a,b) = $\langle$ A$\cdot$ B$\rangle$。

```
S = E(a1,b1) + E(a1,b2) + E(a2,b1) - E(a2,b2)
```

定域实在论预测：$|S| \leq 2$。量子力学预测：对 Bell 态 $|$Φ^+$\rangle$，取最优角度时 $|S| = 2\sqrt2 \approx 2.83$。

**为什么违反**

```
角度: a1=0°, a2=45°, b1=22.5°, b2=67.5°

量子预测:
  E(a1,b1) = -cos(45°) = -0.707
  E(a1,b2) = -cos(135°) = +0.707
  E(a2,b1) = -cos(-45°) = -0.707
  E(a2,b2) = -cos(45°) = -0.707

  S = -0.707 + 0.707 + (-0.707) - (-0.707) = -2.828
  |S| = 2.83 > 2  ← 违反！
```

**违反 Bell 不等式的含义**

| 实验结果 | 结论 |
|:---:|------|
| |S|$\leq 2$| 定域实在论成立，量子力学错了 |
| |S| > 2 | 至少放弃一个：定域性或实在性 |

所有实验得 $|S| \approx 2.83$——定域实在论被排除。

**核心洞察**

纠缠态的相关性不能用"事先商量好"来解释。Alice 和 Bob 的粒子在测量前没有确定的自旋值——它们是一个不可分的整体。Bell 不等式不是量子力学的漏洞，是**经典直觉的边界**。

**结论**：量子力学没有摧毁因果律，而是把因果从"A 导致 B"的朴素直觉升级成了更精密的数学结构——从对易子到退相干到量子因果图。因果不是被否定了，而是被**推广**了。


## 3. Bloch 球的物理意义

### 3.1 几何表示

|psi> = cos(theta/2)|0> + e^{i phi} sin(theta/2)|1>

theta 和 phi 对应球面上的两个角度。
北极=|0>, 南极=|1>, 赤道=等权叠加。

### 3.2 混合态 -> 球内

rho = (I + r.sigma)/2

|r|=1 纯态(球面), |r|=0 完全混合(球心), 0<|r|<1 混合态(球内)

### 3.3 密度矩阵——纯态与混合态

密度矩阵是量子态最通用的描述。纯态只是它的特殊情况。

**纯态 = 信息完整**


|\psi\rangle = \alpha|0\rangle + \beta|1\rangle

\rho = |\psi\rangle\langle\psi| = [[|\alpha|^2,  \alpha\beta* ],
              [\alpha*\beta,   |\beta|^2]]

\rho^2 = \rho (幂等检验)

纯态 \neq 单个粒子。100 个粒子全部在同一个 |\psi\rangle \rightarrow 仍是纯态。
关键不是粒子数量，是你知不知道每个粒子在哪个态。


**混合态 = 信息缺失**


50% |0\rangle + 50% |1\rangle  (不知道每个是哪个)

\rho = 0.5|0\rangle\langle0| + 0.5|1\rangle\langle1| = [[0.5, 0  ],
                                [0,   0.5]]

\rho^2 \neq \rho。对角元仍是概率，非对角元消失——系综不能干涉。


**对角元相同 $\neq$ 态相同**

纯态 |+\rangle:        \rho = [[0.5, 0.5],    非对角元 \neq 0, \rho^2 = \rho
                      [0.5, 0.5]]

混合 50/50:      \rho = [[0.5, 0  ],    非对角元 = 0, \rho^2 \neq \rho
                      [0,   0.5]]


测到 $|0\rangle$ 的概率都是 0.5——对角元看不出来。非对角元是量子相干的指纹。

**两种根本不同的概率**

纯态的概率来自量子测量（单个粒子），混合态的概率来自经典无知（系综）：


纯态 |+\rangle:     一个粒子在 (|0\rangle+|1\rangle)/\sqrt2
              测到 |0\rangle 概率 0.5 \leftarrow 量子不确定性
              \rho = [[0.5, 0.5], [0.5, 0.5]]   非对角元 = 量子相干

混合:         100个粒子, 50在|0\rangle, 50在|1\rangle
              随机抽一个测到 |0\rangle 概率 0.5 \leftarrow 经典无知
              \rho = [[0.5, 0], [0, 0.5]]   非对角元 = 0


对角元相同看不出区别——非对角元才是两者的分界线。

**Bloch 表示**


\rho = (I + r\cdot\sigma)/2

|r| = 1 \rightarrow 纯态 (在球面上)
|r| < 1 \rightarrow 混合态 (在球内)
|r| = 0 \rightarrow 完全混合 (球心, \rho = I/2)


**退相干 = 非对角元死亡**


初态 (纯态):    \rho = [[0.5, 0.5],
                     [0.5, 0.5]]    相干态

  ↓ 退相干

终态 (混合):    \rho \rightarrow [[0.5, 0  ],
                     [0,   0.5]]    相干性消失


### 3.4 量子操作 = 旋转

Pauli-X/Z = 绕轴旋转180度, Hadamard = 绕(x+z)/sqrt(2)旋转180度

### 3.5 退相干 = 球心收缩

Bloch向量向球心收缩 = 非对角元衰减 = 退相干可视化

### 3.5 在 Agent 中

qubit(pi/3, pi/4) | bloch_dm(0.6, 0, 0.8) | bloch_len(rho)

---

## 4. 为什么需要 Wigner 函数

### 4.1 经典相空间的缺失

经典力学：粒子在相空间 (x, p) 中有一个确定的点，概率分布 $\rho(x,p)$ 描述统计系综。你可以同时知道位置和动量。

量子力学：不确定性原理禁止同时精确确定 x 和 p。不存在 P(x,p) ——不存在一个同时给出 x 和 p 的正定概率分布。

**核心问题**：


经典:  \rho(x,p) \geq 0,  \int\int \rho dx dp = 1       (合法概率)
量子:  不存在同时的 x-p 概率分布             (不确定性原理)


那还能不能在相空间中"看"量子态？

### 4.2 Wigner 函数的定义

Wigner (1932) 找到了最接近的答案——准概率分布：

$$
W(x,p) = 1/(\pi\hbar) \int_-\infty^\infty \langle x+y|\rhô|x-y\rangle e^-2ipy/\hbar dy
$$

对密度矩阵做傅里叶变换，映射到相空间 (x, p) 平面。

**性质**：

| 性质 | 公式 | 意义 |
|------|------|------|
| 归一化 | ∬ W dx dp = 1 | 总概率为 1 |
| 实值 | `W ∈ ℝ` | 可解释 |
| 边缘分布 |$\int W dp = \langle x$|$\rhô$|$x\rangle$| 积掉 p 得位置概率 |
| 边缘分布 |$\int W dx = \langle p$|$\rhô$|$p\rangle$| 积掉 x 得动量概率 |
| 期望值 |$\langle\hat{A}\rangle = ∬ W(x,p) A_W dx dp$| Weyl 对应 |
| 范围 |$-2/\pi\hbar \leq W \leq +2/\pi\hbar$| 有界 |

### 4.3 负值——Wigner 函数的核心

W(x,p) 可以取负值。这不是 bug，是 feature——它是区分经典和量子的标志。


相干态:   W(x,p) \geq 0  处处为正   (最"经典"的量子态)
真空:     W(x,p) \geq 0  原点高斯峰
Fock |1\rangle: W(0,0) < 0  原点负值！  (纯量子效应)
猫态:     W 有正负交替的干涉条纹   (非经典性的标志)


> 负值 = 非经典性的充分（非必要）标志

负值意味着不能把 W(x,p) 解释为经典概率——正是这个"不合法"暴露了量子力学的本质。

### 4.4 Wigner 函数回答的问题

| 问题 | Wigner 的回答 |
|------|------|
| 量子态在相空间中"长什么样"？ | W(x,p) 给出一张完整的图 |
| 这个态是经典的还是量子的？ | 有负值 = 非经典 |
| 态经历了什么动力学？ | 时间演化 W(x,p,t) 可视化 |
| 退相干如何发生？ |$干涉条纹逐渐消失 \rightarrow 负值消失$|
| 压缩在哪里？ | W(x,p) 椭圆直接显示压缩方向和幅度 |

### 4.5 与其他相空间表示对比

Wigner 不是唯一的准概率分布。三种最常见：


P 函数 (Glauber-Sudarshan)
  最接近经典——相干态 P = \delta 函数
  但高度奇异，非经典态 P 比 \delta 更奇异

Wigner 函数
  居中——负值可见，始终有限，光滑
  最常用的相空间可视化工具

Q 函数 (Husimi)
  恒正（反卷积了 Wigner）
  丢失了干涉细节——最"模糊"


| 表示 | 正定性 | 奇异性 | 适用场景 |
|------|:---:|:---:|------|
| P | ✗ | 高度奇异 | 理论分析 |
| **Wigner** | ✗ (有负值) | ✓ 有限 | **可视化首选** |
| Q | ✓ (恒正) | ✓ 光滑 | 实验层析 |

### 4.6 在 Quantum Agent 中使用

```python
from src.qm import cat, fock, coherent
from src.viz import wigner, plot_wigner

# 猫态 Wigner 函数——双峰 + 干涉条纹 (负值!)
psi_cat = cat(30, alpha=2.0, phi=0)
x, p, W = wigner(psi_cat, N_grid=61, xlim=(-4, 4), ylim=(-4, 4))

W.min()   # -0.43 ← 负值 = 非经典性
W.max()   # +0.64

plot_wigner(x, p, W, save='cat_wigner.png')

# 对比相干态——全正
psi_coh = coherent(30, alpha=2.0)
x2, p2, W2 = wigner(psi_coh, N_grid=61)
W2.min()   # ~0 ← 无负值，最经典
```

**Agent 交互模式**：

```
⚛ > psi = cat(30, 2.0, 0)
⚛ > x, p, W = wigner(psi)
⚛ > plot_wigner(x, p, W)
⚛ > W.min()           → -0.43
```

---

## 5. 退相干的机制

### 5.1 退相干是什么

退相干 = 量子系统与环境纠缠 $\rightarrow$ 相位关系不可逆泄漏 $\rightarrow$ 叠加态表观坍缩为经典概率混合。

它不是量子力学的附加假设，而是 Schrödinger 方程在开放系统中的**必然结果**。


封闭系统:   |\psi_S\rangle = \alpha|0\rangle + \beta|1\rangle           (纯叠加态)
           \rho_S = |\alpha|^2|0\rangle\langle0| + |\beta|^2|1\rangle\langle1| + \alpha\beta*|0\rangle\langle1| + \alpha*\beta|1\rangle\langle0|
                                              ↑ 量子相干项

开放系统:   系统 + 环境 \rightarrow 纠缠 \rightarrow 对环境求迹
           \rho_S \rightarrow |\alpha|^2|0\rangle\langle0| + |\beta|^2|1\rangle\langle1|    (对角——相干项消失)


### 5.2 通用机制：系统-环境纠缠

设系统初态 $|\psi_S\rangle = \alpha|0\rangle + \beta|1\rangle$，环境初态 $|E_0\rangle$。总初态为直积：

$$
|Ψ(0)\rangle = (\alpha|0\rangle + \beta|1\rangle) \otimes |E_0\rangle
$$

相互作用 H_int 使不同系统态"标记"环境为不同态：


|0\rangle|E_0\rangle \rightarrow |0\rangle|E_0(t)\rangle
|1\rangle|E_0\rangle \rightarrow |1\rangle|E_1(t)\rangle

|Ψ(t)\rangle = \alpha|0\rangle|E_0(t)\rangle + \beta|1\rangle|E_1(t)\rangle      (系统-环境纠缠态)


对环境求迹得约化密度矩阵：


\rho_S(t) = Tr_E[|Ψ(t)\rangle\langleΨ(t)|]
       = |\alpha|^2|0\rangle\langle0| + |\beta|^2|1\rangle\langle1| + \alpha\beta*\langle E_1|E_0\rangle|0\rangle\langle1| + \alpha*\beta\langle E_0|E_1\rangle|1\rangle\langle0|
                                     ↑ 退相干因子


**关键**：相干项被环境重叠积分 $\langle$ E_1$|E_0\rangle$ 压低。环境正交性越大（"which-path"信息越完整），退相干越快。


\langle E_1|E_0\rangle \rightarrow 0  时:   \rho_S \rightarrow |\alpha|^2|0\rangle\langle0| + |\beta|^2|1\rangle\langle1|   (完全退相干)


### 5.3 退相干时间尺度

对于宏观物体的质心位置：


\tau_dec \sim \tau_r \times (\lambda_th / \Delta x)^2

\tau_r     = 环境散射的弛豫时间
\lambda_th    = 环境粒子的热 de Broglie 波长
\Delta x      = |x - x'| = 叠加态两个分支的空间分离


| 系统 |$\Delta x$|$\tau_dec$|
|------|:---:|:---:|
| 电子 (原子尺度) | ~10^{-10} m | ~10^6 年 |
| 分子 (介观) | ~10^{-6} m | ~10^{-3} 秒 |
| 尘埃颗粒 | ~10^{-3} m | ~10^{-20} 秒 |
| 猫 (宏观) | ~0.1 m | ~10^{-40} 秒 |

> **退相干是量子-经典过渡的最快已知过程**——宏观叠加态在飞秒内消失。

不是"大物体不服从量子力学"，而是大物体退相干太快以至于从未被观测到叠加态。

### 5.4 主要退相干通道

#### 振幅阻尼 (Amplitude Damping)

能量从系统流向环境。典型：腔光子通过损耗镜泄漏、原子自发辐射。


Lindblad:  Ľ[\sqrt\gamma \hat{a}]

|n\rangle \rightarrow 混合态 (n>0)
|0\rangle \rightarrow |0\rangle (基态不变)
【相干态特殊】|\alpha\rangle \rightarrow |\alpha e^-\gamma t/2\rangle  (保持纯态，只缩小振幅！)


#### 相位阻尼 (Phase Damping / Pure Dephasing)

能量不损失，但相位随机化。典型：原子在气体中弹性碰撞。


Lindblad:  Ľ[\sqrt\gamma_\phi N̂]

|n\rangle \rightarrow |n\rangle                   (能量不变)
\rho_nₘ \rightarrow \rho_nₘ e^-\gamma_\phi(n-m)^2t  (非对角元衰减)
相干态不再纯——相位信息被洗掉


#### 热噪声 (Thermal Bath)

系统耦合到有限温度热库。同时包含振幅阻尼和相位阻尼。


终态: \rho_ss = 热态    (如果 \hat{H} 无驱动)
      \langle n\rangle_ss = n_th = 1/(e^\hbar\omega/kT - 1)


| 通道 | 能量交换 | 相位破坏 | 相干态表现 |
|------|:---:|:---:|------|
| 振幅阻尼 | ✓ | ✗ (部分) | 保持相干，振幅衰减 |
| 相位阻尼 | ✗ | ✓ | 相位信息被洗掉 |
| 热噪声 | ✓ | ✓ | 变成热混合态 |

### 5.5 退相干 vs 坍缩

退相干解决了"为什么看不到宏观叠加"，但**没有解决测量问题**：

| 退相干解释了 | 退相干没解释 |
|---|---|
| 非对角元为什么衰减 | 为什么会得到这一个结果 |
| 为什么 Wigner 函数负值消失 | 哪一个分支被"选中" |
| 经典概率如何涌现 | 量子到经典的过渡的终极原因 |


退相干后:  \rho = |\alpha|^2|0\rangle\langle0| + |\beta|^2|1\rangle\langle1|
问题:     这是"概率 |\alpha|^2 处于 |0\rangle，|\beta|^2 处于 |1\rangle"
          还是"世界分叉了，你只在一个分支里"？

退相干: "为什么看不到叠加"    \leftarrow 已解决
测量问题: "为什么看到这一个"  \leftarrow 仍需诠释


### 5.6 在 Wigner 函数中的表现

退相干在相空间中的可视化：


初态 (猫态):   W(x,p) 有双峰 + 振荡干涉条纹 (负值)
退相干过程:   干涉条纹逐渐衰减 \rightarrow 负值消失
终态 (混合):   W(x,p) = ½W_|\alpha\rangle + ½W_|-\alpha\rangle  (两个正峰，无条纹)


Wigner 函数的**负值**就是量子相干性的相空间指纹——退相干 = 负值消失。

### 5.7 在 Quantum Agent 中模拟

```python
import numpy as np
from src.qm import FockBasis, cat, coherent_dm, mesolve

# 猫态 + 振幅阻尼
fb = FockBasis(N=50)
H = fb.hamiltonian()                     # 谐振子哈密顿量
rho0 = cat(30, alpha=2.0, phi=0)         # 初态: 偶猫态
rho0 = np.outer(rho0, rho0.conj())       # 纯态密度矩阵

t = np.linspace(0, 10, 50)
gamma = 0.2                              # 退相干速率

result = mesolve(H, rho0, t,
    c_ops=[np.sqrt(gamma) * fb.a],       # 振幅阻尼
    e_ops=[fb.n_op])

n_t = np.real(result['expect'][0])
n_t[0], n_t[-1]   # 初始平均光子数 → 衰减后的光子数
```

**Wigner 演化可视化**：

```python
from src.viz import wigner, plot_wigner

# 退相干前的 Wigner (有负值)
x0, p0, W0 = wigner(rho0, N_grid=61)
# 退相干后的 Wigner (负值消失)
x1, p1, W1 = wigner(result['states'][-1], N_grid=61)

plot_wigner(x0, p0, W0, save='before_decoherence.png')
plot_wigner(x1, p1, W1, save='after_decoherence.png')
```

---

---

## 6. Noether 定理——对称性与守恒律

### 6.1 核心方程

物理学最美的桥梁：每一个连续对称性对应一个守恒量。


连续对称性  ⇔  守恒流 j^\mu  ⇔  守恒荷 Q

\partial_\mu j^\mu = 0          Q = \int j^0 d^3x          dQ/dt = 0


**怎么来的**：作用量 `S = $\int$ L d^4x` 在场变换 `$\phi \rightarrow \phi + \delta\phi$` 下不变：


\delta S = 0  \Rightarrow  \partial_\mu j^\mu = 0

其中 j^\mu = (\partial L/\partial(\partial_\mu \phi)) \cdot \delta\phi - J^\mu
           ↑ 共轭动量              ↑ 变换的边界项


### 6.2 一张表说清楚

| 对称性 | 守恒流 | 守恒荷 Q |
|--------|--------|---------|
|$时间平移 t\rightarrow t+\epsilon$|能量-动量张量 $T_0^\mu$| **能量 E** |
|$空间平移 x\rightarrow x+\epsilon$|$T_i^\mu$| **动量 p** |
| 空间旋转 |$角动量张量 M^\mu\nu$| **角动量 L** |
|$U(1) 规范 \psi\rightarrow e^{i\alpha}\psi$|$j^\mu = \psī\gamma^\mu\psi$| **电荷** |
| SU(3) 色规范 | 8 个胶子流 | **色荷** |

### 6.3 为什么美

这个定理不依赖于具体的 Lagrangian——只要作用量有对称性，守恒律自动产生。

**经典力学**：动量守恒、能量守恒、角动量守恒——全是 Noether 定理的特例。

**量子场论**：规范对称性 $\rightarrow$ 电荷守恒。QED 的电荷守恒 = U(1) 对称性。QCD 的色荷守恒 = SU(3) 对称性。

### 6.4 在 Quantum Agent 中验证

时间平移对称 $\rightarrow$ 能量守恒——SSFM 演化：

```
⚛ > res = evolve_ssfm(psi0, grid, dt=0.01, t_max=5)
⚛ > res['energy'][0], res['energy'][-1]
0.2500  0.2500    # 时间平移对称 → 能量守恒
```

空间平移对称 $\rightarrow$ 动能守恒（自由粒子）：


⚛ > psi0 = gaussian_wavepacket(grid, p0=2.0)
⚛ > res = evolve_ssfm(psi0, grid)  # V=0
# \langle p\rangle 始终不变 — 空间平移对称 \rightarrow 动量守恒


### 6.5 对称性破缺

| 类型 | 机制 | 例子 |
|------|------|------|
| 显式破缺 | Lagrangian 本身不对称 |$加了外势 \rightarrow 动量不守恒$|
| 自发破缺 | Lagrangian 对称，基态不对称 | Higgs 机制、超导 |
| 反常 | 经典对称，量子上破坏 |$轴矢流反常 (\pi^0 \rightarrow \gamma\gamma)$|

**自发破缺 + Noether = Goldstone 定理**：每个自发破缺的连续对称性产生一个无质量玻色子。

### 6.6 第一个例子：概率流密度

Noether 定理最直接的量子实例——波函数 U(1) 对称 $\rightarrow$ 概率守恒。

**连续性方程**

薛定谔方程自带：


\partial\rho/\partial t + \nabla\cdot j = 0

\rho = |\psi|^2                          概率密度
j = (\hbar/2mi)(\psi*\nabla\psi - \psi\nabla\psi*) = \hbar/m Im(\psi*\nabla\psi)   概率流密度


**来自薛定谔方程的推导**

薛定谔方程及其复共轭：

$$
i\hbar \partial\psi/\partial t  = -(\hbar^2/2m)\nabla^2\psi + V\psi       (1)
-i\hbar \partial\psi*/\partial t = -(\hbar^2/2m)\nabla^2\psi* + V\psi*     (2)
$$

(1) $\times \psi*$ − (2) $\times \psi$：


i\hbar \psi*\partial\psi/\partial t + i\hbar \psi\partial\psi*/\partial t = -(\hbar^2/2m)(\psi*\nabla^2\psi - \psi\nabla^2\psi*)

左边:  i\hbar \partial(\psi*\psi)/\partial t = i\hbar \partial\rho/\partial t
右边:  -(\hbar^2/2m) \nabla\cdot(\psi*\nabla\psi - \psi\nabla\psi*)


整理：

$$
\partial\rho/\partial t + \nabla\cdot[ (\hbar/2mi)(\psi*\nabla\psi - \psi\nabla\psi*) ] = 0
         ↑
         j = (\hbar/2mi)(\psi*\nabla\psi - \psi\nabla\psi*) = (\hbar/m) Im(\psi*\nabla\psi)
$$

这就是连续性方程 `$\partial\rho/\partial$ t + $\nabla\cdot$ j = 0`。概率不能凭空产生或消失——只能流动。

**来源**：波函数相位旋转 `$\psi \rightarrow$ e^{i$\alpha}\psi$` 是 U(1) 对称性。拉氏量不变 $\rightarrow$ Noether 给出守恒流 `j^$\mu = (\rho,$ j)`。

**物理直觉**

| 态 |$\rho$| j |
|---|------|------|
| 静止高斯波包 | 包络 |$j \approx 0$|
| 运动波包 e^{ipx} | 包络 |$j = (\hbar k/m) \rho = v\cdot\rho$|
| 驻波 cos(kx) | 振荡 | j = 0 (无净流动) |
| 平面波 e^{ikx} | 常数 |$j = \hbar k/m (均匀流)$|

概率流 = "概率往哪走"的矢量场。和水流、电流完全相同的数学结构。

**和电荷守恒的平行**


概率:  \partial_\mu j^\mu_prob = 0    (U(1) 波函数相位)
电荷:  \partial_\mu j^\mu_EM = 0      (U(1) 规范对称)
能量:  \partial_\mu T^{0\mu} = 0      (时间平移对称)


都是 Noether 定理的同一个公式套不同对称性。

---

## 7. 封闭系统的三种等价表述

Heisenberg 1925, Schrödinger 1926, Feynman 1948 — 量子力学有三种完全等价的数学表述。Dyson 1949 给出了三者等价的严格证明。在 Quantum Agent 中，这是 `sesolve()`、`evolve_ssfm()` 和 `PathIntegralMC` 三条路径，对应三种基底选择。

### 7.1 Heisenberg 路线 — 矩阵力学 \rightarrow sesolve


方程:       i\hbar d|\psi\rangle/dt = H|\psi\rangle             抽象态矢量
解:         |\psi(t)\rangle = exp(-iHt/\hbar) |\psi(0)\rangle
计算:       H = U\Lambda U^\dagger \rightarrow |\psi(t)\rangle = U e^{-i\Lambda t/\hbar} U^\dagger|\psi(0)\rangle


**物理脉络:**

1925 年 Heisenberg/Born/Jordan 提出量子力学的矩阵形式——可观测量是矩阵，能级是本征值。态的演化完全由矩阵代数决定，不涉及空间坐标 x。

**代码实现 (dynamics.py:17-58):**

```python
# 1) 对角化 H (Hermitian → 实本征值)
eigvals, eigvecs = np.linalg.eigh(H)

# 2) 对每个时刻: 基变换 → 相位旋转 → 逆变换
for t in tlist:
    U_diag = np.exp(-1j * eigvals * t / hbar)       # e^{-iλₙt/ℏ}
    psi_t = eigvecs @ (U_diag * (eigvecs.conj().T @ psi0))
    #           ^^^^^^^^   ^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^
    #              U       diag(...)        U^\dagger |ψ₀⟩
```

**对角化推导 — Taylor 展开全过程**

从 $H = U\Lambda U^\dagger$ 到 $|\psi(t)\rangle = U e^{-i\Lambda t/\hbar} U^\dagger|\psi(0)\rangle$，每一步都是恒等式：


① 形式解:  |\psi(t)\rangle = e^{-iHt/\hbar} |\psi(0)\rangle              薛定谔方程的形式积分

② 对角化:  H = U\Lambda U^\dagger                                   Hermitian \rightarrow 存在幺正对角化

③ Taylor:  e^{-iHt/\hbar} = \sum_n (-it/\hbar)^n H^n / n!        矩阵指数的定义

④ 代入 H=U\Lambda U^\dagger, 利用 U^\daggerU=I:
           H^n = (U\Lambda U^\dagger)^n = U \Lambda^n U^\dagger                  中间的 U^\daggerU 全部消掉
           证明: (U\Lambda U^\dagger)^2 = U\Lambda U^\daggerU\Lambda U^\dagger = U\Lambda^2U^\dagger
                (U\Lambda U^\dagger)^3 = U\Lambda^2U^\daggerU\Lambda U^\dagger = U\Lambda^3U^\dagger
                ...归纳

⑤ 代回级数:
           \Sigma(-it/\hbar)^n H^n/n! = \Sigma(-it/\hbar)^n U\Lambda^nU^\dagger/n!
                              = U [\Sigma(-it/\hbar)^n \Lambda^n/n!] U^\dagger
                              = U e^{-i\Lambda t/\hbar} U^\dagger          ↑ 对角阵\rightarrow对每个对角元求指数

⑥ 对角矩阵指数 = 逐元取指数:
           \Lambda=diag(E_1,...,E_N) \rightarrow e^{-i\Lambda t/\hbar}=diag(e^{-iE_1t/\hbar},...,e^{-iE_Nt/\hbar})


核心洞察：对角化把 N 个耦合的 ODE 变成 N 个独立的 ODE。


原始:        i\hbar d|\psi\rangle/dt = H|\psi\rangle           \leftarrow N 分量耦合
对角化后:    i\hbar dc_n/dt = E_n c_n         \leftarrow N 个独立一阶方程
           \rightarrow c_n(t) = c_n(0) e^{-iE_nt/\hbar}


**关键性质:**

| 性质 | 表现 | 原因 |
|------|------|------|
| 能量守恒 |$\sigma_E ~ 10^{15}$| 对角化是恒等式，无近似 |
| 不含时 H | 完美演化 | exp(-iHt) 精确 |
| 含时 H | 不支持 | 矩阵指数不交换 |
| 空间坐标 | 不存在 | 基底是 {|

**算法剖析 — 对角化即精确解:**

矩阵力学的核心洞察：$exp(-iHt)=Uexp(-i\Lambda t) U^\dagger$ 不是近似，是恒等式。对易 H 的矩阵指数等价于在本征基下每个分量独立旋转 $e^{-i\lambda_nt/\hbar}$。

为什么是精确的？因为 U 是幺正矩阵——基变换。$\Lambda$ 是对角矩阵——本征值。在能量本征基 ${|E_n\rangle}$ 下，H 是对角的，演化方程退化为 N 个独立的一阶 ODE：

```
d cₙ/dt = -iEₙ cₙ/ℏ  →  cₙ(t) = cₙ(0) e^{-iEₙt/ℏ}
```

唯一的数值误差来自浮点精度，而非算法近似。代价是 O(N^3) 对角化——只能用于小 Hilbert 空间 (N ~ 100)。

**验证数据:**

H = 2a^\dagger a, |\psi_0\rangle = |\alpha=3\rangle, N=50
E_0 = 18.0000, \sigma_E = 1.3\times 10^{-15}           \leftarrow 能量守恒到机器精度
t = T/4:  |\langle\psi_0|\psi(t)\rangle| = 0.000115 (理论 0.000123)
t = T/2:  |\langle\psi_0|\psi(t)\rangle| = 0.000000 (理论 0.000000)
t = T:    |\langle\psi_0|\psi(t)\rangle| = 0.995228 (理论 1.000000)


谐振子周期 T = 2$\pi/\omega$: 相干态完美回归。T/2 处态正交——相干态 $\pi$ 相移后投影为零。

### 7.2 Schrödinger 路线 — 波动力学 \rightarrow evolve_ssfm


方程:       i\hbar \partial\psi/\partial t = -\hbar^2/2m \partial^2\psi/\partial x^2 + V(x)\psi
                                      \hat{T}               \hat{V}
解:         \psi(x,t) 通过 PDE 数值求解
计算:       SSFM — Strang 拆分 + FFT


**物理脉络:**

1926 年 Schrödinger 提出波动力学——物质是波，$\psi(x)$ 在空间中弥散。演化是 PDE：波包扩散、干涉、隧穿。微分算符 -$\hbar^2/2m \partial^2/\partial$ x^2 代表动能，V(x) 代表势能。

**核心难题: $\hat{T}$ 和 $\hat{V}$ 不对易**


经典:  T(p) + V(x) \rightarrow 直接相加
量子:  \hat{T} 和 \hat{V} 不对易 \rightarrow e^{-i(\hat{T}+\hat{V})\Delta t/\hbar} \neq e^{-i\hat{T}\Delta t/\hbar} \cdot e^{-i\hat{V}\Delta t/\hbar}


这就是 BCH (Baker-Campbell-Hausdorff) 公式的核心限制:

```
e^A e^B = e^{A+B+½[A,B]+¼([A,[A,B]]+[B,[B,A]])+...}
```

$\hat{T}$ 和 $\hat{V}$ 的对易子 [$\hat{T}, \hat{V}] \neq 0$ 意味着泰勒展开的交叉项不消失——直接拆分会引入一阶误差 O($\Delta$ t)。

**解决方案: Strang 拆分**

$$
e^{-i(\hat{T}+\hat{V})\Delta t/\hbar} = e^{-i\hat{V}\Delta t/2\hbar} \cdot e^{-i\hat{T}\Delta t/\hbar} \cdot e^{-i\hat{V}\Delta t/2\hbar} + O(\Delta t^3)
$$

对称拆分消除了一阶对易子项 ½[$\hat{T}, \hat{V}]$——正反两半相消，只留三阶及更高阶的误差。这就是为什么 SSFM 是 O($\Delta$ t^3) 而非 O($\Delta$ t)。

**为什么 FFT？**

在坐标空间 $\hat{V}$ 是对角的（势能 V(x) 在每个点独立作用），但 $\hat{T} = -\hbar^2/2m \partial^2/\partial$ x^2 不是。而在动量空间:


\hat{T} \psĩ(k) = (\hbar^2k^2/2m) \psĩ(k)    \leftarrow 对角! 就是乘法


所以策略是：势能半步(实空间) $\rightarrow$ FFT $\rightarrow$ 动能全步(动量空间) $\rightarrow$ IFFT $\rightarrow$ 势能半步(实空间)。FFT 完成基变换，O(N log N)，比对角化 O(N^3) 快得多——N 可以上万。

**代码实现 (wave.py:83-97):**

```python
pe_half = np.exp(-0.5j * Vx * dt / hbar)       # e^{-iVΔt/2ℏ}
ke_full = np.exp(-1.0j * hbar * k**2 * dt / (2*mass))  # e^{-iTΔt/ℏ}

for step in range(1, n_steps + 1):
    psi = pe_half * psi                         # ① 势能半步 (实空间)
    psi_k = np.fft.fft(psi)                     # ② FFT → 动量空间
    psi_k *= ke_full                            # ③ 动能全步
    psi = np.fft.ifft(psi_k)                    # ④ IFFT → 实空间
    psi = pe_half * psi                         # ⑤ 势能半步
```

**验证数据 — 自由粒子扩散:**

\sigma_param = 1.5, p_0 = 3.0
\Delta x(0)   = 1.0607  (理论 \sigma/\sqrt2 = 1.0607)
\Delta x(3.0) = 1.7584  (理论 1.7584)     \leftarrow 完美吻合!
E = 4.611111111517, \sigma_E = 3.4\times 10^{-15}


自由粒子的扩散公式:


\Delta x(t) = (\sigma/\sqrt2) \sqrt{1 + \hbar^2t^2/m^2\sigma^4}    \leftarrow 波包必然扩散


扩散的根本原因: 波包包含着不同动量的叠加——每个 k 以不同速度移动。初始 $\Delta$ x$\cdot\Delta$ p $\geq \hbar/2$ 意味着动量有分布，而位置扩散是动量分布的必然结果。

**谐振子重生 (revival):**


T = 2\pi/\omega = \pi (\omega=2)
t = T/2:   |\langle\psi_0|\psi\rangle| = 0.999935     \leftarrow 几乎回归
t = T:     |\langle\psi_0|\psi\rangle| = 0.999739     \leftarrow 完整重生
t = 2T:    |\langle\psi_0|\psi\rangle| = 0.998959     \leftarrow 再重生


谐振子的波包每周期 T 回归原位——因为能级等间距，所有分量 $e^{-iE_nt/\hbar}$ 在 t = nT 时同步回归。

### 7.3 两种表示的统一 — 1926 Schrödinger 等价性

两条路解的是**同一个方程** i$\hbar \partial/\partial$ t $|\psi\rangle =$ H$|\psi\rangle$，只是基底选择不同:


矩阵力学:
  基底  = {|0\rangle, |1\rangle, |2\rangle, ...}     Fock 基 (离散)
  态    = \sum_n c_n |n\rangle             列向量
  算符  = 矩阵 H_{mn} = \langle m|H|n\rangle
  演化  = 对角化 \rightarrow 矩阵指数 \rightarrow 投影

波动力学:
  基底  = {|x\rangle}                    位置本征基 (连续)
  态    = \int dx \psi(x) |x\rangle           波函数
  算符  = 微分算符 \hat{H}_x = -\hbar^2/2m \partial^2/\partial x^2 + V(x)
  演化  = PDE \rightarrow SSFM/FFT \rightarrow 波函数序列


**等价性的桥: $\langle$ x$|n\rangle = \phi_n(x)**$

这是谐振子的本征函数——Fock 基向量 $|n\rangle$ 在坐标表象中的波函数。通过它可以在两种表示间任意转换:

```python
# 谐振子: 两种方法应当给出相同的 ⟨x⟩(t), ⟨p⟩(t), Δx(t)
# Fock:  sesolve(Ĥ∝a^\dagger a, |ψ₀⟩=coherent)
# Wave:  evolve_ssfm(V=½mω²x², ψ₀=gaussian)
# → 结果等价
```

**选择指南:**

| 问题 | 方法 | 原因 |
|------|:---:|------|
| 光子统计 g^2(0) 随时间变吗？ | sesolve | Fock 基自然，N~50 |
| 粒子从势垒隧穿的概率？ | evolve_ssfm | 需要坐标空间，N~2048 |
| 密度矩阵退相干动力学？ | mesolve | Lindblad 主方程 |
| 基态能量？ | PIMC | 路径积分 Monte Carlo |
| 自旋在磁场中进动？ | sesolve |$2\times2 矩阵，秒出$|
| 双缝干涉条纹形成？ | evolve_ssfm | 2D SSFM |
| 能量本征谱？ | 对角化 H | 直接本征值问题 |

**两条路线的特征对比:**

矩阵力学 (sesolve)        波动力学 (evolve_ssfm)
基底                |n\rangle Fock 态              |x\rangle 连续坐标
算符                矩阵 H_{mn}              微分算符 \hat{H}_x
数值方法            对角化 O(N^3)              FFT O(N log N)
N 限制              ~100                       ~10^4
能量守恒            精确 (10^{-15})              数值 (~10^{-12})
含时 H              不支持                    不支持 (需扩展)
恰当场景            小N 量子态                散射/隧穿/干涉
| 历史渊源            | Heisenberg 1925          | Schrödinger 1926 |

### 7.4 问题分类 — 什么场景用什么路线

两条路线不只是"换个基底"，它们各自对应根本不同的**物理问题结构**。

**矩阵力学 (sesolve) 适合的问题**

| 问题类型 | 实例 | 为什么 |
|----------|------|--------|
| 光子统计 | g^2(0), Mandel Q, 光子数分布 | Fock 基是光子数的自然语言——a^\dagger a 是对角的，g^2 表达式在 Fock 基下最简 |
| 能谱分析 |$本征值 En, 能隙 \Delta E, 简并度$| 对角化 H 直接输出本征值——这就是能谱 |
| 量子态保真度 |$\langle\psi_1$|$\psi_2\rangle, 态距离, 正交性$| 内积在离散基下是向量点积，天然 |
| 相空间准概率 |$Wigner(\alpha), Husimi Q, Glauber P$|$三种表象都在 Fock 基下通过移位算符 \hat{D}(\alpha) 自然定义$|
| 自旋/量子比特 | Bloch 球, 纠缠, 量子门 |$2\times2 或 2^n\times2^n 矩阵——小维度对角化瞬时完成$|
| 量子光学 | JC 模型, Rabi 振荡, 光子阻塞 |$原子+腔场=复合 Hilbert 空间，对角化求本征态\rightarrow Rabi 频率$|
| 对称性分析 | CASIMIR 算符, 不可约表示 | 群论在离散基底上自然——表示矩阵 |
| 微扰论 |$能级移动 \delta E, 跃迁振幅 \langle m$|V|$n\rangle$| 矩阵元在 Fock 基下是清晰的数，对角化前或后都可算 |

**共性:** 这些问题的物理量**不依赖空间坐标 x**。光子数、自旋投影、能级、纠缠度——都是"内部自由度"，在 Fock/自旋基下即为对角或稀疏矩阵。

---

**波动力学 (evolve_ssfm) 适合的问题**

| 问题类型 | 实例 | 为什么 |
|----------|------|--------|
| 散射 | 势垒隧穿, 透射率/反射率 vs E |$需要 \psi(x) 在坐标空间分成入射/反射/透射三区——坐标自然$|
| 干涉 | 双缝, Aharonov-Bohm, 量子擦除 | 干涉条纹在屏幕上——本质是 |$\psi(x)$|^2 的空间结构 |
| 束缚态 | 有限深势阱的波函数形状 | 能量可能连续也可能离散——对角化只能出离散谱 |
| 波包动力学 | 自由扩散, 谐振子振荡, 势阱捕获 |$\langle x\rangle(t), \Delta x(t), 概率流 j(x,t)——全部在坐标基下定义$|
| 隧穿时间 | Hartman 效应, 粒子在势垒内多久 |$需要追踪波包的 \langle x\rangle(t) 穿过势垒——坐标演化必需$|
| 量子混沌 | Sinai 台球, 能级统计 | 需要边界条件——复杂几何只能上网格 |
| 开放系统 | 波导, 量子点接触, 透射本征道 | 边界条件(入射/出射)在坐标空间最自然 |
| 含时驱动 | Floquet 态, 周期驱动下的隧穿 | 含时 V(x,t) 在 Fock 基下根本不是矩阵——只能 SSFM 逐时步 |
| 非线性 QM | Gross-Pitaevskii 方程 (BEC) | 非线性项 g|$\psi$|^2 在坐标空间是乘法——对角化完全失效 |

**共性:** 这些问题的物理量**依赖空间坐标 x**。波包位置、干涉条纹、隧穿概率——都是"外部自由度"，坐标基下即为直接可观测量。

---

**边界案例 — 两者都可以，选择取决于侧重点**

| 系统 | 矩阵力学角度 | 波动力学角度 |
|------|-------------|-------------|
| 谐振子 |$对角化 H=\omega(a^\dagger a+½): 本征谱精确, 光子统计天然$|$SSFM(V=½m\omega^2x^2): 波包振荡可视化, revival 动态$|
| 势垒 | 转移矩阵法: 透射率 T(E) 解析 | SSFM: T(E) 数值, 含时隧穿过程动画 |
| 量子 Walk | 硬币态+位置 = 直积空间, 对角化演化 | SSFM 等效连续极限, 大尺度扩散 |
| 双阱 | 对称/反对称能级差 = 隧穿频率 |$SSFM: 波包在两个阱之间振荡 \rightarrow 视觉直观$|

**决策树:**

问题是关于"光子/自旋/能级"的?
  Yes $\rightarrow$ 空间坐标无关 $\rightarrow$ sesolve (Fock 对角化)
  No  ↓
问题是关于"波包/条纹/隧穿"的?
  Yes $\rightarrow$ 空间坐标依赖 $\rightarrow$ evolve_ssfm (SSFM)
  No  ↓
问题涉及密度矩阵 + 环境相互作用?
  Yes $\rightarrow$ mesolve (Lindblad RK4)
  No  ↓
只需要基态?
  Yes $\rightarrow$ PIMC (路径积分 MC)

---

**物理直觉的根源 — 对易子决定基底**

为什么谐振子两种方法都适用？因为 [\hat{T}, \hat{V}] \propto [\hat{p}^2, \hat{x}^2] = 2i\hbar(\hat{p}\hat{x} + \hat{x}\hat{p})——不对易，但谐振子的本征态恰好是两种基底的桥梁：Fock 态 |n\rangle 在坐标表象中就是已知的 Hermite 函数 \phi_n(x)。这种"已知的桥"使得两种表示可以互相转换。

而隧穿问题为什么不适合对角化？因为势垒的连续谱——散射态不是束缚态，H 在有限 Fock 截断下没有对应的本征态。坐标空间的 SSFM 自然地处理散射边界条件：波包从左边来，一部分反射、一部分透射——这在 Fock 基下甚至无法表述。

### 7.5 澄清 — 不是"矩阵力学 vs 波动力学"，是"基底选择"

说"矩阵力学不适合隧穿"容易引起误解。准确的说法是：

> **不是矩阵力学不行，是 Fock 基不行。换基底，矩阵力学仍然适用。**

Fock 基 ${|n\rangle}$   ——— 谐振子天然，散射不搭 ──$\rightarrow$ "不能用"
平面波基 ${|k\rangle}$  ——— 散射天然，束缚态不搭 ──$\rightarrow$ S 矩阵理论
坐标基 ${|x\rangle}$    ——— 两者都行但需离散化 ──$\rightarrow$ 有限差分 N$\times$ N 矩阵

```

**特例：一维方势垒的透射率，用转移矩阵法（还是矩阵力学）**

```python
M_left  = [[1,     1   ],    # x < 0: 入射 + 反射
           [ik_1,  -ik_1]]
M_right = [[1,     1   ],    # x > a: 透射
           [ik_1,  -ik_1]]
M_barrier = [[cosh($\kappa$ a), sinh($\kappa$ a)/$\kappa],$  # 0 < x < a
             [$\kappa$ sinh($\kappa$ a), cosh($\kappa$ a)]]

T(E) = $|t|^2 \leftarrow 2\times2$ 矩阵乘完就出结果，全程矩阵力学


这甚至比 SSFM 更方便——不需要跑演化，直接算。关键在于基底换成了"边界处的入射/反射/透射波"，而不是 Fock 态。

**波动力学同样有盲区**

对称地，坐标基 {|x\rangle} 也不是万能的：

| 问题 | 坐标基为何别扭 | Fock/离散基为何自然 |
|------|---------------|-------------------|
| 光子统计 g^2(0) |$需算 \langle\psi$|a^\dagger^2a^2|$\psi\rangle——坐标基要反推回 Fock，绕远路$| Fock 基下是向量点积，一行 |
| 自旋/量子比特 | 自旋没有 x 坐标——"自旋在哪儿？"无意义 |$2\times2 或 2^n\times2^n 矩阵，原生$|
| 纠缠度 concurrence |$坐标空间需两粒子联合概率密度 \int\int$|$\psi(x_1,x_2)$|^2，算部分迹极繁琐 | 直积空间部分迹，einsum 秒出 |
| Wigner 函数 |$W(x,p) = (1/\pi\hbar)\int\langle x+y$|$\rho$|$x-y\rangle e^{-2ipy/\hbar}dy——每个 (x,p) 点都要积分$| Fock 基下位移算符公式，高效得多 |
| 能谱 |$坐标基需解微分方程 -\hbar^2/2m \psi''+V\psi=E\psi 逐本征值 \rightarrow 边值问题$| 对角化 H 一次返回全部本征值 |

**根本原因 — 每一种基底都有它的"舒适区"**


基底选择 = 把哪些算符变成对角矩阵

Fock 基 {$|n\rangle}$:  a^$\dagger$ a 对角  $\rightarrow$  光子数天然  $\rightarrow$  能谱/统计/Wigner 自然
                   $\hat{x}$ 稠密   $\rightarrow$  位置信息散落在所有分量

坐标基 {$|x\rangle}$:  V(x) 对角  $\rightarrow$  势能天然  $\rightarrow$  散射/干涉/隧穿自然
                   a^$\dagger$ a 稠密 $\rightarrow$  光子数信息散落在所有分量

动量基 {$|k\rangle}$:  $\hat{p}$ 对角    $\rightarrow$  动能天然  $\rightarrow$  平面波散射理论自然
                   V(x) 稠密 $\rightarrow$  势能信息散落在所有分量


傅里叶变换就是在 {|x\rangle} 和 {|k\rangle} 之间切换对角的工具。对角化就是在 {|n\rangle} 和 {|E_n\rangle} 之间切换。**所有数值方法归根结底就是"在算符对角的基底下做乘法，在不对角的基底下承受稠密矩阵"——FFT 和对角化只是两种切换方式。**

### 7.6 Feynman 路线 — 路径积分

**核心直觉 — 粒子同时走所有路径**

经典力学说粒子走一条最小作用量轨道。Feynman 说：粒子**同时走所有可能的路径**，每条贡献一个相位 e^{iS[x]/\hbar}。相消干涉抹掉非经典路径，留下经典轨道附近的最强贡献。


传播子 K(b,a) = $\int$ Dx(t) exp(i S[x]/$\hbar)$
                  ^^^^^^
                所有路径的等权叠加


这等价于量子力学的第三种完整表述。Dyson 证明三者给出完全相同的物理预言。

**虚时间旋转 — 量子 \leftrightarrow 经典桥梁**

实时间 e^{iS/\hbar} 是振荡相位，Monte Carlo 无法采样（符号问题）。Wick 旋转 t \rightarrow -i\tau：


实时间: e^{iS/$\hbar}           \leftarrow$ 复数振荡，MC 无法采样
虚时间: e^{-S_E/$\hbar}         \leftarrow$ 实正 Boltzmann 权重，完美适配 MC

Z = $\int$ Dx($\tau)$ exp(-S_E[x]/$\hbar)$


d 维量子系统 \rightarrow (d+1) 维经典统计力学系统，精确映射。\beta \rightarrow \infty 投出基态：


Z = Tr(e^{-$\beta$ H}) = $\sum_n$ e^{-$\beta$ E_n} $\rightarrow$ e^{-$\beta$ E_0}  ($\beta\rightarrow\infty)$
PIMC 采样的路径自然来自基态分布
```

**离散化 (path_integral.py)**

虚时间切 N 片，路径变为经典链 x₁,...,x_N，周期边界：

```
S_E = $\sum_i$ [½m(x_{i+1}-x_i)^2/($\hbar^2\Delta\tau) + \Delta\tau$ V(x_i)]
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^
         动能弹簧 (耦合相邻片)         势能局域
```

单点 Metropolis：每次只扰动一个时间片→接受率 ~30-50%。全路径同时扰动→接受率 ~0.1%。

**验证 (谐振子，E₀_theory=0.500):**

```
$\beta=20,$ N=200:  E_0 = 0.494 $\pm 0.003   (\Delta = 0.006)$
|\psi_0(x)|^2 保真度: 0.946
```

**PIMC 定位:** 无截断误差、天然处理连续谱、可扩展到多体——但比对角化慢，且实时演化有符号问题。

### 7.7 三条路线的统一

Heisenberg / Schrödinger / Feynman 是同一个封闭量子系统的三种等价数学语言：

```
表述         核心方程              数值方法          直觉
──────────────────────────────────────────────────────────
Heisenberg   i$\hbar|\psi$̇$\rangle =$ H$|\psi\rangle$         对角化 (精确)     能谱/跃迁/统计
Schrödinger  i$\hbar\partial\psi/\partial$ t = $\hat{H}\psi(x)$     SSFM (FFT)       波包/干涉/隧穿
Feynman      K = $\int$ Dx e^{iS/$\hbar}$    PIMC (MC)        所有路径求和
```

**三者等价的数学桥:**

```
Heisenberg $\leftarrow\rightarrow$ Schrödinger:  基底变换 $\langle$ x$|n\rangle = \phi_n(x)$
                             同一方程在不同基底下

Heisenberg $\leftarrow\rightarrow$ Feynman:      Trotter 公式
                             e^{-iHt} = lim_{N$\rightarrow\infty} (e^{-iTt/N}$ e^{-iVt/N})^N
                             每步插入完备基 $\int$ dx$|x\rangle\langle$ x|

Schrödinger $\leftarrow\rightarrow$ Feynman:     传播子满足 Schrödinger 方程
                             (i$\hbar\partial/\partial$ t - $\hat{H})$ K(x,t; x_0,0) = i$\hbar\delta(x-x_0)\delta(t)$
```

**选择指南 (封闭系统):**

| 场景 | 路线 | 原因 |
|------|:---:|------|
| 光子统计/能谱/自旋 | Heisenberg | 离散基自然，N~100 |
| 波包/散射/干涉/隧穿 | Schrödinger | 坐标基自然，N~10⁴ |
| 基态能量/波函数/连续谱 | Feynman | 无截断误差，MC 采样 |

**共同本质 — 都是 H 的谱分解:**

```
对角化:   H = $\sum_n$ E_n $|n\rangle\langle$ n|         $\rightarrow$ 本征态 = Fock 态
SSFM:     H = $\int$ dk E_k $|k\rangle\langle$ k|        $\rightarrow$ 本征态 = 平面波 e^{ikx}
PIMC:     Z = Tr(e^{-$\beta$ H}) $\rightarrow$ e^{-$\beta$ E_0} $\rightarrow$ 基态 = Boltzmann 极限
```

三条路线覆盖了封闭量子系统演化的全部疆域——对角化的精确、SSFM 的速度、PIMC 的无截断，互补而非竞争。

---

## 8. Lindblad 主方程 — 开放量子系统

**Lindblad 不是第四种等价表述，而是推广。**

```
封闭系统的三种表述:
  Heisenberg ≡ Schrödinger ≡ Feynman     (Dyson 1949 严格证明)
  全部是幺正演化 $\rightarrow$ Tr($\rho^2)=$常数 $\rightarrow$ dS/dt=0 $\rightarrow$ 可逆

Lindblad 主方程:
  包含封闭部分 -i/$\hbar[H,\rho]$，增加环境耦合项 D[$\rho]$
  非幺正 $\rightarrow$ Tr($\rho^2)$减小 $\rightarrow$ dS/dt $\geq 0 \rightarrow$ 不可逆
  
  当 c_ops=[] 时退化为 Liouville–von Neumann 方程
  ≡ 封闭系统密度矩阵形式 ≡ Heisenberg/Schrödinger/Feynman
```

### 8.1 从封闭到开放

§7 所有内容建立在封闭系统假设上。物理世界中系统永远与环境耦合：

```
耦合到环境的后果:
  • 叠加态 $\rightarrow$ 经典概率混合   (退相干)
  • 激发态 $\rightarrow$ 基态            (能量耗散)
  • 纯态 $\rightarrow$ 混合态            (熵增加)
  • 量子信息 $\rightarrow$ 经典信息      (不可逆)
```

封闭系统的幺正演化无法描述这些——需要推广。

### 8.2 Lindblad 形式 (1976)

Gorini-Kossakowski-Sudarshan-Lindblad 定理：最一般的量子 Markov 主方程是：

```
d$\rho/dt = -i/\hbar$ [H, $\rho] + \sum_k (L_k \rho$ L_k^$\dagger -$ ½{L_k^$\daggerL_k, \rho})$
       ^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       幺正部分       耗散部分 (Lindblad 超算符)
```

逐项解读：

```
L_k $\rho$ L_k^$\dagger$     : "量子跳跃" — 环境测量了系统
                 纯态 $\rightarrow$ 混合态，对角项保留，非对角项衰减

½{L_k^$\daggerL_k, \rho}$  : "反作用" — 确保 Tr($\rho)=1$
                 量子版本的细致平衡，没有它概率不守恒

$\sum_k$            : 多个独立环境通道同时作用
                 光子自发辐射 + 声子散射 + 热噪声 + ...
```

### 8.3 三种基本退相干通道

| 通道 | L | 物理过程 | 效应 |
|------|---|---------|------|
| 振幅阻尼 | $a$ (湮灭算符) | 自发辐射，能量弛豫 | $|1\rangle \rightarrow |0\rangle$，$T_1$ |
| 相位阻尼 | $a^\dagger a$ (粒子数算符) | 弹性散射，无能量损失 | 相干消失，$T_2$ |
| 热耗散 | $\sqrt{\bar{\gamma}} a^\dagger, \sqrt{\gamma} a$ | 与热库热平衡 | $\rightarrow \rho_{th}$ |

振幅阻尼：\|1⟩⟨1\| → (1-γ)\|1⟩⟨1\| + γ\|0⟩⟨0\|——激发态衰减。

相位阻尼：α\|0⟩+β\|1⟩ → \|α\|²\|0⟩⟨0\|+\|β\|²\|1⟩⟨1\|——叠加态变混合态，对角项保留，非对角项消失。

**T₂ ≤ 2T₁ 普适关系:** 振幅阻尼同时致能量弛豫(T₁)和退相干(T₂)，相位阻尼只致退相干。纯退相干率 γ_φ = 1/T₂ − 1/(2T₁) ≥ 0 → T₂ 不可能超过 2T₁。

### 8.4 代码实现 (dynamics.py)

```python
# Lindblad 超算符右侧 — RK4 每步调用
def lindblad_rhs(H, rho, c_ops, hbar):
    drho = -1j/hbar * (H @ rho - rho @ H)       # 幺正 = Liouville–von Neumann
    for L in c_ops:
        LdL = L.conj().T @ L
        drho += L @ rho @ L.conj().T              # 量子跳跃
                - 0.5 * (LdL @ rho + rho @ LdL)   # 反作用 (迹守恒)
    return drho

# RK4 积分
def mesolve(H, rho0, tlist, c_ops, e_ops, hbar):
    rho = rho0.copy()
    for each timestep:
        result['states'].append(rho.copy())
        if not last step:
            rho = rk4(lambda r: lindblad_rhs(H, r, c_ops, hbar), rho, dt)


稳态求解 (steadystate)：d\rho/dt=0 \rightarrow 构造 Liouville 超算符 L (M^2\times M^2)\rightarrow 线性方程 L\cdot vec(\rho)=0 \rightarrow 加迹归一化约束 \rightarrow 直接求解。

### 8.5 从封闭到开放的完整路线


问题 $\rightarrow$ 封闭系统? $\rightarrow$ Yes $\rightarrow$ 能谱/统计?  $\rightarrow$ sesolve (对角化)
                        波包/隧穿?  $\rightarrow$ evolve_ssfm (SSFM)
                        基态能量?   $\rightarrow$ PathIntegralMC (PIMC)
        ↓ No (开放系统)
      退相干/耗散/不可逆 $\rightarrow$ mesolve (Lindblad RK4)
      稳态性质?          $\rightarrow$ steadystate()
```

**为什么需要 Lindblad？**

封闭系统的三条路线都假设 H 是系统的全部——没有外部世界。但：

```
对角化:   $|\psi(t)\rangle =$ e^{-iHt} $|\psi_0\rangle      \rightarrow$ 永远幺正，永远可逆
SSFM:     $\psi(x,t)$ 由 $\hat{H}$ 唯一决定          $\rightarrow$ 波包永远不会自发衰减
PIMC:     Z = Tr(e^{-$\beta$ H})              $\rightarrow$ 只有封闭系统的热平衡

要描述:
  原子自发辐射  $\rightarrow$  需要光子浴 (环境)
  量子比特退相干 $\rightarrow$  需要电磁环境
  激光冷却     $\rightarrow$  需要光和原子耦合
  ...全部超出封闭系统框架
```

Lindblad 补上这个缺口——让量子力学可以描述"系统+环境"中系统的有效演化。它不是第四种等价表述，而是对前三种的推广：从封闭到开放，从幺正到非幺正，从可逆到不可逆。

---

## 附录：参考

- Peskin & Schroeder, *An Introduction to Quantum Field Theory*
- Weinberg, *The Quantum Theory of Fields*, Vol. 1
- Pearl, *Causality* (经典因果推断)
- Chiribella et al., *Quantum computations without definite causal structure* (量子 switch)
- Zurek, W. H., *Decoherence, einselection, and the quantum origins of the classical*, Rev. Mod. Phys. 75, 715 (2003)
- Joos, E. et al., *Decoherence and the Appearance of a Classical World in Quantum Theory*
- Wigner, E. P., *On the Quantum Correction for Thermodynamic Equilibrium*, Phys. Rev. 40, 749 (1932)
- Quantum Agent 源码: `src/qm/`, `src/qft/`

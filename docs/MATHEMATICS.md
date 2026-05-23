# Quantum Agent — 数学模型

> QuTiP 风格量子力学函数库的完整数学基础

---

## 目录

1. [Fock 空间与算符](#1-fock-空间与算符)
2. [量子态](#2-量子态)
3. [算符代数与可观测量](#3-算符代数与可观测量)
4. [光子统计](#4-光子统计)
5. [时间演化](#5-时间演化)
6. [相空间分布](#6-相空间分布)
7. [数值实现](#7-数值实现)

---

## 1. Fock 空间与算符

### 1.1 Fock 基

量子谐振子的本征态构成 Fock 空间的正交归一基：

$$\boxed{\hat{N}|n\rangle = n|n\rangle, \quad n = 0,1,2,\ldots, N-1}$$

在数值计算中截断到有限维 $N$（默认 $N=50$）。

### 1.2 产生与湮灭算符

$$\boxed{\begin{aligned}
\hat{a}|n\rangle &= \sqrt{n}\,|n-1\rangle \\
\hat{a}^\dagger|n\rangle &= \sqrt{n+1}\,|n+1\rangle
\end{aligned}}$$

矩阵元（$N \times N$ 截断）：

$$\boxed{\langle m|\hat{a}|n\rangle = \sqrt{n}\,\delta_{m,n-1}, \quad \langle m|\hat{a}^\dagger|n\rangle = \sqrt{n+1}\,\delta_{m,n+1}}$$

### 1.3 正则对易关系

$$\boxed{[\hat{a}, \hat{a}^\dagger] = \hat{I}, \quad [\hat{a}, \hat{a}] = [\hat{a}^\dagger, \hat{a}^\dagger] = 0}$$

在 $N$ 维截断下，$[\hat{a}, \hat{a}^\dagger] = I - N|N-1\rangle\langle N-1|$，最后一个对角元有截断误差。

### 1.4 坐标与动量算符

定义特征长度 $x_0 = \sqrt{\hbar/m\omega}$ 和特征动量 $p_0 = \sqrt{m\hbar\omega}$：

$$\boxed{\begin{aligned}
\hat{x} &= \frac{x_0}{\sqrt{2}}(\hat{a} + \hat{a}^\dagger) \\
\hat{p} &= i\frac{p_0}{\sqrt{2}}(\hat{a}^\dagger - \hat{a})
\end{aligned}}$$

正则对易关系：

$$\boxed{[\hat{x}, \hat{p}] = i\hbar\hat{I}}$$

在低能子空间（$n \ll N$）精确成立。

### 1.5 数算符

$$\boxed{\hat{N} = \hat{a}^\dagger\hat{a}, \quad \hat{N}|n\rangle = n|n\rangle}$$

### 1.6 宇称算符

$$\boxed{\hat{\Pi} = (-1)^{\hat{N}}, \quad \hat{\Pi}|n\rangle = (-1)^n|n\rangle}$$

在 Wigner 函数计算中起核心作用。

### 1.7 位移算符

$$\boxed{\hat{D}(\alpha) = e^{\alpha\hat{a}^\dagger - \alpha^*\hat{a}}}$$

性质：
- 幺正性：$\hat{D}^\dagger(\alpha)\hat{D}(\alpha) = \hat{I}$
- 位移作用：$\hat{D}^\dagger(\alpha)\hat{a}\hat{D}(\alpha) = \hat{a} + \alpha$
- 作用于真空产生相干态：$\hat{D}(\alpha)|0\rangle = |\alpha\rangle$

数值上通过对角化 $X = \alpha\hat{a}^\dagger - \alpha^*\hat{a}$ 计算矩阵指数：
$$\hat{D}(\alpha) = V\,\text{diag}[e^{\lambda_i}]\,V^{-1}$$

其中 $V$ 是 $X$ 的本征矢矩阵。

---

## 2. 量子态

### 2.1 Fock 态

$$\boxed{|n\rangle = (0,\ldots,0,1,0,\ldots,0)^T}$$

密度矩阵：$\rho_n = |n\rangle\langle n|$

### 2.2 相干态

$$\boxed{|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{N-1}\frac{\alpha^n}{\sqrt{n!}}|n\rangle}$$

**性质**：
- 位移真空态：$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$
- 湮灭算符的本征态：$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$
- 最小不确定态：$\Delta x \cdot \Delta p = \hbar/2$
- 光子数分布：Poisson $P(n) = e^{-|\alpha|^2}|\alpha|^{2n}/n!$
- 平均光子数：$\langle\hat{N}\rangle = |\alpha|^2$

### 2.3 压缩真空态

$$\boxed{|\zeta\rangle = \hat{S}(\zeta)|0\rangle, \quad \hat{S}(\zeta) = e^{(\zeta^*\hat{a}^2 - \zeta\hat{a}^{\dagger 2})/2}}$$

其中 $\zeta = re^{i\theta}$。

**Fock 基展开**（仅偶光子数非零）：

$$\boxed{|\zeta\rangle = \frac{1}{\sqrt{\cosh r}}\sum_{m=0}^{\lfloor N/2\rfloor} \frac{\sqrt{(2m)!}}{2^m m!}(-e^{i\theta}\tanh r)^m|2m\rangle}$$

**性质**：
- 压缩一个正交分量的量子噪声
- 平均光子数：$\langle\hat{N}\rangle = \sinh^2 r$
- $g^2(0) > 1$（光子聚束）

### 2.4 热态

$$\boxed{\rho_{\text{th}} = \sum_{n=0}^{N-1}\frac{\bar{n}^n}{(\bar{n}+1)^{n+1}}|n\rangle\langle n|}$$

其中 $\bar{n} = \langle\hat{N}\rangle$ 是平均热光子数。

**性质**：
- Bose-Einstein 分布：$P(n) = \bar{n}^n / (\bar{n}+1)^{n+1}$
- 最大熵态（对给定的平均光子数）
- $g^2(0) = 2$（热光子聚束）
- Mandel $Q = \bar{n}$

### 2.5 薛定谔猫态

$$\boxed{|\psi_{\text{cat}}\rangle = \frac{|\alpha\rangle + e^{i\phi}|-\alpha\rangle}{\sqrt{2(1 + e^{-2|\alpha|^2}\cos\phi)}}}$$

- $\phi = 0$：偶猫态（仅偶光子数）
- $\phi = \pi$：奇猫态（仅奇光子数）
- 宏观叠加态，$|\alpha| \gg 1$ 时两个分量在相空间几乎正交

### 2.6 纯态诊断

**保真度**：
$$\boxed{F(|\psi_1\rangle, |\psi_2\rangle) = |\langle\psi_1|\psi_2\rangle|^2 \in [0,1]}$$

**纯度**：
$$\boxed{\gamma(\rho) = \text{Tr}[\rho^2] \in [1/N, 1]}$$

$\gamma = 1$ 纯态，$\gamma < 1$ 混合态。

---

## 3. 算符代数与可观测量

### 3.1 对易子与反对易子

$$\boxed{[\hat{A}, \hat{B}] = \hat{A}\hat{B} - \hat{B}\hat{A}, \quad \{\hat{A}, \hat{B}\} = \hat{A}\hat{B} + \hat{B}\hat{A}}$$

### 3.2 期望值

纯态：
$$\boxed{\langle\hat{O}\rangle = \langle\psi|\hat{O}|\psi\rangle = \sum_{m,n} c_m^* c_n O_{mn}}$$

混合态：
$$\boxed{\langle\hat{O}\rangle = \text{Tr}[\rho\hat{O}]}$$

### 3.3 方差

$$\boxed{\text{Var}(\hat{O}) = \langle\hat{O}^2\rangle - \langle\hat{O}\rangle^2}$$

不确定度：$\Delta O = \sqrt{\text{Var}(\hat{O})}$

---

## 4. 光子统计

### 4.1 光子数分布

$$P(n) = |\langle n|\psi\rangle|^2 \quad\text{(纯态)}, \quad P(n) = \rho_{nn} \quad\text{(密度矩阵)}$$

### 4.2 平均光子数

$$\boxed{\langle\hat{N}\rangle = \langle\hat{a}^\dagger\hat{a}\rangle = \sum_n n P(n)}$$

### 4.3 二阶关联函数 $g^{(2)}(0)$

$$\boxed{g^{(2)}(0) = \frac{\langle\hat{a}^\dagger\hat{a}^\dagger\hat{a}\hat{a}\rangle}{\langle\hat{a}^\dagger\hat{a}\rangle^2}}$$

**物理意义**：

| 值 | 含义 | 例子 |
|:---:|------|------|
| $g^{(2)} = 1$ | Poisson 统计，无关联 | 相干态 |
| $g^{(2)} = 2$ | 热光子聚束 | 热态 |
| $g^{(2)} < 1$ | 反聚束（非经典光） | Fock 态 $|n\rangle$，$g^{(2)} = 1 - 1/n$ |
| $g^{(2)} = 0$ | 完美单光子 | $|1\rangle$ |
| $g^{(2)} > 1$ | 光子聚束 | 压缩真空，热态 |

### 4.4 Mandel Q 参数

$$\boxed{Q = \frac{\langle\Delta\hat{N}^2\rangle - \langle\hat{N}\rangle}{\langle\hat{N}\rangle} = \langle\hat{N}\rangle(g^{(2)}(0) - 1)}$$

| 值 | 统计性质 |
|:---:|------|
| $Q = 0$ | Poisson（相干态） |
| $Q > 0$ | 超 Poisson / 经典 |
| $Q < 0$ | 亚 Poisson / 量子 |

---

## 5. 时间演化

### 5.1 Schrödinger 方程

$$\boxed{i\hbar\frac{d}{dt}|\psi(t)\rangle = \hat{H}|\psi(t)\rangle}$$

**形式解**（不含时 $\hat{H}$）：

$$|\psi(t)\rangle = e^{-i\hat{H}t/\hbar}|\psi(0)\rangle$$

数值实现：对角化 $\hat{H} = U\Lambda U^\dagger$，则

$$\boxed{|\psi(t)\rangle = U\,\text{diag}[e^{-iE_k t/\hbar}]\,U^\dagger|\psi(0)\rangle}$$

适用于小 $N$（$\lesssim 500$，对角化 $\mathcal{O}(N^3)$）。

### 5.2 Lindblad 主方程

开放量子系统的时间演化：

$$\boxed{\frac{d\rho}{dt} = -\frac{i}{\hbar}[\hat{H}, \rho] + \sum_k \gamma_k\mathcal{D}[\hat{L}_k]\rho}$$

其中 Lindblad 耗散超算符：

$$\boxed{\mathcal{D}[\hat{L}]\rho = \hat{L}\rho\hat{L}^\dagger - \frac{1}{2}\{\hat{L}^\dagger\hat{L}, \rho\}}$$

- $\hat{L}_k$：坍缩算符（如 $\hat{a}$ 表示光子衰减）
- $\gamma_k$：衰减速率

**物理过程示例**：

| 坍缩算符 | 过程 |
|----------|------|
| $\hat{a}$ | 光子衰减（腔损耗） |
| $\hat{N}$ | 纯退相（dephasing） |
| $\hat{x}$ | 位置测量反作用 |
| $\hat{a}^2$ | 双光子吸收 |

**RK4 数值积分**：

$$\begin{aligned}
k_1 &= f(\rho_t) \\
k_2 &= f(\rho_t + \tfrac{\Delta t}{2}k_1) \\
k_3 &= f(\rho_t + \tfrac{\Delta t}{2}k_2) \\
k_4 &= f(\rho_t + \Delta t\,k_3) \\
\rho_{t+\Delta t} &= \rho_t + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}$$

其中 $f(\rho)$ 是 Lindblad 方程右侧。

### 5.3 稳态求解

$$\frac{d\rho_{ss}}{dt} = 0$$

将密度矩阵拉直为向量：$\text{vec}(\rho) \in \mathbb{C}^{N^2}$

Liouville 超算符 $\mathcal{L}$ 满足：

$$\frac{d}{dt}\text{vec}(\rho) = \mathcal{L}\,\text{vec}(\rho)$$

$\mathcal{L}$ 的显式构造：

$$\boxed{\mathcal{L} = -\frac{i}{\hbar}(\hat{H}\otimes\hat{I} - \hat{I}\otimes\hat{H}^T) + \sum_k\gamma_k\left[\hat{L}_k\otimes\hat{L}_k^* - \frac{1}{2}(\hat{L}_k^\dagger\hat{L}_k\otimes\hat{I} + \hat{I}\otimes(\hat{L}_k^\dagger\hat{L}_k)^T)\right]}$$

稳态条件 $\mathcal{L}\,\text{vec}(\rho_{ss}) = 0$ + 迹约束 $\text{Tr}[\rho_{ss}]=1$，直接求解线性系统。

---

## 6. 相空间分布

### 6.1 Wigner 函数

Wigner 准概率分布是量子态在相空间 $(x,p)$ 中的表示：

$$\boxed{W(x,p) = \frac{2}{\pi}\,\text{Tr}\!\left[\rho\,\hat{D}(\alpha)\,\hat{\Pi}\,\hat{D}(-\alpha)\right]}$$

其中 $\alpha = (x + ip)/\sqrt{2}$，$\hat{\Pi} = (-1)^{\hat{N}}$ 是宇称算符。

**性质**：
- 实函数：$W(x,p) \in \mathbb{R}$
- 归一化：$\iint W(x,p)\,dx\,dp = 1$
- 边缘分布：$\int W(x,p)\,dp = |\psi(x)|^2$，$\int W(x,p)\,dx = |\tilde{\psi}(p)|^2$
- **可取负值**（量子性的标志）
- 高斯态的 Wigner 函数为正

**典型 Wigner 函数**：

| 态 | Wigner 函数特征 |
|----|----------------|
| 真空 $|0\rangle$ | 原点处的高斯峰 |
| 相干态 $|\alpha\rangle$ | 位移的高斯峰（$\alpha$ 处） |
| Fock 态 $|n\rangle$ | 环状结构，$n$ 个负值环 |
| 薛定谔猫态 | 两个高斯峰 + 干涉条纹（负值） |
| 压缩真空 | 压扁的高斯椭圆 |
| 热态 | 展宽的高斯峰 |

### 6.2 Husimi Q 函数

$$\boxed{Q(\alpha) = \frac{1}{\pi}\langle\alpha|\rho|\alpha\rangle}$$

- 恒为非负：$Q(\alpha) \geq 0$
- 与 Wigner 的关系：$Q$ 是 $W$ 的高斯平滑
- 相干态 $|\beta\rangle$ 的 Q 函数：$Q(\alpha) = \frac{1}{\pi}e^{-|\alpha-\beta|^2}$

---

## 7. 数值实现

### 7.1 截断维度选择

Fock 空间截断 $N$ 决定了计算的精确度。

| 物理场景 | 推荐 $N$ |
|----------|:---:|
| 低激发（$\langle\hat{N}\rangle < 2$） | 10-20 |
| 中等激发（$\langle\hat{N}\rangle < 10$） | 30-50 |
| 强场（$\langle\hat{N}\rangle < 50$） | 100-200 |
| 猫态（$|\alpha| < 3$） | 30-50 |

截断误差 $\propto P(N-1)$，即最高 Fock 态的概率。

### 7.2 矩阵指数

$\hat{D}(\alpha) = \exp(X)$ 通过对角化 $X$ 计算：

$$X = V\Lambda V^{-1} \implies e^X = V\,\text{diag}[e^{\lambda_i}]\,V^{-1}$$

对于非正规矩阵 $X$（如位移算符的生成元），使用 `np.linalg.eig` 而非 `eigh`。

### 7.3 稳态求解的稳定性

Liouville 超算符 $\mathcal{L}$ 是奇异矩阵（有一个零本征值对应稳态）。添加迹归一化约束后求解：

$$\begin{pmatrix} \mathcal{L}_{1:N^2-1} \\ \text{vec}(\hat{I})^T \end{pmatrix} \text{vec}(\rho) = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$$

求解后做后处理确保 $\rho$ 的物理性（厄米、正定、迹为 1）。

### 7.4 对易子截断效应

在 $N$ 维截断 Fock 空间中，$[\hat{x}, \hat{p}] = i\hbar\hat{I}$ 仅在 $n \ll N$ 的子空间精确成立。测试中在 $N-5$ 维子空间验证：

$$\|[\hat{x}, \hat{p}]_{\text{sub}} - i\hbar I_{\text{sub}}\|_F < 10^{-14}$$

---

## 附录：符号速查

| 符号 | 含义 | 程序中 |
|------|------|:---:|
| $\hat{a}, \hat{a}^\dagger$ | 湮灭/产生算符 | `fb.a`, `fb.a_dag` |
| $\hat{x}, \hat{p}$ | 坐标/动量 | `fb.x`, `fb.p` |
| $\hat{N}$ | 数算符 | `fb.n_op` |
| $\hat{\Pi}$ | 宇称算符 | `fb.parity` |
| $\hat{D}(\alpha)$ | 位移算符 | `fb.displacement(alpha)` |
| $\hat{H}_{\text{HO}}$ | 谐振子哈密顿量 | `fb.hamiltonian()` |
| $N$ | Fock 空间截断 | `FockBasis(N)` |
| $\alpha$ | 相干态振幅 | `coherent(N, alpha)` |
| $\zeta = re^{i\theta}$ | 压缩参数 | `squeezed(N, zeta)` |
| $\bar{n}$ | 平均热光子数 | `thermal_dm(N, n_th)` |
| $g^{(2)}(0)$ | 二阶关联 | `g2(state, fb)` |
| $Q$ | Mandel 参数 | `mandel_q(state, fb)` |
| $W(x,p)$ | Wigner 函数 | `wigner(state)` |
| $Q(\alpha)$ | Husimi Q 函数 | `qfunc(state)` |

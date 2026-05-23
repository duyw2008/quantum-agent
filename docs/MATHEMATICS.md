# Quantum Agent — 数学模型

本文档完整涵盖 quantum_agent 中使用的所有数学模型、数值方法和物理理论。

---

## 目录

1. [含时薛定谔方程 (TDSE)](#1-含时薛定谔方程-tdse)
2. [Split-Step Fourier Method (SSFM)](#2-split-step-fourier-method-ssfm)
3. [Crank-Nicolson Method (CN)](#3-crank-nicolson-method-cn)
4. [空间离散化与 FFT](#4-空间离散化与-fft)
5. [波函数与可观测量](#5-波函数与可观测量)
6. [矩阵力学](#6-矩阵力学)
7. [势函数](#7-势函数)
8. [本征值问题](#8-本征值问题)
9. [数值稳定性与误差分析](#9-数值稳定性与误差分析)

---

## 1. 含时薛定谔方程 (TDSE)

### 1.1 基本形式

一维含时薛定谔方程（原子单位 ℏ = mₑ = e = 1）：

$$\boxed{i\frac{\partial}{\partial t}\psi(x,t) = \hat{H}\psi(x,t) = \left[-\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi(x,t)}$$

在程序中还原了 ℏ 和 m 参数，通用形式为：

$$i\hbar\frac{\partial\psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2\psi}{\partial x^2} + V(x)\psi$$

### 1.2 形式解

对于不含时哈密顿量，形式解为：

$$\psi(t) = e^{-i\hat{H}t/\hbar}\,\psi(0) = \hat{U}(t)\,\psi(0)$$

其中 $\hat{U}(t) = e^{-i\hat{H}t/\hbar}$ 是幺正时间演化算符（$\hat{U}^\dagger\hat{U} = \hat{I}$）。

### 1.3 守恒量

- **概率守恒**: $\frac{d}{dt}\int|\psi|^2 dx = 0$，即 $\|\psi(t)\| = \|\psi(0)\|$
- **能量守恒**: $\langle\hat{H}\rangle = \text{const}$（对于不含时 H）

---

## 2. Split-Step Fourier Method (SSFM)

### 2.1 Trotter-Suzuki 分解

核心思想：将演化算符分解为动能和势能部分的乘积。

$$\hat{U}(\Delta t) = e^{-i(\hat{T}+\hat{V})\Delta t/\hbar}$$

一阶 Trotter 分解（误差 $\mathcal{O}(\Delta t^2)$）：

$$e^{-i(\hat{T}+\hat{V})\Delta t/\hbar} = e^{-i\hat{V}\Delta t/\hbar}\,e^{-i\hat{T}\Delta t/\hbar} + \mathcal{O}(\Delta t^2)$$

二阶对称 Trotter 分解（误差 $\mathcal{O}(\Delta t^3)$，程序使用此版本）：

$$\boxed{e^{-i(\hat{T}+\hat{V})\Delta t/\hbar} = e^{-i\hat{V}\Delta t/2\hbar}\,e^{-i\hat{T}\Delta t/\hbar}\,e^{-i\hat{V}\Delta t/2\hbar} + \mathcal{O}(\Delta t^3)}$$

### 2.2 算法步骤

**Step 1** — 半步势能演化（坐标空间）：
$$\psi_1(x) = \exp\left[-\frac{iV(x)\Delta t}{2\hbar}\right] \psi(x, t)$$

这里的乘法是逐点（element-wise）的，因为 $\hat{V}$ 在坐标表象是对角的。

**Step 2** — 傅里叶变换到动量空间：
$$\tilde{\psi}_1(k) = \mathcal{F}[\psi_1(x)] = \int_{-\infty}^{\infty} \psi_1(x)\,e^{-ikx}\,dx$$

数值实现使用 FFT：
$$\tilde{\psi}_1[k] = \sum_{j=0}^{N-1} \psi_1[j]\,e^{-2\pi i j k / N}$$

**Step 3** — 动能演化（动量空间）：
$$\tilde{\psi}_2(k) = \exp\left[-\frac{i\hbar k^2\Delta t}{2m}\right] \tilde{\psi}_1(k)$$

动量空间动能算符是对角的 $\hat{T}(k) = \hbar^2 k^2 / 2m$，这也是逐点乘法。

**Step 4** — 逆傅里叶变换回坐标空间：
$$\psi_2(x) = \mathcal{F}^{-1}[\tilde{\psi}_2(k)]$$

**Step 5** — 半步势能演化：
$$\psi(x, t + \Delta t) = \exp\left[-\frac{iV(x)\Delta t}{2\hbar}\right] \psi_2(x)$$

### 2.3 复杂度

每一步的复杂度由 FFT 主导：$\mathcal{O}(N\log N)$，其中 $N$ 是网格点数。

### 2.4 适用范围

- ✅ 光滑势函数（无奇点）
- ✅ 长时间演化（误差积累可控）
- ✅ 谱精度（指数收敛于空间导数）
- ❌ 尖锐势边界（需要额外处理）
- ❌ 势能含时（需要更高阶分解）

---

## 3. Crank-Nicolson Method (CN)

### 3.1 离散化

对 TDSE 使用时间中心差分（Crank-Nicolson 格式）：

$$i\hbar\frac{\psi^{n+1} - \psi^n}{\Delta t} = \hat{H}\frac{\psi^{n+1} + \psi^n}{2}$$

其中 $\psi^n = \psi(x, t_n)$，$t_n = n\Delta t$。

重写为矩阵形式：

$$\boxed{\left(\hat{I} + \frac{i\Delta t}{2\hbar}\hat{H}\right)\psi^{n+1} = \left(\hat{I} - \frac{i\Delta t}{2\hbar}\hat{H}\right)\psi^n}$$

### 3.2 空间离散化

使用三点中心差分近似动能算符：

$$\frac{\partial^2\psi}{\partial x^2}\bigg|_{x_j} \approx \frac{\psi_{j-1} - 2\psi_j + \psi_{j+1}}{\Delta x^2}$$

因此哈密顿量矩阵是三对角的：

$$\hat{H}\psi\big|_j = -\alpha\psi_{j-1} + (2\alpha + V_j)\psi_j - \alpha\psi_{j+1}$$

其中 $\alpha = \frac{\hbar^2}{2m\Delta x^2}$，$V_j = V(x_j)$。

### 3.3 矩阵结构

左侧矩阵 $\hat{A} = \hat{I} + \frac{i\Delta t}{2\hbar}\hat{H}$ 的元素：

$$\boxed{\begin{aligned}
A_{jj} &= 1 + \frac{i\Delta t}{2\hbar}(2\alpha + V_j) \\
A_{j,j+1} = A_{j+1,j} &= -\frac{i\Delta t}{2\hbar}\alpha
\end{aligned}}$$

### 3.4 Thomas 算法求解

三对角系统 $\mathbf{A}\mathbf{x} = \mathbf{b}$ 通过 Thomas 算法以 $\mathcal{O}(N)$ 求解：

**前向消元**：
$$\begin{aligned}
c'_0 &= c_0 / d_0 \\
x'_0 &= b_0 / d_0 \\
\text{for } i &= 1,\ldots,N-1: \\
& d'_i = d_i - a_{i-1}c'_{i-1} \\
& c'_i = c_i / d'_i \quad (i < N-1) \\
& x'_i = (b_i - a_{i-1}x'_{i-1}) / d'_i
\end{aligned}$$

**回代**：
$$\begin{aligned}
x_{N-1} &= x'_{N-1} \\
\text{for } i &= N-2,\ldots,0: \\
& x_i = x'_i - c'_i x_{i+1}
\end{aligned}$$

程序中实现为 `_complex_tridiag_solve()`，直接处理复数运算。

### 3.5 性质

- **无条件稳定**: 任意 $\Delta t$ 都不导致指数增长
- **二阶精度**: 截断误差 $\mathcal{O}(\Delta t^2 + \Delta x^2)$
- **幺正性**: 近似保持概率守恒
- **边界条件**: 默认 Dirichlet ($\psi = 0$ at boundary)

---

## 4. 空间离散化与 FFT

### 4.1 网格定义

$$\boxed{\begin{aligned}
x_j &= x_{\min} + j\Delta x, \quad j = 0,1,\ldots,N-1 \\
\Delta x &= \frac{x_{\max} - x_{\min}}{N-1}
\end{aligned}}$$

### 4.2 动量空间网格

FFT 自然定义动量空间网格：

$$\boxed{k_n = \frac{2\pi n}{N\Delta x}, \quad n = -\frac{N}{2},\ldots,0,\ldots,\frac{N}{2}-1}$$

用 `np.fft.fftfreq` 实现，频率间隔 $\Delta k = 2\pi/(N\Delta x)$。

### 4.3 Nyquist 条件

最大可分辨动量：$k_{\max} = \pi / \Delta x$

对应的最大动能：$E_{\max} = \hbar^2 k_{\max}^2 / 2m = \hbar^2\pi^2/(2m\Delta x^2)$

对于精确模拟，需确保 $\Delta x$ 足够小以满足 $E_{\max} \gg \max(V) + \text{typical kinetic energy}$。

### 4.4 FFT 约定

程序使用标准 `numpy.fft` 约定：

$$\tilde{\psi}_k = \sum_{j=0}^{N-1} \psi_j\,e^{-2\pi i j k / N}$$

$$\psi_j = \frac{1}{N} \sum_{k=0}^{N-1} \tilde{\psi}_k\,e^{2\pi i j k / N}$$

归一化：`np.fft.fft` 无归一化因子，`np.fft.ifft` 除以 $N$。

---

## 5. 波函数与可观测量

### 5.1 波函数表示

$$\psi(x) = |\psi(x)|\,e^{i\phi(x)}$$

- 概率密度：$\rho(x) = |\psi(x)|^2$
- 相位：$\phi(x) = \arg[\psi(x)]$

### 5.2 归一化

$$\|\psi\|^2 = \int_{-\infty}^{\infty} |\psi(x)|^2\,dx = 1$$

数值积分使用梯形法则（$\texttt{np.trapezoid}$）：

$$\int_a^b f(x)\,dx \approx \Delta x\left[\frac{f_0 + f_{N-1}}{2} + \sum_{j=1}^{N-2} f_j\right]$$

### 5.3 内积

$$\langle\phi|\psi\rangle = \int \phi^*(x)\,\psi(x)\,dx$$

### 5.4 期望值

位置期望值：
$$\boxed{\langle x\rangle = \int x\,|\psi(x)|^2\,dx}$$

动量期望值（动量空间计算）：
$$\boxed{\langle p\rangle = \frac{\int k\,|\tilde{\psi}(k)|^2\,dk}{\int |\tilde{\psi}(k)|^2\,dk}}$$

更高阶矩：
$$\langle x^2\rangle = \int x^2\,|\psi|^2\,dx, \quad \langle p^2\rangle = \frac{\int k^2\,|\tilde{\psi}|^2\,dk}{\int |\tilde{\psi}|^2\,dk}$$

### 5.5 不确定度

$$\boxed{\Delta x = \sqrt{\langle x^2\rangle - \langle x\rangle^2}, \quad \Delta p = \sqrt{\langle p^2\rangle - \langle p\rangle^2}}$$

海森堡不确定度原理：
$$\boxed{\Delta x \cdot \Delta p \geq \frac{\hbar}{2}}$$

高斯波包饱和此下界（最小不确定态）。

### 5.6 能量期望值

动能（动量空间计算，避免二阶导数）：
$$\langle T\rangle = \frac{\hbar^2}{2m}\frac{\int k^2\,|\tilde{\psi}(k)|^2\,dk}{\int |\tilde{\psi}(k)|^2\,dk}$$

势能：
$$\langle V\rangle = \int V(x)\,|\psi(x)|^2\,dx$$

总能量：
$$E = \langle H\rangle = \langle T\rangle + \langle V\rangle$$

---

## 6. 矩阵力学

### 6.1 数态表象 (Fock Basis)

选择谐振子数态 $|n\rangle$ 作为基：

$$\boxed{\begin{aligned}
\hat{a}|n\rangle &= \sqrt{n}\,|n-1\rangle \\
\hat{a}^\dagger|n\rangle &= \sqrt{n+1}\,|n+1\rangle \\
\hat{N}|n\rangle &= n|n\rangle
\end{aligned}}$$

矩阵元（N 维截断）：
$$\boxed{\begin{aligned}
\langle n|\hat{a}|m\rangle &= \sqrt{m}\,\delta_{n,m-1} \\
\langle n|\hat{a}^\dagger|m\rangle &= \sqrt{m+1}\,\delta_{n,m+1}
\end{aligned}}$$

### 6.2 坐标与动量算符

$$\boxed{\begin{aligned}
\hat{x} &= \sqrt{\frac{\hbar}{2m\omega}}\,(\hat{a} + \hat{a}^\dagger) \\
\hat{p} &= i\sqrt{\frac{m\hbar\omega}{2}}\,(\hat{a}^\dagger - \hat{a})
\end{aligned}}$$

其中 $\omega$ 是参考频率，$x_0 = \sqrt{\hbar/m\omega}$ 是特征长度。

### 6.3 正则对易关系

$$\boxed{[\hat{x}, \hat{p}] = i\hbar\hat{I}, \quad [\hat{a}, \hat{a}^\dagger] = \hat{I}}$$

> **注意**: 在 N 维截断基中，$[\hat{x}, \hat{p}] = i\hbar\hat{I}$ 仅在低能子空间 ($n \ll N$) 近似成立。最后一个基态因截断而偏离，程序测试中在前 $N-5$ 维子空间验证此关系。

### 6.4 谐振子哈密顿量

$$\boxed{\hat{H}_{\text{HO}} = \hbar\omega\left(\hat{a}^\dagger\hat{a} + \frac{1}{2}\right) = \hbar\omega\left(\hat{N} + \frac{1}{2}\right)}$$

解析能级（在完整无穷维空间中）：
$$\boxed{E_n = \hbar\omega\left(n + \frac{1}{2}\right), \quad n = 0, 1, 2, \ldots}$$

### 6.5 一般势的哈密顿量

$$\hat{H} = \frac{\hat{p}^2}{2m} + V(\hat{x})$$

程序通过在坐标基对角化 $\hat{x}$ 构建 $V(\hat{x})$：
1. 对角化 $\hat{x}$ (N×N 矩阵) → 本征值 $x_i$，本征矢 $U$
2. $V(\hat{x}) = U \cdot \text{diag}[V(x_i)] \cdot U^\dagger$
3. $\hat{H} = \hat{p}^2/2m + V(\hat{x})$

### 6.6 时间演化 (矩阵形式)

$$\boxed{|\psi(t)\rangle = e^{-i\hat{H}t/\hbar}\,|\psi(0)\rangle}$$

通过对角化 $\hat{H}$ 计算矩阵指数：
1. 对角化 $\hat{H} = U \cdot \text{diag}[E_i] \cdot U^\dagger$
2. $e^{-i\hat{H}t/\hbar} = U \cdot \text{diag}[e^{-iE_i t/\hbar}] \cdot U^\dagger$

海森堡绘景：
$$\boxed{\hat{O}(t) = e^{i\hat{H}t/\hbar}\,\hat{O}\,e^{-i\hat{H}t/\hbar}}$$

### 6.7 对易子与反对易子

$$\boxed{\begin{aligned}
[\hat{A}, \hat{B}] &= \hat{A}\hat{B} - \hat{B}\hat{A} \\
\{\hat{A}, \hat{B}\} &= \hat{A}\hat{B} + \hat{B}\hat{A}
\end{aligned}}$$

---

## 7. 势函数

### 7.1 无限深势阱

$$V_{\text{well}}(x) = \begin{cases}
0 & |x| < a/2 \\
V_{\max} & |x| \geq a/2
\end{cases}$$

数值实现使用平滑 sigmoid 边界避免无穷大：

$$V(x) = V_{\max}\left[\frac{1}{1+e^{s(x+a/2+\varepsilon)}} + \frac{1}{1+e^{-s(x-a/2-\varepsilon)}}\right]$$

其中 $s = 50/a$ 控制边界锐度，$\varepsilon = 0.001a$ 避免数值溢出。

**解析能级** (在精确无限深极限下)：
$$\boxed{E_n = \frac{\pi^2\hbar^2 n^2}{2ma^2}, \quad n = 1,2,3,\ldots}$$

**解析波函数**：
$$\psi_n(x) = \sqrt{\frac{2}{a}} \times \begin{cases}
\cos\left(\frac{n\pi x}{a}\right) & n\text{ odd} \\
\sin\left(\frac{n\pi x}{a}\right) & n\text{ even}
\end{cases}$$

### 7.2 谐振子

$$\boxed{V_{\text{HO}}(x) = \frac{1}{2}m\omega^2 x^2}$$

特征长度：$a_{\text{HO}} = \sqrt{\hbar/m\omega}$

**解析能级**：
$$\boxed{E_n = \hbar\omega\left(n + \frac{1}{2}\right)}$$

**解析波函数**：
$$\psi_n(x) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4} H_n\!\left(\sqrt{\frac{m\omega}{\hbar}}x\right) e^{-m\omega x^2/2\hbar}$$

其中 $H_n$ 是 Hermite 多项式。

**相干态** (位移真空态，$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$)：
$$\psi_\alpha(x) = \left(\frac{m\omega}{\pi\hbar}\right)^{1/4} \exp\!\left[-\frac{m\omega}{2\hbar}\!\left(x - \sqrt{\frac{2\hbar}{m\omega}}\,\text{Re}(\alpha)\right)^2 + i\sqrt{\frac{2m\omega}{\hbar}}\,\text{Im}(\alpha)\,x\right]$$

相干态在谐振子势中保持形状不变（不弥散），质心做经典谐振动：$\langle x\rangle(t) = x_0\cos(\omega t) + (p_0/m\omega)\sin(\omega t)$。

### 7.3 矩形势垒

$$\boxed{V_{\text{barrier}}(x) = \begin{cases}
V_0 & |x| < w/2 \\
0 & |x| \geq w/2
\end{cases}}$$

$V_0 > 0$ 为势垒，$V_0 < 0$ 为势阱。

**经典转折点**: 若 $E < V_0$，经典粒子无法穿越。量子隧穿允许概率传输。

**WKB 隧穿概率近似** (矩形势垒，$E < V_0$，$\kappa w \gg 1$)：
$$\boxed{T_{\text{WKB}} \approx \exp\!\left[-2\kappa w\right], \quad \kappa = \frac{\sqrt{2m(V_0 - E)}}{\hbar}}$$

程序中的势垒 demo 展示了超越 WKB 近似的精确数值隧穿。

### 7.4 双势阱 (Quartic)

$$\boxed{V_{\text{DW}}(x) = V_0\left[\left(\frac{2x}{a}\right)^2 - 1\right]^2 - V_0}$$

极小值在 $x = \pm a/2$，阱底 $V_{\min} = -V_0$，中心势垒 $V(0) = 0$。

**隧穿劈裂**：
$$\boxed{\Delta E = E_1 - E_0}$$

对称基态 $|g\rangle = (|L\rangle + |R\rangle)/\sqrt{2}$ 和反对称第一激发态 $|u\rangle = (|L\rangle - |R\rangle)/\sqrt{2}$ 之间的能差。

**隧穿周期**：
$$\boxed{T_{\text{tunnel}} = \frac{2\pi\hbar}{\Delta E}}$$

初始局域在左阱的态会在 $T_{\text{tunnel}}/2$ 时间后隧穿到右阱。

### 7.5 Morse 势

$$\boxed{V_{\text{Morse}}(x) = D_e\left[1 - e^{-\alpha(x-x_0)}\right]^2}$$

用于模拟双原子分子的非谐振动。

**解析能级**：
$$\boxed{E_n = \hbar\omega_0\!\left(n + \frac{1}{2}\right) - \frac{\left[\hbar\omega_0(n + \frac{1}{2})\right]^2}{4D_e}}$$

其中 $\omega_0 = \alpha\sqrt{2D_e/m}$ 是等效谐振子频率。

最大束缚态数：$n_{\max} = \lfloor 2D_e/\hbar\omega_0 - 1/2 \rfloor$

### 7.6 软核库仑势 (1D)

$$\boxed{V_{\text{Coulomb}}(x) = -\frac{Z}{\sqrt{x^2 + a^2}}}$$

软化参数 $a > 0$ 避免 $x=0$ 处的奇点。$Z$ 是有效核电荷。

> **注意**: 这是 1D 模型系统，不等价于 3D 氢原子的径向方程。3D 氢原子的解析能级为 $E_n = -Z^2/(2n^2)$，1D 模型只有有限个束缚态。

### 7.7 周期势 (Kronig-Penney 型)

$$\boxed{V_{\text{periodic}}(x) = A\cos(kx)}$$

根据 Bloch 定理产生能带结构。$A$ 控制带隙宽度。

### 7.8 阶梯势

$$\boxed{V_{\text{step}}(x) = V_0\,\Theta(x - x_0) = \begin{cases}
0 & x < x_0 \\
V_0 & x > x_0
\end{cases}}$$

用于散射问题：反射系数 $R$ 和透射系数 $T$ 取决于入射能量 $E$ 与 $V_0$ 的关系。

### 7.9 自由粒子

$$V_{\text{free}}(x) = 0$$

解析解：平面波 $\psi_k(x) = e^{ikx}/\sqrt{2\pi}$，连续谱 $E = \hbar^2 k^2/2m$。

高斯波包的自由演化：
$$\psi(x,t) = \frac{1}{(2\pi\sigma_t^2)^{1/4}}\exp\!\left[-\frac{(x - p_0 t/m)^2}{4\sigma_0\sigma_t} + i\frac{p_0}{\hbar}x - i\frac{p_0^2}{2m\hbar}t\right]$$

其中 $\sigma_t = \sigma_0(1 + i\hbar t/2m\sigma_0^2)$，波包宽度 $\Delta x(t) = \sigma_0\sqrt{1 + (\hbar t/2m\sigma_0^2)^2}$ 随时间增长（量子弥散）。

---

## 8. 本征值问题

### 8.1 定态薛定谔方程

$$\hat{H}\psi_n(x) = E_n\psi_n(x)$$

对于不含时势，$\psi_n(x,t) = \psi_n(x)e^{-iE_n t/\hbar}$。

### 8.2 数值对角化

程序使用两种方法：

**三对角对角化** (适用于有限差分哈密顿量)：
使用 `scipy.linalg.eigh_tridiagonal`，仅需对角线 $d_j$ 和次对角线 $e_j$：

$$\hat{H} = \begin{pmatrix}
d_0 & e_0 & 0 & \cdots \\
e_0 & d_1 & e_1 & \cdots \\
0 & e_1 & d_2 & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}$$

复杂度 $\mathcal{O}(N^2)$，但仅需要 $\mathcal{O}(N)$ 存储。

**满矩阵对角化** (用于一般势在数态基中的哈密顿量)：
使用 `numpy.linalg.eigh`（厄米矩阵），复杂度 $\mathcal{O}(N^3)$。

### 8.3 虚时间演化 (基态搜索)

对于任意初始态 $|\psi_0\rangle$（要求 $\langle\psi_0|0\rangle \neq 0$），虚时间演化投影到基态：

$$\boxed{|\psi(\tau)\rangle = \frac{e^{-\hat{H}\tau}|\psi_0\rangle}{\|e^{-\hat{H}\tau}|\psi_0\rangle\|} \xrightarrow{\tau\to\infty} |0\rangle}$$

在能量本征基中展开：
$$|\psi_0\rangle = \sum_n c_n|n\rangle \implies |\psi(\tau)\rangle \propto \sum_n c_n e^{-E_n\tau}|n\rangle$$

高能分量被指数压制，$\tau$ 足够大时只剩基态。程序实现：
$$\text{for } k = 1,\ldots,K:\; |\psi\rangle \leftarrow e^{-\hat{H}\Delta\tau}|\psi\rangle,\; |\psi\rangle \leftarrow |\psi\rangle/\||\psi\rangle\|$$

其中 $\Delta\tau = \tau_{\text{total}} / K$。

---

## 9. 数值稳定性与误差分析

### 9.1 SSFM 误差

每步局部截断误差：$\mathcal{O}(\Delta t^3)$

全局误差（$N_t = T/\Delta t$ 步）：$\mathcal{O}(\Delta t^2)$

优势来源：动能算符在动量空间严格对角（谱精度），势能算符在坐标空间严格对角。

### 9.2 CN 误差

时间离散：$\mathcal{O}(\Delta t^2)$
空间离散（三点差分）：$\mathcal{O}(\Delta x^2)$

无条件稳定意味着 $\Delta t$ 可以大于 CFL 条件限制，但精度要求 $\Delta t$ 和 $\Delta x$ 适当匹配。

### 9.3 CFL 条件 (SSFM 的稳定性约束)

虽然 SSFM 不是严格的 CFL 限制，但建议：

$$\Delta t \lesssim \frac{\Delta x}{|v_{\max}|}$$

其中 $v_{\max} = k_{\max}/m = \pi/(m\Delta x)$ 是最快动量分量对应的速度。

### 9.4 网格分辨率要求

**空间分辨率**：
- 势的特征尺度：$\Delta x \ll \min(\text{potential width features})$
- 波函数特征尺度：$\Delta x \ll \sigma$（波包宽度）
- Nyquist 条件：$\Delta x < \pi/k_{\max}$

**网格范围**：
- 必须包含波函数的完整支集
- 边界反射应可忽略（$\psi \approx 0$ at boundaries）
- 对束缚态，$x_{\max} \gg$ 经典转折点

### 9.5 范数与能量守恒

程序监视两个守恒量：

**范数守恒**：
$$\delta_N = \max_t \big|1 - \|\psi(t)\|^2\big|$$

SSFM 的范数精确守恒（幺正性），仅受浮点精度影响。CN 的范数 drift 与截断相关。

**能量守恒**（不含时势）：
$$\delta_E = \frac{\max_t E(t) - \min_t E(t)}{\bar{E}}$$

SSFM 能量 drift 来源：
- Trotter 分解误差 ($\propto \Delta t^2$)
- 势的梯度大时，Trotter 分解精度下降
- 谐振子外围（$\langle V\rangle$ 大）比阱底误差更大

### 9.6 实践建议

| 物理场景 | 推荐方法 | 推荐 $\Delta t$ | 推荐 $N$ |
|----------|----------|:---:|:---:|
| 自由粒子 | SSFM | 0.005–0.02 | 512–2048 |
| 谐振子 | SSFM | 0.001–0.01 | 512–1024 |
| 势垒隧穿 | SSFM | 0.001–0.005 | 1024–4096 |
| 双势阱 | SSFM 或 CN | 0.001–0.005 | 1024–2048 |
| 刚性势 (尖锐边界) | CN (精确守恒) | 0.01–0.05 | 256–512 |
| 长时间 (数百周期) | CN (无条件稳定) | 0.01–0.05 | 256–512 |

---

## 附录：符号表

| 符号 | 含义 | 程序中的默认值 |
|------|------|:---:|
| $\hbar$ | 约化普朗克常数 | 1.0 (a.u.) |
| $m$ | 粒子质量 | 1.0 (a.u.) |
| $\Delta t$ | 时间步长 | 0.001 |
| $\Delta x$ | 空间步长 | $(x_{\max}-x_{\min})/(N-1)$ |
| $N$ | 网格点数 | 1024 |
| $x_{\min}, x_{\max}$ | 空间范围 | [-10, 10] |
| $N$ (矩阵) | 数态基截断 | 50 |
| $\omega$ | 谐振子频率 | 1.0 |
| $a_{\text{HO}}$ | 特征长度 | $\sqrt{\hbar/m\omega}$ |
| $\alpha$ | CN 常数 | $\hbar^2/(2m\Delta x^2)$ |

# 量子力学物理基础

## 含时薛定谔方程 (TDSE)

一维含时薛定谔方程:

$$i\hbar\frac{\partial}{\partial t}\psi(x,t) = \left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi(x,t)$$

或写成:

$$i\hbar\frac{\partial\psi}{\partial t} = \hat{H}\psi$$

其中 $\hat{H} = \hat{T} + \hat{V}$ 是哈密顿量。

## 数值方法

### Split-Step Fourier Method (SSFM)

基于 Trotter-Suzuki 分解:

$$e^{-i\hat{H}\Delta t/\hbar} = e^{-i\hat{V}\Delta t/2\hbar} e^{-i\hat{T}\Delta t/\hbar} e^{-i\hat{V}\Delta t/2\hbar} + \mathcal{O}(\Delta t^3)$$

算法:
1. 半步势能演化 (坐标空间): $\psi_1 = e^{-iV(x)\Delta t/2\hbar}\psi$
2. FFT 到动量空间: $\tilde{\psi}_1 = \mathcal{F}[\psi_1]$
3. 动能演化 (动量空间): $\tilde{\psi}_2 = e^{-i\hbar k^2\Delta t/2m}\tilde{\psi}_1$
4. 逆 FFT: $\psi_2 = \mathcal{F}^{-1}[\tilde{\psi}_2]$
5. 半步势能演化: $\psi(t+\Delta t) = e^{-iV(x)\Delta t/2\hbar}\psi_2$

**优势**: 谱精度、O(N log N) 复杂度
**限制**: 要求光滑势函数、周期边界

### Crank-Nicolson Method (CN)

隐式离散化:

$$\frac{i\hbar}{\Delta t}(\psi^{n+1} - \psi^n) = \frac{1}{2}\hat{H}(\psi^{n+1} + \psi^n)$$

重写为:

$$(I + \frac{i\Delta t}{2\hbar}\hat{H})\psi^{n+1} = (I - \frac{i\Delta t}{2\hbar}\hat{H})\psi^n$$

使用三对角近似 + Thomas 算法求解。

**优势**: 无条件稳定、精确守恒
**限制**: 需要解矩阵方程 (但 Thomas 算法 O(N))

## 矩阵力学

### 正则对易关系

$$[\hat{x}, \hat{p}] = i\hbar$$

$$[\hat{a}, \hat{a}^\dagger] = 1$$

### 数态表象 (Fock basis)

$$\hat{a}|n\rangle = \sqrt{n}|n-1\rangle$$
$$\hat{a}^\dagger|n\rangle = \sqrt{n+1}|n+1\rangle$$
$$\hat{N}|n\rangle = n|n\rangle$$

坐标和动量:

$$\hat{x} = \sqrt{\frac{\hbar}{2m\omega}}(\hat{a} + \hat{a}^\dagger)$$
$$\hat{p} = i\sqrt{\frac{m\hbar\omega}{2}}(\hat{a}^\dagger - \hat{a})$$

### 谐振子

$$\hat{H} = \hbar\omega(\hat{a}^\dagger\hat{a} + \frac{1}{2})$$

能级: $E_n = \hbar\omega(n + \frac{1}{2})$, $n = 0,1,2,...$

基态波函数: $\psi_0(x) = (\frac{m\omega}{\pi\hbar})^{1/4} e^{-m\omega x^2/2\hbar}$

不确定度: $\Delta x \cdot \Delta p = \hbar/2$ (最小不确定态)

### 相干态

位移真空态: $|\alpha\rangle = e^{-|\alpha|^2/2}\sum_n \frac{\alpha^n}{\sqrt{n!}}|n\rangle$

在谐振子势中保持形状不变，质心做经典运动 (Ehrenfest 定理)。

## 势函数

### 无限深势阱

$$V(x) = \begin{cases} 0 & |x| < a/2 \\ \infty & |x| \geq a/2 \end{cases}$$

解析能级: $E_n = \frac{\pi^2\hbar^2 n^2}{2ma^2}$, $n = 1,2,3,...$

### 谐振子

$$V(x) = \frac{1}{2}m\omega^2 x^2$$

解析能级: $E_n = \hbar\omega(n + 1/2)$

### 势垒隧穿

$$V(x) = \begin{cases} V_0 & |x| < w/2 \\ 0 & \text{elsewhere} \end{cases}$$

WKB 隧穿概率: $T \approx \exp(-2\int_{x_1}^{x_2} \sqrt{2m(V(x)-E)}/\hbar \,dx)$

对于矩形势垒 ($E < V_0$): $T \approx e^{-2\kappa w}$, $\kappa = \sqrt{2m(V_0-E)}/\hbar$

### 双势阱 (Quartic)

$$V(x) = V_0\left[\left(\frac{2x}{a}\right)^2 - 1\right]^2 - V_0$$

极小值在 $x = \pm a/2$, 阱底 $V = -V_0$, 中心势垒 $V(0) = 0$。

隧穿劈裂 $\Delta E = E_1 - E_0$ 决定隧穿周期 $T = 2\pi\hbar/\Delta E$。

### Morse 势 (双原子分子)

$$V(x) = D_e\left[1 - e^{-\alpha(x-x_0)}\right]^2$$

解析能级: $E_n = \hbar\omega_0(n + \frac{1}{2}) - \frac{[\hbar\omega_0(n + \frac{1}{2})]^2}{4D_e}$

其中 $\omega_0 = \alpha\sqrt{2D_e/m}$。

## 不确定度原理

$$\Delta x \cdot \Delta p \geq \frac{\hbar}{2}$$

高斯波包饱和下界 (最小不确定态)。

## 参考

- Sakurai, *Modern Quantum Mechanics*
- Tannor, *Introduction to Quantum Mechanics: A Time-Dependent Perspective*
- Feit, Fleck, & Steiger, *J. Comput. Phys.* 47, 412 (1982) — SSFM 方法

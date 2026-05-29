# Quantum Agent 学习路径

> 边仿真边学量子力学——从波函数到量子场论

---

## 总览：四阶段递进

```
第一阶段 ──→ 第二阶段 ──→ 第三阶段 ──→ 第四阶段
波函数图像    Fock 空间      自旋 & 纠缠   量子场论
(可视化)     (量子光学)     (量子信息)    (QFT 基础)
```

每阶段：跑 demo → 改参数 → 验证理论 → 自己写 .qms。

---

## 第一阶段：波函数图像（1-2 天）

**核心问题**：粒子的波函数如何随时间演化？

### 实验 1：自由粒子量子弥散

```bash
python agent.py --run scripts/free_particle.qms
```

打开 `scripts/output/animations/free_particle_spreading.mp4`。

**观察**：波包随时间展宽。为什么？不确定性原理——初始位置较确定 → 动量不确定 → 不同动量分量以不同速度传播。

**动手改**：编辑 `scripts/free_particle.qms`，改 `sigma=0.5` vs `sigma=3.0`。弥散速度差多少倍？对照 `τ = 2mσ²/ℏ`。

### 实验 2：不确定性原理

```bash
python agent.py --run scripts/heisenberg_uncertainty.qms
```

4 面板动画：位置空间、动量空间、Δx(t) & Δp(t)、Δx·Δp ≥ ℏ/2。

**观察**：t=0 时 Δx·Δp = ℏ/2（最小不确定态），随后持续增大。"FORBIDDEN" 区域永远不可达。

### 实验 3：测量坍缩

```bash
python agent.py --run scripts/measurement_collapse.qms
```

**观察**：宽波包 → 测量 → 坍缩为窄波包 → 坍缩后弥散加速 100×。

**动手**：改 `sigma_narrow=0.1`，弥散加速倍数变为多少？

### 实验 4：双阱势隧穿

```bash
python agent.py --run scripts/double_well.qms
```

**观察**：粒子在两个势阱之间振荡——经典粒子需要翻越势垒，量子的波函数直接隧穿。

### 第一阶段 Check

```python
⚛ > grid = WaveGrid(-20, 20, 1024)
⚛ > psi0 = gaussian_wavepacket(grid, x0=-5, p0=3, sigma=1)
⚛ > res = evolve_ssfm(psi0, grid, dt=0.01, t_max=10, snapshots=200)
⚛ > animate res output/my_wave.mp4
```

---

## 第二阶段：Fock 空间与量子光学（2-3 天）

**核心问题**：光子是什么？

### 实验 5：Fock 空间初探

```python
⚛ > demo      # 跑内建演示
⚛ > fb = FockBasis(30)
⚛ > fb.a[:4,:4]    # 湮灭算符
⚛ > fb.a_dag[:4,:4] # 产生算符
```

**关键**：`a|n⟩ = √n|n-1⟩`，`a†|n⟩ = √(n+1)|n+1⟩`。Fock 态 |n⟩ 表示 n 个光子/量子激发。

### 实验 6：对易子验证

```python
⚛ > C = commutator(fb.x, fb.p)
⚛ > np.linalg.norm(C[:10,:10] - 1j*np.eye(10), 'fro')
# → ~4e-15  [x̂,p̂] = iħ ✓
⚛ > help commutator
```

### 实验 7：五种量子态的 Wigner 函数

```bash
python agent.py --run scripts/wigner_gallery.qms
```

| 态 | Wigner 特征 | g²(0) | 经典/量子 |
|---|------|:---:|:---:|
| 真空 |0⟩ | 原点高斯 | — | |
| 相干态 |α⟩ | 位移高斯，全正 | 1.0 | 最经典 |
| 热态 ρ_th | 展宽高斯，全正 | 2.0 | 经典 |
| 压缩真空 |ζ⟩ | 压扁椭圆 | >1 | 量子 |
| 猫态 | 双峰+负值干涉条纹 | — | 纯量子 |

**核心洞察**：`W < 0 = 非经典性的标志`。

### 实验 8：光子统计

```python
⚛ > psi = coherent(30, 2.0)
⚛ > g2(psi)              # 1.0 — Poisson
⚛ > help g2              # g²(0) = ⟨a†a†aa⟩/⟨a†a⟩²
```

### 第二阶段 Check

能解释：Fock 态 |n⟩、产生/湮灭算符、相干态为什么是 Poisson 统计、Wigner 负值意味着什么。

---

## 第三阶段：自旋与纠缠（1-2 天）

**核心问题**：两体量子系统有什么新现象？

### 实验 9：Pauli 代数

```python
⚛ > PauliX
⚛ > sx, sy, sz = PauliX @ PauliY, PauliZ
⚛ > sx @ sy - sy @ sx    # [σₓ, σ_y] = 2i σ_z
```

### 实验 10：Bloch 球

```python
⚛ > psi = qubit(np.pi/3, np.pi/4)     # θ=60°, φ=45°
⚛ > rho = bloch_dm(0.6, 0, 0.8)       # 混合态
⚛ > bloch_length(rho)                  # 纯态=1, 混合<1
```

### 实验 11：Bell 态与并发度

```python
⚛ > from src.qm import bell_states, concurrence, entropy_vn
⚛ > bells = bell_states()
⚛ > concurrence(bells['Φ⁺'])          # 1.0 — 最大纠缠
⚛ > concurrence(np.array([1,0,0,0]))  # 0.0 — 可分离
```

**关键**：纠缠不是"一个粒子立刻影响另一个"——而是"两个粒子构成一个不可分的整体"。

### 实验 12：张量积与偏迹

```python
⚛ > rho_A = np.array([[0.6,0],[0,0.4]])
⚛ > rho_B = np.array([[0.5,0],[0,0.5]])
⚛ > rho_AB = tensor(rho_A, rho_B)
⚛ > partial_trace(rho_AB, dims=(2,2), keep=0)  # 还原 rho_A
```

---

## 第四阶段：量子场论入门（3-5 天）

**核心问题**：场是更基本的——粒子是场的激发。

### 实验 13：标量场对易子

```python
⚛ > sf = ScalarField(mass=1.0)
⚛ > sf.commutator(0, 2)          # [φ̂(x), π̂(y)] = iδ(x-y)
⚛ > help scalarfield
```

### 实验 14：格点 φ⁴ 理论

```python
⚛ > lat = LatticePhi4(N_sites=8, mass=1.0, coupling=0.5)
⚛ > E0, psi0 = lat.ground_state()
⚛ > E0                              # 基态能量
```

### 实验 15：路径积分 Monte Carlo

```bash
python agent.py --run scripts/pimc_demo.qms
```

**核心**：PIMC 验证谐振子基态 E₀=0.5。量子力学 = 经典统计力学(虚时间)。

### 实验 16：Feynman 振幅

```python
⚛ > from src.qft import feynman_amplitude_phi4_2to2
⚛ > amp = feynman_amplitude_phi4_2to2(s=10, t=-2, u=-8, coupling=1.0)
```

---

## 学习技巧

### 1. 善用 `help <function>`

```python
⚛ > help coherent    # 公式 + 物理解释
⚛ > help g2          # 24 个函数有公式
```

### 2. 读完知识手册对应章节

每个实验对应 `docs/KNOWLEDGE_HANDBOOK.md` 的一卷：

| 实验阶段 | 知识手册卷 |
|---------|:---:|
| 第一阶段 (波函数) | — |
| 第二阶段 (Fock) | 第〇卷：为什么需要 Fock 空间 |
| 第二阶段 (Wigner) | 第三卷：为什么需要 Wigner 函数 |
| 第三阶段 (自旋) | — |
| 第四阶段 (QFT) | 第一卷：闵氏空间与 Fock 空间的缝合 |
| (进阶) 因果 | 第二卷：因果关系在 QM 中的表现 |
| (进阶) 退相干 | 第四卷：退相干的机制 |

### 3. 动笔——改写 .qms 脚本

每个脚本都是学习的起点。改参数、改势函数、加新的可观测量。从改数字开始，到写新的分析逻辑。

```python
# 从 harmonic_oscillator.qms 开始
# 改 alpha、改 Fock 截断、加新的测量
⚛ > run scripts/harmonic_oscillator.qms
# 打开文件编辑 → 再运行 → 观察变化
```

### 4. 查文档

```bash
# 完整函数参考
less docs/USER_GUIDE.md

# 8 步入门教程
less docs/TUTORIAL.md

# 物理基础索引
less docs/PHYSICS.md
```

---

## 时间表建议

| 天数 | 内容 | 关键动作 |
|:---:|------|------|
| 1 | 自由弥散 + 不确定性原理 | 跑 demo, 改 σ, 看动画 |
| 2 | 测量坍缩 + 双阱隧穿 | 写第一个自己的势函数 |
| 3 | Fock 空间 + 对易子 | `demo`, `fb.a`, `commutator` |
| 4 | 五种态 + Wigner | `wigner_gallery.qms`, 理解负值 |
| 5 | 光子统计 g² | 对比 coherent/thermal/fock |
| 6 | Pauli 代数 + Bloch 球 | `PauliX`, `qubit`, `bloch_dm` |
| 7 | Bell 态 + 并发度 | `bell_states`, `concurrence` |
| 8-10 | 标量场 → φ⁴ → PIMC → Feynman | 知手卷第一卷, `pimc_demo.qms` |

---

## 黄金原则

```
跑 demo → 改参数 → 验证理论 → 自己写脚本
   ↓          ↓          ↓           ↓
  看          问         算          造
```

每一步不跳过。看完动画问"为什么"，算完数字对理论，最后造出自己的分析。

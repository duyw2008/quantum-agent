# Quantum Agent 教程

> 30 分钟从零到第一个量子动画

---

## 第 1 步：启动 Agent

```bash
cd ~/quantum_agent
python agent.py
```

```
⚛ ~/quantum_agent > 
```

直接输入 Python 表达式，函数已预加载。

---

## 第 2 步：第一个量子态

```python
# 相干态
⚛ > psi = coherent(30, 2.0)
⚛ > mean_photon(psi)
4.0                    # ⟨n⟩ = |α|²

# 查公式
⚛ > help coherent
  ═══ coherent ═══
  |α⟩ = e-|α|²/2 ∑ₙ (αⁿ)/(√(n!)) |n⟩
```

---

## 第 3 步：对易子验证

```python
⚛ > C = commutator(fb.x, fb.p)
⚛ > np.linalg.norm(C[:10,:10] - 1j*np.eye(10), 'fro')
4.2e-15                # [x̂,p̂] = iħ ✓
```

---

## 第 4 步：光子统计

```python
⚛ > g2(psi)             # 相干态 Poisson
1.0

⚛ > rho = thermal_dm(30, 1.5)
⚛ > g2(rho)             # 热态聚束
2.0
```

---

## 第 5 步：Wigner 函数

```python
⚛ > x, p, W = wigner(psi, N_grid=61)
⚛ > plot_wigner(x, p, W)
⚛ > W.min()             # 相干态 W>0
~0
```

试试猫态——有负值：
```python
⚛ > cat_even = cat(30, 2.0, 0)
⚛ > x2, p2, W2 = wigner(cat_even)
⚛ > W2.min()            # 负值 = 非经典！
-0.43
```

---

## 第 6 步：波函数动画

```python
⚛ > grid = WaveGrid(-30, 30, 1024)
⚛ > psi0 = gaussian_wavepacket(grid, x0=-8, p0=2.0, sigma=1.5)
⚛ > res = evolve_ssfm(psi0, grid, dt=0.01, t_max=8.0, snapshots=200)
⚛ > animate res output/wave.mp4
```

打开 `output/wave.mp4` — 波包弥散动画。

---

## 第 7 步：复杂势函数

用 PotentialBuilder 链式组合：

```python
⚛ > V = (PotentialBuilder(grid)
...      .harmonic(0.3)
...      .barrier(0, 5, 1.5)
...      .gaussian(-5, -3, 1.5)
...      .gaussian(5, -3, 1.5)
...      .build())
⚛ > V.plot()           # 预览
⚛ > V.summary()         # 组件清单
```

---

## 第 8 步：运行脚本

已有 10 个预制脚本：

```bash
python agent.py --run scripts/heisenberg_uncertainty.qms
```

```
⚛ > run scripts/double_well.qms      # Tab 补全可用
```

---

## 下一步

- 读 [USER_GUIDE.md](USER_GUIDE.md) 完整函数参考
- 读 [KNOWLEDGE_HANDBOOK.md](KNOWLEDGE_HANDBOOK.md) 五卷深度物理
- 写自己的 `.qms` 脚本

```
⚛ > help           命令列表
⚛ > help wigner    24 个函数的物理公式
⚛ > demo           Fock 基演示
⚛ > test           自检
```

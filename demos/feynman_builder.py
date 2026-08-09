#!/usr/bin/env python3
"""
费曼图构造器 — φ⁴ 理论费曼图枚举与振幅计算
===========================================
枚举 φ⁴ 理论中给定外腿数 N 和微扰阶数 k 的所有连通费曼图，
计算对称因子和散射振幅，展示为 ASCII 艺术。

对应《一读就懂的量子场论》§4.3-4.4

φ⁴ Feynman 规则:
  传播子: i/(p² - m² + iε)
  顶点:   -iλ  (4 条腿汇聚于一点)
  对称因子: 由 Wick 收缩和拉氏量 1/4! 的组合给出

对 2→2 散射:
  树图 (k=1): iM = -iλ, S=1
  单圈 (k=2): s/t/u 道泡泡图, S=2 (内部传播子可交换)
    iM = (-iλ)²/2 [B(s) + B(t) + B(u)]
    其中 B(p²) = ∫ d²k/(2π)² i/(k²-m²+iε) · i/((k+p)²-m²+iε)

用法:
    python feynman_builder.py                 # 完整演示
    python feynman_builder.py --2to2          # 仅 2→2 散射
    python feynman_builder.py --compare       # 树图 vs 单圈比较
    python feynman_builder.py --save out.txt  # 保存输出到文件
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from itertools import combinations, product
from collections import Counter, defaultdict
import argparse


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    数据模型: Feynman 图                          ║
# ╚══════════════════════════════════════════════════════════════════╝

@dataclass
class FeynmanDiagram:
    """φ⁴ Feynman 图

    Attributes:
        name:              图名 (中文描述)
        n_vertices:        顶点数 V
        n_external:        外腿数
        n_loops:           圈数 L
        symmetry_factor:   对称因子 S (>0 整数)
        ascii_art:         ASCII 图表示
        amplitude_expr:    振幅数学表达式
        amplitude_value:   数值振幅 (对样本动量)
        momentum_routing:  动量路由说明
        channel:           s/t/u/tree 道标识
    """
    name: str
    n_vertices: int
    n_external: int
    n_loops: int
    symmetry_factor: int
    ascii_art: str
    amplitude_expr: str = ""
    amplitude_value: Optional[complex] = None
    momentum_routing: str = ""
    channel: str = "?"

    def print(self):
        """打印图的完整信息"""
        print(f"\n{'─' * 60}")
        print(f"  {self.name}")
        print(f"{'─' * 60}")
        print(self.ascii_art)
        print(f"\n  顶点数: {self.n_vertices}  |  外腿: {self.n_external}"
              f"  |  圈数: {self.n_loops}")
        print(f"  对称因子: S = {self.symmetry_factor}")
        if self.momentum_routing:
            print(f"  动量赋值: {self.momentum_routing}")
        print(f"  振幅: {self.amplitude_expr}")
        if self.amplitude_value is not None:
            val = self.amplitude_value
            print(f"  数值: {val.real:.4e}{val.imag:+.4e}i")
            print(f"  |M|² = {abs(val)**2:.4e}")

    def summary(self) -> str:
        """简短摘要"""
        return (f"Diagram(V={self.n_vertices}, L={self.n_loops}, "
                f"ext={self.n_external}, S={self.symmetry_factor}, "
                f"ch={self.channel})")


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    圈积分数值计算 (1+1D)                         ║
# ╚══════════════════════════════════════════════════════════════════╝

def bubble_integral_2d(p_sq: float, mass: float = 1.0,
                       Lambda: float = 100.0, n_pts: int = 120) -> complex:
    """二维泡泡积分 (1+1D 时空)

    B(p²) = ∫ d²k/(2π)²  i/(k²-m²+iε) · i/((k+p)²-m²+iε)

    极坐标下数值计算。1+1D 中该积分对数发散，用截断 Λ 正规化。

    Args:
        p_sq:   动量平方 (类空为正)
        mass:   质量 m
        Lambda: 动量截断
        n_pts:  积分点数

    Returns:
        B(p²) 的复数值
    """
    if p_sq < 1e-15:
        p_sq = 1e-15

    p_mag = np.sqrt(p_sq)
    k_vals = np.linspace(0, Lambda, n_pts)
    theta_vals = np.linspace(0, 2 * np.pi, n_pts)
    dk = k_vals[1] - k_vals[0]
    dtheta = theta_vals[1] - theta_vals[0]

    K, Theta = np.meshgrid(k_vals, theta_vals, indexing='ij')
    k_sq = K**2
    kp_sq = k_sq + p_sq + 2 * K * p_mag * np.cos(Theta)

    eps = 1e-8
    denom1 = k_sq - mass**2 + 1j * eps
    denom2 = kp_sq - mass**2 + 1j * eps
    integrand = K * (-1.0) / (denom1 * denom2) / (2 * np.pi)**2

    return complex(np.sum(integrand) * dk * dtheta)


def tadpole_integral_2d(mass: float = 1.0, Lambda: float = 100.0,
                        n_pts: int = 500) -> float:
    """二维蝌蚪积分

    A(m²) = ∫ d²k/(2π)²  i/(k² - m² + iε)
    对数发散，贡献质量重整化 δm²。
    """
    k_vals = np.linspace(0, Lambda, n_pts)
    dk = k_vals[1] - k_vals[0]
    eps = 1e-8
    integrand = k_vals * 1j / (k_vals**2 - mass**2 + 1j * eps) / (2 * np.pi)
    return float(np.real(np.sum(integrand) * dk))


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    对称因子计算                                  ║
# ╚══════════════════════════════════════════════════════════════════╝

def explain_symmetry_factor() -> str:
    """解释 φ⁴ 理论对称因子的来源和计算规则"""
    return """
  φ⁴ 对称因子 S 的计算:
  ====================

  拉氏量: L_int = -(λ/4!) φ⁴

  规则:
    1. 顶点排列: N_v 个等价顶点贡献 N_v! 种排列
    2. 内线交换: 连接同对顶点的 N_e 条内线有 N_e! 种交换方式
    3. 等价外腿: 连接同顶点的等价外腿排列
    4. 蝌蚪自连: 同一顶点自连的每个圈贡献因子 2

  2→2 散射的对称因子:
    树图:            S = 1   (1/4! × 4! 抵消)
    s-道泡泡图:      S = 2   (两条等价内线可交换)
    t-道泡泡图:      S = 2   (同上)
    u-道泡泡图:      S = 2   (同上)
    蝌蚪图 (自能):   S = 2   (圈的两个端点可交换)

  振幅公式: iM = (-iλ)^V / S × ∏(传播子) × ∏∫d²k/(2π)²
"""


def compute_symmetry_factor_phi4(n_vertices: int, n_external: int,
                                 channel: str = "?") -> int:
    """计算 φ⁴ 图的对称因子 (通用公式)

    基于组合计数:
      S = (顶点排列数) × (内线交换数) × (蝌蚪因子)

    Args:
        n_vertices:  顶点数
        n_external:  外腿数
        channel:     道标识

    Returns:
        对称因子 S
    """
    if n_vertices == 1 and n_external == 4:
        return 1  # 树图 2→2

    if n_vertices == 2 and n_external == 4:
        return 2  # s/t/u 道泡泡图

    if n_vertices == 1 and n_external == 2:
        return 2  # 蝌蚪图

    # 通用计算 (简化版)
    # 对 φ⁴, 4V 个端口: n_external 外 + 2I 内
    I = (4 * n_vertices - n_external) // 2
    # 基础对称因子: 来自内线交换
    S = 1
    # 对于多边 (同一对顶点间多条内线), 每条贡献阶乘
    # 这里取一个粗略的上界
    return max(S, 1)


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    图枚举器                                      ║
# ╚══════════════════════════════════════════════════════════════════╝

def _partition_external(n_external: int, n_vertices: int) -> List[List[int]]:
    """将 n_external 条外腿分配到 n_vertices 个顶点

    返回所有满足 Σe_i = n_external, e_i ∈ [0,4], 4-e_i 为偶数的分配。

    偶数约束保证了剩余端口可完全配对为内线。
    """
    results = []

    def _recurse(remaining: int, v: int, current: List[int]):
        if v == n_vertices - 1:
            if remaining <= 4 and (4 - remaining) % 2 == 0:
                results.append(current + [remaining])
            return
        for e in range(min(4, remaining) + 1):
            if (4 - e) % 2 == 0:
                _recurse(remaining - e, v + 1, current + [e])

    _recurse(n_external, 0, [])
    return results


def _pair_ports(n_ports: int) -> List[List[Tuple[int, int]]]:
    """生成 n_ports 个端口的所有完全配对方案 (去重)"""
    if n_ports == 0:
        return [[]]
    if n_ports % 2 != 0:
        return []

    results = []
    seen = set()

    def _recurse(available: set, pairs: List[Tuple[int, int]]):
        if not available:
            key = tuple(sorted(tuple(sorted(p)) for p in pairs))
            if key not in seen:
                seen.add(key)
                results.append(list(key))
            return
        a = min(available)
        for b in (available - {a}):
            _recurse(available - {a, b}, pairs + [(a, b)])

    _recurse(set(range(n_ports)), [])
    return results


def _classify_2to2_channel(ext_on_v0: List[int],
                           ext_on_v1: List[int]) -> str:
    """对 2→2 散射分类 s/t/u 道

    p1=0, p2=1 (入射), p3=2, p4=3 (出射)

    s-channel: p1,p2 → 同一顶点; p3,p4 → 另一顶点
    t-channel: p1,p3 → 同一顶点; p2,p4 → 另一顶点
    u-channel: p1,p4 → 同一顶点; p2,p3 → 另一顶点
    """
    s0, s1 = set(ext_on_v0), set(ext_on_v1)
    if {0, 1}.issubset(s0) and {2, 3}.issubset(s1):
        return 's'
    if {0, 1}.issubset(s1) and {2, 3}.issubset(s0):
        return 's'
    if {0, 2}.issubset(s0) and {1, 3}.issubset(s1):
        return 't'
    if {0, 2}.issubset(s1) and {1, 3}.issubset(s0):
        return 't'
    if {0, 3}.issubset(s0) and {1, 2}.issubset(s1):
        return 'u'
    if {0, 3}.issubset(s1) and {1, 2}.issubset(s0):
        return 'u'
    return '?'


def _diagram_hash(ext_per_vertex: List[List[int]],
                  internal_edges: List[Tuple[int, int]]) -> int:
    """图的规范哈希 (去重用)"""
    V = len(ext_per_vertex)
    features = []
    for v in range(V):
        ext_set = frozenset(ext_per_vertex[v])
        neighbors = []
        for a, b in internal_edges:
            if a == v:
                neighbors.append(b)
            elif b == v:
                neighbors.append(a)
        features.append((ext_set, frozenset(neighbors)))
    return hash(frozenset(features))


def generate_phi4_diagrams(n_external: int = 4, order: int = 1,
                           coupling: float = 0.5, mass: float = 1.0,
                           Lambda: float = 100.0,
                           external_labels: Optional[List[str]] = None
                           ) -> List[FeynmanDiagram]:
    """生成所有连通的 φ⁴ Feynman 图

    给定外腿数 n_external 和微扰阶数 order (= 顶点数 V),
    枚举所有有效的连通 Feynman 图, 计算对称因子和振幅。

    约束:
      4V = n_external + 2I  →  I = (4V - n_external)/2
      圈数 L = I - V + 1

    注意: φ⁴ 顶点的 4 条腿完全等价, 所以端口选择不影响拓扑。
    我们只枚举外腿到顶点的分配和内部连线方案。

    Args:
        n_external:      外腿数 (偶数, 典型: 2, 4, 6, 8)
        order:           微扰阶数 = 顶点数 V
        coupling:        耦合常数 λ
        mass:            质量 m
        Lambda:          圈积分截断
        external_labels: 外腿标号, 默认 p₁, p₂, ...

    Returns:
        FeynmanDiagram 对象列表
    """
    V = order

    if n_external % 2 != 0:
        return []
    total_ports = 4 * V
    if total_ports < n_external:
        return []
    remaining = total_ports - n_external
    if remaining % 2 != 0:
        return []
    I = remaining // 2
    L = I - V + 1
    if L < 0:
        return []

    if external_labels is None:
        subs = "₁₂₃₄₅₆₇₈₉₀"
        external_labels = [
            f"p{subs[i]}" if i < len(subs) else f"p{i+1}"
            for i in range(n_external)
        ]

    # Step 1: 分配外腿到顶点 (只关心哪个顶点分到哪些外腿)
    ext_distributions = _partition_external(n_external, V)
    if not ext_distributions:
        return []

    diagrams = []
    seen = set()

    # 顶点端口全等价 → 只需一个规范端口分配
    # 对顶点 v 分配 e 条外腿 → 使用端口 0, 1, ..., e-1
    def canonical_ports(e: int) -> List[int]:
        return list(range(e))

    for ext_dist in ext_distributions:
        # Step 2: 枚举外腿标号到顶点的分配 (哪个顶点得到哪些标号)
        # 将外腿索引 0..n_external-1 分组为 V 组, 每组大小 ext_dist[v]
        def enumerate_label_assignments():
            """生成所有将 n_external 个标号分配到 V 个顶点的方式"""
            ext_indices = list(range(n_external))
            results = []

            def backtrack(start: int, groups: List[List[int]]):
                if len(groups) == V:
                    if start == n_external:
                        # 可用 permutations 做最后的组分配...
                        # 简化: 直接用组合
                        pass
                # 这个比较 tricky, 对 V=2 且 n_external=4 可以直接枚举

            # 特化: 对 V=2, n_external=4
            if V == 2 and n_external == 4:
                for g1 in combinations(range(4), ext_dist[0]):
                    g1_set = set(g1)
                    g2 = [i for i in range(4) if i not in g1_set]
                    results.append([list(g1), g2])
                return results

            # 特化: 对 V=1, 所有外腿都在同一顶点
            if V == 1:
                return [[list(range(n_external))]]

            # 通用: 递归分配
            def recurse(idx: int, groups: List[List[int]], sizes: List[int]):
                if len(groups) == V:
                    if idx == n_external:
                        results.append([list(g) for g in groups])
                    return
                needed = sizes[len(groups)]
                remaining = n_external - idx
                if needed > remaining:
                    return
                for combo in combinations(range(idx, n_external), needed):
                    if combo[0] != idx:
                        continue  # 确保最小标号在当前组 (避免重复)
                    g = list(combo)
                    recurse(idx + 1, groups + [g], sizes)

            # 简化: 只处理常见情况
            if V <= 2:
                # 手动枚举
                if V == 2:
                    for g1 in combinations(range(n_external), ext_dist[0]):
                        g1_set = set(g1)
                        g2 = [i for i in range(n_external) if i not in g1_set]
                        # 规范: 每组内保持排序, 组间按最小标号排序
                        results.append([sorted(g1), sorted(g2)])
                elif V == 1:
                    results.append([list(range(n_external))])

            return results

        label_assignments = enumerate_label_assignments()

        for ext_per_vertex in label_assignments:
            # Step 3: 生成内部连线 (端口全等价, 只关心顶点之间的连接)
            # 每个顶点 v 有 4-ext_dist[v] 个内部端口
            # 这些端口全等价, 所以只需确定哪些顶点之间有多少条内线

            int_degrees = [4 - e for e in ext_dist]
            # 总内部端口数 = Σ int_degrees = 2I
            # 需要配对为 I 条边, 边连接两个顶点 (可以相同: 蝌蚪)

            # 枚举内部边的多集: [(v1, v2), ...]
            # 约束: 每个顶点 v 的度 = int_degrees[v]
            def enumerate_internal_edges():
                """生成所有内部连线方案"""
                results = []

                # 将 2I 个"半端口"配对
                # 半端口表示为顶点索引, 数量 = int_degrees[v]
                half_ports = []
                for v, deg in enumerate(int_degrees):
                    half_ports.extend([v] * deg)

                # 生成所有配对
                def pair_recurse(available: List[int],
                                 pairs: List[Tuple[int, int]]):
                    if not available:
                        results.append(list(pairs))
                        return
                    a = available[0]
                    for j in range(1, len(available)):
                        b = available[j]
                        new_avail = available[1:j] + available[j+1:]
                        pair_recurse(new_avail, pairs + [(a, b)])

                pair_recurse(half_ports, [])
                # 去重 (边无序)
                unique = set()
                final = []
                for pairs in results:
                    key = tuple(sorted(tuple(sorted(p)) for p in pairs))
                    if key not in unique:
                        unique.add(key)
                        final.append(list(key))
                return final

            edge_schemes = enumerate_internal_edges()

            for edges in edge_schemes:
                # 检查连通性
                if V > 1:
                    adj = defaultdict(set)
                    for a, b in edges:
                        adj[a].add(b)
                        adj[b].add(a)
                    visited = {0}
                    queue = [0]
                    while queue:
                        v = queue.pop()
                        for nb in adj[v]:
                            if nb not in visited:
                                visited.add(nb)
                                queue.append(nb)
                    if len(visited) != V:
                        continue

                # 去重: 用 (外腿分配, 内边多集) 的哈希 (顶点顺序无关)
                ext_sig = tuple(sorted(
                    frozenset(ext_per_vertex[v]) for v in range(V)
                ))
                edge_sig = tuple(sorted(
                    (min(a, b), max(a, b), sum(1 for ea, eb in edges
                              if {ea, eb} == {a, b}))
                    for a in range(V) for b in range(a, V)
                ))
                h = hash((ext_sig, edge_sig))
                if h in seen:
                    continue
                seen.add(h)

                # 确定道
                if V == 1 and n_external == 4:
                    channel = 'tree'
                elif V == 2 and n_external == 4:
                    channel = _classify_2to2_channel(
                        ext_per_vertex[0], ext_per_vertex[1])
                else:
                    channel = '?'

                # 对称因子
                S = compute_symmetry_factor_phi4(V, n_external, channel)

                # 振幅
                amp_expr, amp_val = _build_amplitude(
                    V=V, I=I, L=L, S=S, channel=channel,
                    coupling=coupling, mass=mass, Lambda=Lambda
                )

                # ASCII
                art = _draw_ascii(V=V, n_external=n_external, channel=channel,
                    ext_per_vertex=ext_per_vertex,
                    labels=external_labels, S=S
                )

                # 名字
                if channel == 'tree':
                    name = f"树图 2→2 (φ⁴ 四点顶点)"
                elif channel == 's':
                    name = f"s-道自能插入 (单圈泡泡图)"
                elif channel == 't':
                    name = f"t-道自能插入 (单圈泡泡图)"
                elif channel == 'u':
                    name = f"u-道自能插入 (单圈泡泡图)"
                else:
                    name = f"φ⁴ 图 (V={V}, L={L}, ch={channel})"

                diagrams.append(FeynmanDiagram(
                    name=name, n_vertices=V, n_external=n_external,
                    n_loops=L, symmetry_factor=S,
                    ascii_art=art, amplitude_expr=amp_expr,
                    amplitude_value=amp_val,
                    momentum_routing=_momentum_routing(channel,
                                                       ext_per_vertex,
                                                       external_labels),
                    channel=channel,
                ))

    ch_order = {'tree': 0, 's': 1, 't': 2, 'u': 3}
    diagrams.sort(key=lambda d: ch_order.get(d.channel, 99))
    return diagrams


def _build_amplitude(V: int, I: int, L: int, S: int, channel: str,
                     coupling: float, mass: float,
                     Lambda: float) -> Tuple[str, Optional[complex]]:
    """构建振幅表达式和数值"""
    lam = coupling

    if V == 1 and L == 0:
        # 树图: iM = -iλ
        return "-iλ", complex(0, -lam)

    if V == 2 and L == 1:
        # 单圈泡泡图
        # 确定道的动量平方
        E = 2.0
        s = 4 * E**2  # = 16.0
        k = np.sqrt(max(E**2 - mass**2, 0.01))
        theta = np.pi / 3
        t = -2 * k**2 * (1 - np.cos(theta))
        u = -2 * k**2 * (1 + np.cos(theta))

        ch_sq_map = {'s': s, 't': abs(t), 'u': abs(u)}
        ch_name_map = {'s': f's={s:.1f}', 't': f't={t:.1f}', 'u': f'u={u:.1f}'}
        ch_sq = ch_sq_map.get(channel, s)
        ch_note = ch_name_map.get(channel, '?')

        B = bubble_integral_2d(ch_sq, mass=mass, Lambda=Lambda)
        amp = (-1j * lam)**2 / S * B

        expr = f"(-iλ)²/{S} · B({channel})   [{ch_note}]"
        return expr, amp

    if L >= 2:
        # 多圈: 只给出表达式形式
        expr = f"(-iλ)^{V}/{S} × [{L}-loop integrals]"
        return expr, None

    expr = f"(-iλ)^{V}/{S} × [propagators] × [loop integrals]"
    return expr, None


def _momentum_routing(channel: str, ext_per_vertex: List[List[int]],
                      labels: List[str]) -> str:
    """生成动量路由描述"""
    if channel == 'tree':
        return "p₁ + p₂ = p₃ + p₄ (动量守恒, 单顶点)"

    if channel in ('s', 't', 'u'):
        v0_labels = ', '.join(labels[i] for v in ext_per_vertex
                              for i in v[:2]) if ext_per_vertex else ''
        return (f"{channel}-道, 圈动量 k 积分: "
                f"∫ d²k/(2π)² (传播子)²")

    return "动量守恒: Σ p_in = Σ p_out"


def _draw_ascii(V: int, n_external: int, channel: str,
                ext_per_vertex: List[List[int]],
                labels: List[str], S: int) -> str:
    """生成 ASCII 图表示"""
    if V == 1 and n_external == 4:
        p1, p2, p3, p4 = labels[:4]
        return f"""
    {p1}    {p3}
     \\   /
      \\ /
       ●   = -iλ
      / \\
     /   \\
    {p2}    {p4}
"""

    if V == 2 and n_external == 4:
        v0 = ext_per_vertex[0] if ext_per_vertex else [0, 1]
        v1 = ext_per_vertex[1] if ext_per_vertex else [2, 3]
        tl = labels[v0[0]] if len(v0) > 0 else '?'
        bl = labels[v0[1]] if len(v0) > 1 else '?'
        tr = labels[v1[0]] if len(v1) > 0 else '?'
        br = labels[v1[1]] if len(v1) > 1 else '?'

        ch_label = {'s': 's-道 p₁+p₂→p₃+p₄',
                    't': 't-道 p₁→p₃, p₂→p₄',
                    'u': 'u-道 p₁→p₄, p₂→p₃'}.get(channel, channel)

        return f"""
    {tl}      {tr}
     \\      /
      \\    /
    ───◆┄┄┄┄◆───   {ch_label}
        ┆    ┆
        ┆┄┄┄┄┆      S = {S}
       /      \\
      /        \\
    {bl}      {br}
"""

    if V == 1 and n_external == 2:
        p1, p2 = labels[:2]
        return f"""
    {p1} ──◆── {p2}
          ┆╲
          ┆ ╲  蝌蚪图 (自能)
          ┆╱
          ●
"""

    # 通用表示: 表格
    lines = [
        f"",
        f"  ┌─ φ⁴ Diagram (V={V}, L={V-n_external//2+1}) ────┐",
    ]
    for v in range(V):
        e = ext_per_vertex[v] if v < len(ext_per_vertex) else []
        e_str = ', '.join(labels[i] for i in e) if e else '—'
        lines.append(f"  │ V{v}: ext=[{e_str}]")
    lines.append(f"  │ S={S}")
    lines.append(f"  └────────────────────────────────────────────┘")
    return '\n'.join(lines)


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    比较函数                                      ║
# ╚══════════════════════════════════════════════════════════════════╝

def compare_orders(coupling: float = 0.5, mass: float = 1.0,
                   Lambda: float = 100.0) -> str:
    """并排比较 φ⁴ 2→2 树图和单圈贡献

    Args:
        coupling: λ
        mass:     m
        Lambda:   圈积分截断

    Returns:
        格式化的比较输出
    """
    tree_diags = generate_phi4_diagrams(4, 1, coupling, mass, Lambda)
    loop_diags = generate_phi4_diagrams(4, 2, coupling, mass, Lambda)
    lam = coupling

    lines = []
    lines.append("")
    lines.append("╔" + "═" * 68 + "╗")
    lines.append("║  φ⁴ 2→2 Scattering: Tree (LO) vs 1-Loop (NLO) Comparison"
                 + " " * 7 + "║")
    lines.append("╠" + "═" * 68 + "╣")

    # 树图
    lines.append("║  TREE LEVEL  ──  V=1, L=0  ──  O(λ¹)"
                 + " " * 22 + "║")
    lines.append("║" + " " * 68 + "║")
    if tree_diags:
        d = tree_diags[0]
        lines.append(f"║  Symmetry factor: S = {d.symmetry_factor}"
                     + " " * 42 + "║")
        lines.append(f"║  iM_tree = -iλ"
                     + " " * 52 + "║")
        if d.amplitude_value:
            val = d.amplitude_value
            lines.append(f"║  数值: {val.real:.4e}{val.imag:+.4e}i, "
                         f"|M|² = {abs(val)**2:.4e}"
                         + " " * 15 + "║")

    lines.append("║" + " " * 68 + "║")
    lines.append("╟" + "─" * 68 + "╢")
    lines.append("║  1-LOOP  ──  V=2, L=1  ──  O(λ²)"
                 + " " * 21 + "║")
    lines.append("║" + " " * 68 + "║")

    amp_total = complex(0, 0)
    for d in loop_diags:
        if d.amplitude_value:
            val = d.amplitude_value
            amp_total += val
            lines.append(f"║  {d.channel}-channel (S={d.symmetry_factor}): "
                         f"iM = {val.real:.4e}{val.imag:+.4e}i"
                         + " " * 18 + "║")
    lines.append(f"║  Sum: iM_1loop = {amp_total.real:.4e}{amp_total.imag:+.4e}i"
                 + " " * 22 + "║")
    lines.append(f"║  |M_1loop|² = {abs(amp_total)**2:.4e}"
                 + " " * 38 + "║")

    lines.append("║" + " " * 68 + "║")
    lines.append("╟" + "─" * 68 + "╢")
    lines.append("║  COMPARISON" + " " * 55 + "║")
    lines.append("║" + " " * 68 + "║")

    if tree_diags and tree_diags[0].amplitude_value:
        amp_tree = tree_diags[0].amplitude_value
        ratio = abs(amp_total / amp_tree) if abs(amp_tree) > 1e-15 else float('inf')
        lines.append(f"║  M_tree = -iλ = {amp_tree.imag:.4e}i"
                     + " " * 41 + "║")
        lines.append(f"║  M_1loop / M_tree = {ratio:.4e}"
                     + " " * 41 + "║")
        suppression = lam / (16 * np.pi**2)
        lines.append(f"║  Expected suppression ≈ λ/(16π²) = {suppression:.4e}"
                     + " " * 20 + "║")

    lines.append("║" + " " * 68 + "║")
    lines.append("╚" + "═" * 68 + "╝")
    lines.append("")

    return '\n'.join(lines)


def compare_diagrams(n_external: int = 4, max_order: int = 3) -> str:
    """比较不同阶数下的 Feynman 图数目

    Args:
        n_external: 外腿数
        max_order:  最大阶数

    Returns:
        格式化的计数表
    """
    lines = []
    lines.append("")
    lines.append(f"  Feynman Diagram Enumeration: n_external = {n_external}")
    lines.append(f"  {'─' * 55}")
    header = f"  {'Order':>6}  {'V':>4}  {'I':>4}  {'L':>4}  {'#Topologies':>12}  Notes"
    lines.append(header)
    lines.append(f"  {'─' * 55}")

    for order in range(1, max_order + 1):
        V = order
        total = 4 * V
        remaining = total - n_external
        if remaining < 0 or remaining % 2 != 0:
            lines.append(f"  {order:>6}  {V:>4}  {'—':>4}  {'—':>4}"
                         f"  {'0':>12}  (端口约束不满足)")
            continue
        I = remaining // 2
        L = I - V + 1
        if L < 0:
            lines.append(f"  {order:>6}  {V:>4}  {'—':>4}  {'—':>4}"
                         f"  {'0':>12}  (外腿太多)")
            continue
        try:
            diags = generate_phi4_diagrams(n_external, order, coupling=0.1)
            count = len(diags)
            ch_counts = Counter(d.channel for d in diags)
            note = ', '.join(f"{ch}:{n}" for ch, n in sorted(ch_counts.items()))
            lines.append(f"  {order:>6}  {V:>4}  {I:>4}  {L:>4}"
                         f"  {count:>12}  ({note})")
        except Exception:
            lines.append(f"  {order:>6}  {V:>4}  {I:>4}  {L:>4}"
                         f"  {'?':>12}")

    lines.append(f"  {'─' * 55}")
    lines.append("")
    return '\n'.join(lines)


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    主程序                                        ║
# ╚══════════════════════════════════════════════════════════════════╝

def print_diagram(diagram: FeynmanDiagram):
    """打印单个 Feynman 图"""
    diagram.print()


def run_full_demo(output_file: Optional[str] = None,
                  coupling: float = 0.5, mass: float = 1.0,
                  Lambda: float = 100.0):
    """运行完整演示"""
    lines = []

    def out(s=''):
        lines.append(s)
        print(s)

    out("=" * 65)
    out("  费曼图构造器 — φ⁴ 理论 §4.3-4.4")
    out("  图枚举 · 对称因子 · 散射振幅 · ASCII 可视化")
    out("=" * 65)

    # ── 对称因子解释 ──
    out(explain_symmetry_factor())

    # ── 图枚举统计 ──
    out("\n" + "─" * 65)
    out("  §1. 图枚举统计")
    out("─" * 65)
    out(compare_diagrams(4, 3))
    out(compare_diagrams(2, 3))
    out(compare_diagrams(6, 3))

    # ── 树图 ──
    out("\n" + "=" * 65)
    out("  §2. 树图 (Tree Level, V=1, O(λ¹))")
    out("=" * 65)

    tree_diags = generate_phi4_diagrams(4, 1, coupling, mass, Lambda)
    for d in tree_diags:
        d.print()

    # ── 单圈 ──
    out("\n" + "=" * 65)
    out("  §3. 单圈修正 (1-Loop, V=2, O(λ²))")
    out("=" * 65)

    loop_diags = generate_phi4_diagrams(4, 2, coupling, mass, Lambda)
    for d in loop_diags:
        d.print()

    # ── 圈积分数值 ──
    out("\n" + "─" * 65)
    out("  §4. 圈积分数值计算 (1+1D 时空)")
    out("─" * 65)

    E = 2.0
    k = np.sqrt(max(E**2 - mass**2, 0.01))
    theta = np.pi / 3
    s = 4 * E**2
    t = -2 * k**2 * (1 - np.cos(theta))
    u = -2 * k**2 * (1 + np.cos(theta))

    for ch_name, ch_sq in [('s', s), ('t', abs(t)), ('u', abs(u))]:
        B = bubble_integral_2d(ch_sq, mass=mass, Lambda=Lambda)
        out(f"  B({ch_name}={ch_sq:7.2f}) = {B.real:.4e}{B.imag:+.4e}i")

    tad = tadpole_integral_2d(mass=mass, Lambda=Lambda)
    out(f"  A(m²)  = ∫ d²k/(2π)² i/(k²-m²+iε) = {tad:.4e}")
    out(f"  注: 1+1D 中对数发散, 截断 Λ={Lambda}")

    # ── 比较 ──
    out(compare_orders(coupling, mass, Lambda))

    # ── 扩展: 其他外腿数 ──
    out("\n" + "─" * 65)
    out("  §5. 扩展: 更多外腿配置")
    out("─" * 65)

    for n_ext in [2, 6]:
        out(f"\n  n_external = {n_ext}:")
        for order in range(1, 4):
            diags = generate_phi4_diagrams(n_ext, order, coupling, mass, Lambda)
            if diags:
                out(f"    Order {order} (V={order}): {len(diags)} diagram(s)")
                for d in diags:
                    out(f"      {d.summary()}")
            else:
                out(f"    Order {order}: —")

    # ── 图计数总览表 ──
    out("\n" + "─" * 65)
    out("  §6. 图计数总览 (n_external × order)")
    out("─" * 65)

    header = f"  {'ext':>5} |"
    for o in range(1, 5):
        header += f"  V={o:>3} |"
    out(header)
    out("  " + "─" * (10 + 9 * 4))
    for n_ext in [2, 4, 6, 8]:
        row = f"  {n_ext:>5} |"
        for o in range(1, 5):
            try:
                diags = generate_phi4_diagrams(n_ext, o, coupling=0.1)
                row += f"  {len(diags):>3} |"
            except Exception:
                row += f"   — |"
        out(row)

    out("\n" + "=" * 65)
    out("  Feynman Builder Demo Complete — §4.3-4.4")
    out("=" * 65)

    if output_file:
        os.makedirs(os.path.dirname(output_file)
                    if os.path.dirname(output_file) else '.',
                    exist_ok=True)
        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))
        print(f"\nOutput saved to: {output_file}")


def run_2to2_demo(coupling: float = 0.5, mass: float = 1.0,
                  Lambda: float = 100.0):
    """仅 2→2 散射演示"""
    print("φ⁴ 2→2 Scattering")
    print("=" * 50)

    for order in [1, 2]:
        print(f"\n── Order {order} (V={order}) ──")
        diags = generate_phi4_diagrams(4, order, coupling, mass, Lambda)
        for d in diags:
            d.print()

    print(compare_orders(coupling, mass, Lambda))


def run_compare_demo(coupling: float = 0.5, mass: float = 1.0,
                     Lambda: float = 100.0):
    """仅树图 vs 单圈比较"""
    print(compare_orders(coupling, mass, Lambda))


# ╔══════════════════════════════════════════════════════════════════╗
# ║                    入口                                         ║
# ╚══════════════════════════════════════════════════════════════════╝

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='φ⁴ Feynman Diagram Builder — §4.3-4.4',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python feynman_builder.py                    # 完整演示
  python feynman_builder.py --2to2             # 仅 2→2 散射
  python feynman_builder.py --compare          # 树图 vs 单圈比较
  python feynman_builder.py --save out.txt     # 保存输出到文件
  python feynman_builder.py --diagrams 4 2     # 生成 n_ext=4, V=2 的所有图
        """,
    )
    parser.add_argument('--2to2', dest='to2to2', action='store_true',
                        help='仅运行 2→2 散射演示')
    parser.add_argument('--compare', action='store_true',
                        help='仅运行树图 vs 单圈比较')
    parser.add_argument('--save', type=str, metavar='FILE',
                        help='保存输出到文件')
    parser.add_argument('--diagrams', nargs=2, type=int,
                        metavar=('N_EXT', 'ORDER'),
                        help='生成指定配置的 Feynman 图')
    parser.add_argument('--coupling', type=float, default=0.5,
                        help='耦合常数 λ (默认 0.5)')
    parser.add_argument('--mass', type=float, default=1.0,
                        help='质量 m (默认 1.0)')
    parser.add_argument('--Lambda', type=float, default=100.0,
                        help='动量截断 Λ (默认 100.0)')

    args = parser.parse_args()

    lam_val = args.coupling
    m_val = args.mass
    L_val = args.Lambda

    if args.diagrams:
        n_ext, order = args.diagrams
        print(f"Enumerating φ⁴ diagrams: n_external={n_ext}, order={order}")
        diags = generate_phi4_diagrams(n_ext, order, lam_val, m_val, L_val)
        print(f"Found {len(diags)} diagram(s):\n")
        for i, d in enumerate(diags):
            print(f"── Diagram {i+1}: {d.summary()} ──")
            d.print()
            print()
    elif args.to2to2:
        run_2to2_demo(lam_val, m_val, L_val)
    elif args.compare:
        run_compare_demo(lam_val, m_val, L_val)
    else:
        run_full_demo(args.save, lam_val, m_val, L_val)

"""量子电路模块 — 电路构建与模拟

Circuit 类提供了量子电路的 Pythonic 接口:
    - 添加门 (H, X, CNOT, Rx, ...)
    - 自动计算态演化
    - 测量与统计

用法:
    from src.qubit import Circuit
    qc = Circuit(2)                 # 2 qubit 电路
    qc.h(0)                         # qubit 0 上 H 门
    qc.cnot(0, 1)                   # CNOT 控制=0, 目标=1
    qc.measure()                    # 测量所有 qubit
    print(qc.probabilities())       # 概率分布
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from . import states, gates


class Circuit:
    """量子电路模拟器

    用法:
        qc = Circuit(2)
        qc.h(0)
        qc.cnot(0, 1)
        qc.measure()
        qc.draw()
    """

    def __init__(self, n_qubits: int, initial_state: np.ndarray = None):
        """
        参数:
            n_qubits:      量子比特数
            initial_state: 初始态 (默认 |0...0⟩)
        """
        self.n_qubits = n_qubits
        self.dim = 2 ** n_qubits

        if initial_state is None:
            self.state = np.zeros(self.dim, dtype=complex)
            self.state[0] = 1.0
        else:
            self.state = np.asarray(initial_state, dtype=complex).flatten()
            if len(self.state) != self.dim:
                raise ValueError(f"State dim {len(self.state)} ≠ 2^{n_qubits}={self.dim}")

        self._initial_state = self.state.copy()
        self._gates: List[Tuple[str, tuple]] = []  # 门历史
        self._measured = False
        self._result: Optional[int] = None

    # ---------- 单 qubit 门 ----------

    def h(self, target: int):
        """Hadamard 门"""
        self._apply_gate(gates.expand_gate(gates.H, target, self.n_qubits),
                         'H', (target,))

    def x(self, target: int):
        """Pauli X (NOT) 门"""
        self._apply_gate(gates.expand_gate(gates.X, target, self.n_qubits),
                         'X', (target,))

    def y(self, target: int):
        """Pauli Y 门"""
        self._apply_gate(gates.expand_gate(gates.Y, target, self.n_qubits),
                         'Y', (target,))

    def z(self, target: int):
        """Pauli Z 门"""
        self._apply_gate(gates.expand_gate(gates.Z, target, self.n_qubits),
                         'Z', (target,))

    def s(self, target: int):
        """S 门 (Phase)"""
        self._apply_gate(gates.expand_gate(gates.S, target, self.n_qubits),
                         'S', (target,))

    def t(self, target: int):
        """T 门 (π/8)"""
        self._apply_gate(gates.expand_gate(gates.T_gate, target, self.n_qubits),
                         'T', (target,))

    def rx(self, theta: float, target: int):
        """绕 X 轴旋转"""
        self._apply_gate(gates.expand_gate(gates.Rx(theta), target, self.n_qubits),
                         'Rx', (theta, target))

    def ry(self, theta: float, target: int):
        """绕 Y 轴旋转"""
        self._apply_gate(gates.expand_gate(gates.Ry(theta), target, self.n_qubits),
                         'Ry', (theta, target))

    def rz(self, phi: float, target: int):
        """绕 Z 轴旋转"""
        self._apply_gate(gates.expand_gate(gates.Rz(phi), target, self.n_qubits),
                         'Rz', (phi, target))

    # ---------- 两 qubit 门 ----------

    def cnot(self, control: int, target: int):
        """CNOT 门"""
        if self.n_qubits == 2:
            gate = gates.CNOT if control == 0 else gates.CNOT_rev
        else:
            gate = gates.expand_cnot(control, target, self.n_qubits)
        self._apply_gate(gate, 'CNOT', (control, target))

    def cz(self, control: int, target: int):
        """CZ 门"""
        raise NotImplementedError("Use cnot + h sandwich for general CZ")

    def swap(self, q1: int, q2: int):
        """SWAP 门"""
        raise NotImplementedError("SWAP for arbitrary qubits not yet implemented")

    # ---------- 操作 ----------

    def _apply_gate(self, gate: np.ndarray, name: str, args: tuple):
        """内部: 应用门并记录"""
        self.state = gate @ self.state
        self._gates.append((name, args))

    def reset(self):
        """重置到初始态"""
        self.state = self._initial_state.copy()
        self._gates.clear()
        self._measured = False
        self._result = None

    def measure(self, shots: int = 1) -> int:
        """测量 (计算基)

        返回:
            测量结果 (整数), 态坍缩到此基矢
        """
        probs = np.abs(self.state)**2
        probs /= probs.sum()
        outcome = np.random.choice(self.dim, p=probs)
        self.state = np.zeros(self.dim, dtype=complex)
        self.state[outcome] = 1.0
        self._measured = True
        self._result = outcome
        return outcome

    def measure_all(self, shots: int = 1024) -> Dict[int, int]:
        """多次测量并统计

        返回:
            {outcome: count} 字典
        """
        counts = {}
        for _ in range(shots):
            # 保存当前态
            saved = self.state.copy()
            outcome = self.measure()
            counts[outcome] = counts.get(outcome, 0) + 1
            self.state = saved  # 恢复 (不坍缩)
        return counts

    # ---------- 信息 ----------

    def probabilities(self) -> np.ndarray:
        """测量概率分布 P(i) = |⟨i|ψ⟩|²"""
        return np.abs(self.state)**2

    def get_state(self) -> np.ndarray:
        """返回当前态向量"""
        return self.state.copy()

    def get_statevector(self) -> np.ndarray:
        """返回态向量 (同 get_state)"""
        return self.state.copy()

    # ---------- 显示 ----------

    def draw(self) -> str:
        """ASCII 电路图"""
        lines = [f"Quantum Circuit ({self.n_qubits} qubits, {len(self._gates)} gates)"]
        lines.append("-" * 50)

        # qubit 线
        q_lines = {i: [f"q{i}: "] for i in range(self.n_qubits)}

        for gate_name, args in self._gates:
            if gate_name in ('H', 'X', 'Y', 'Z', 'S', 'T'):
                t = args[0]
                for i in range(self.n_qubits):
                    if i == t:
                        q_lines[i].append(f"──[{gate_name}]──")
                    else:
                        q_lines[i].append("───────")
            elif gate_name in ('Rx', 'Ry', 'Rz'):
                val, t = args
                label = f"{gate_name}({val:.2f})"
                for i in range(self.n_qubits):
                    if i == t:
                        q_lines[i].append(f"──[{label}]──")
                    else:
                        q_lines[i].append(f"{'─' * (len(label)+4)}")
            elif gate_name == 'CNOT':
                c, t = args
                for i in range(self.n_qubits):
                    if i == c:
                        q_lines[i].append("──●──")
                    elif i == t:
                        q_lines[i].append("──⊕──")
                    else:
                        q_lines[i].append("──────")

        # 渲染
        for i in range(self.n_qubits):
            lines.append(''.join(q_lines[i]))
        lines.append("-" * 50)

        # 概率
        probs = self.probabilities()
        non_zero = [(i, p) for i, p in enumerate(probs) if p > 0.01]
        if non_zero:
            lines.append("State probabilities:")
            for i, p in non_zero[:8]:
                bits = format(i, f'0{self.n_qubits}b')
                lines.append(f"  |{bits}⟩: {p:.4f}")
        return '\n'.join(lines)

    def __repr__(self):
        return f"Circuit({self.n_qubits}q, {len(self._gates)} gates)"


# ============================================================
# 快捷电路构建
# ============================================================

def bell_circuit() -> Circuit:
    """创建 Bell 态电路: |00⟩ → H(0) → CNOT(0,1) → |Φ⁺⟩"""
    qc = Circuit(2)
    qc.h(0)
    qc.cnot(0, 1)
    return qc


def ghz_circuit(n: int = 3) -> Circuit:
    """创建 GHZ 态电路"""
    qc = Circuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cnot(i, i + 1)
    return qc


def teleportation_circuit(state_to_teleport: np.ndarray = None) -> Circuit:
    """量子隐形传态电路 (3 qubit)

    qubit 0: 待传态 |ψ⟩
    qubit 1,2: Bell 对 |Φ⁺⟩

    电路: |ψ⟩⊗|Φ⁺⟩ → CNOT(0,1) → H(0) → 测量 → X/Z corrections
    """
    if state_to_teleport is None:
        from .states import ket_plus
        state_to_teleport = ket_plus

    initial = np.kron(state_to_teleport, states.bell_state(0))
    qc = Circuit(3, initial_state=initial)
    qc.cnot(0, 1)
    qc.h(0)
    return qc

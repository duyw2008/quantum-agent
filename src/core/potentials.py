"""势函数模块

支持的势函数:
    - InfiniteWell:      无限深势阱 (平滑边界)
    - Harmonic:          谐振子 V(x) = ½ m ω² x²
    - PotentialBarrier:  矩形势垒 V(x) = V₀ (|x| < w/2), 0 elsewhere
    - FiniteWell:        有限深势阱
    - DoubleWell:        双势阱 (quartic)
    - Morse:             Morse 势 V(x) = Dₑ [1 - exp(-α(x-x₀))]²
    - Coulomb1D:         一维软核库仑势 V(x) = -Z / sqrt(x² + a²)
    - Periodic:          周期势 V(x) = A cos(kx)
    - Custom:            自定义势函数
    - Step:              阶梯势
"""

import numpy as np
from typing import Callable, Optional, Dict, Any
from abc import ABC, abstractmethod
from enum import Enum


class PotentialType(Enum):
    INFINITE_WELL = "infinite_well"
    HARMONIC = "harmonic"
    POTENTIAL_BARRIER = "potential_barrier"
    FINITE_WELL = "finite_well"
    DOUBLE_WELL = "double_well"
    MORSE = "morse"
    COULOMB_1D = "coulomb_1d"
    PERIODIC = "periodic"
    STEP = "step"
    ZERO = "zero"
    CUSTOM = "custom"


# ============================================================
# 抽象基类
# ============================================================

class Potential(ABC):
    """势函数抽象基类 V(x)"""

    def __init__(self, name: str = "Potential"):
        self.name = name

    @abstractmethod
    def __call__(self, x: np.ndarray) -> np.ndarray:
        """V(x): 势函数调用接口"""
        ...

    @property
    def label(self) -> str:
        """用于图例或显示的字符串标签"""
        return self.name

    @staticmethod
    def _x_range_info(x: np.ndarray) -> Dict[str, float]:
        """提取网格的基本信息"""
        return {
            'dx': x[1] - x[0],
            'n_points': len(x),
            'x_min': x[0],
            'x_max': x[-1],
        }


# ============================================================
# 具体势函数
# ============================================================

class InfiniteWell(Potential):
    """无限深势阱 (平滑边界)

    V(x) = 0  for |x| < a/2
    V(x) → ∞ for |x| ≥ a/2  (用大值近似)

    参数:
        width: 阱宽 a
        wall_height: 壁高度 (默认 1e6)
        steepness: 边界平滑度 (越大越陡)
    """

    def __init__(self, width: float = 2.0, wall_height: float = 1e6,
                 steepness: float = None):
        super().__init__(f"Infinite Well (a={width})")
        self.width = width
        self.wall_height = wall_height
        self.steepness = steepness if steepness is not None else 50.0 / width

    def __call__(self, x: np.ndarray) -> np.ndarray:
        half = self.width / 2
        eps = self.width * 0.001  # 微小偏移避免数值溢出
        V = np.zeros_like(x)
        # 左壁: 在 x < -half 时上升，sigmoid: 1/(1+exp(+s*(x+half)))
        V += self.wall_height / (1 + np.exp(self.steepness * (x + half + eps)))
        # 右壁: 在 x > +half 时上升，sigmoid: 1/(1+exp(-s*(x-half)))
        V += self.wall_height / (1 + np.exp(-self.steepness * (x - half - eps)))
        return V


class Harmonic(Potential):
    """谐振子势

    V(x) = ½ m ω² x²

    解析解能量: E_n = ℏω(n + ½),  n = 0,1,2,...
    特征长度: a_ho = sqrt(ℏ / mω)
    """

    def __init__(self, omega: float = 1.0, mass: float = 1.0):
        super().__init__(f"Harmonic (ω={omega})")
        self.omega = omega
        self.mass = mass

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return 0.5 * self.mass * self.omega**2 * x**2

    @property
    def characteristic_length(self) -> float:
        """特征长度 a_ho = √(ℏ/mω)"""
        return np.sqrt(1.0 / (self.mass * self.omega))


class PotentialBarrier(Potential):
    """矩形势垒 / 势阱

    V(x) = V0 (|x| < w/2), 0 elsewhere
    V0 > 0: 势垒 (barrier)
    V0 < 0: 势阱 (well)
    """

    def __init__(self, height: float = 10.0, width: float = 1.0):
        label = f"{'Barrier' if height > 0 else 'Well'} (V₀={height}, w={width})"
        super().__init__(label)
        self.height = height
        self.width = width

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return np.where(np.abs(x) < self.width / 2, self.height, 0.0)


class FiniteWell(Potential):
    """有限深势阱 (负值表示阱)

    V(x) = -|depth| (|x| < w/2), 0 elsewhere
    """

    def __init__(self, depth: float = 10.0, width: float = 1.0):
        super().__init__(f"Finite Well (depth={depth}, w={width})")
        self.depth = abs(depth)
        self.width = width

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return np.where(np.abs(x) < self.width / 2, -self.depth, 0.0)


class DoubleWell(Potential):
    """双势阱 (quartic 形式)

    V(x) = V₀ * [(x/a)² - 1]² - V₀
    极小值在 x = ±a, 阱底为 V = -V₀
    中心势垒高度为 V(0) = 0
    参数 separation = 2a (两阱之间的距离)
    """

    def __init__(self, separation: float = 2.0, depth: float = 10.0):
        super().__init__(f"Double Well (sep={separation}, depth={depth})")
        self.separation = separation
        self.depth = depth
        self.a = separation / 2  # 每个阱的位置

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.depth * ((x / self.a)**2 - 1)**2 - self.depth


class Morse(Potential):
    """Morse 势 (双原子分子)

    V(x) = Dₑ [1 - exp(-α(x - x₀))]²

    解析解能级: E_n = ℏω₀(n + ½) - [ℏω₀(n + ½)]² / (4Dₑ)
    其中 ω₀ = α √(2Dₑ/m)
    """

    def __init__(self, depth: float = 10.0, alpha: float = 1.0,
                 x0: float = 0.0, mass: float = 1.0):
        super().__init__(f"Morse (Dₑ={depth}, α={alpha})")
        self.depth = depth
        self.alpha = alpha
        self.x0 = x0
        self.mass = mass

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.depth * (1 - np.exp(-self.alpha * (x - self.x0)))**2

    @property
    def harmonic_omega(self) -> float:
        """等效谐振子频率 ω₀ = α √(2Dₑ/m)"""
        return self.alpha * np.sqrt(2 * self.depth / self.mass)

    def analytic_energy(self, n: int) -> float:
        """解析解能级 (n = 0, 1, 2, ...)"""
        omega0 = self.harmonic_omega
        e_harmonic = omega0 * (n + 0.5)
        return e_harmonic - e_harmonic**2 / (4 * self.depth)


class Coulomb1D(Potential):
    """一维软核库仑势

    V(x) = -Z / √(x² + a²)
    a > 0 是软化参数，避免在 x=0 处的奇点
    """

    def __init__(self, Z: float = 1.0, softening: float = 1.0):
        super().__init__(f"Coulomb 1D (Z={Z}, a={softening})")
        self.Z = Z
        self.softening = softening

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return -self.Z / np.sqrt(x**2 + self.softening**2)


class Periodic(Potential):
    """周期势 (Kronig-Penney 模型基础)

    V(x) = A cos(kx)
    用于模拟晶体中的电子运动
    """

    def __init__(self, amplitude: float = 1.0, k: float = 2 * np.pi):
        super().__init__(f"Periodic (A={amplitude}, k={k})")
        self.amplitude = amplitude
        self.k = k

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.amplitude * np.cos(self.k * x)


class StepPotential(Potential):
    """阶梯势

    V(x) = V0 * Θ(x - x0)    (单位阶跃)
    用于散射问题
    """

    def __init__(self, height: float = 5.0, x0: float = 0.0):
        super().__init__(f"Step (V₀={height}, x₀={x0})")
        self.height = height
        self.x0 = x0

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return np.where(x > self.x0, self.height, 0.0)


class ZeroPotential(Potential):
    """自由粒子: V(x) = 0"""

    def __init__(self):
        super().__init__("Zero (free particle)")

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)


class CustomPotential(Potential):
    """自定义势函数

    示例:
        V = CustomPotential(lambda x: 0.5 * x**4, name="Quartic")
    """

    def __init__(self, func: Callable[[np.ndarray], np.ndarray],
                 name: str = "Custom"):
        super().__init__(name)
        self._func = func

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self._func(x)


# ============================================================
# 工厂函数
# ============================================================

def create_potential(ptype: PotentialType, **kwargs) -> Potential:
    """工厂函数：根据类型创建势函数

    用法:
        V = create_potential(PotentialType.HARMONIC, omega=2.0, mass=1.0)
        V = create_potential(PotentialType.INFINITE_WELL, width=3.0)
        V = create_potential(PotentialType.CUSTOM,
                             func=lambda x: 0.5*x**4, name="Quartic")
    """
    ptype_map = {
        PotentialType.INFINITE_WELL: InfiniteWell,
        PotentialType.HARMONIC: Harmonic,
        PotentialType.POTENTIAL_BARRIER: PotentialBarrier,
        PotentialType.FINITE_WELL: FiniteWell,
        PotentialType.DOUBLE_WELL: DoubleWell,
        PotentialType.MORSE: Morse,
        PotentialType.COULOMB_1D: Coulomb1D,
        PotentialType.PERIODIC: Periodic,
        PotentialType.STEP: StepPotential,
        PotentialType.ZERO: ZeroPotential,
        PotentialType.CUSTOM: CustomPotential,
    }
    cls = ptype_map[ptype]
    return cls(**kwargs)

#!/usr/bin/env python3
"""
Quantum Agent — 量子力学智能体

交互式量子力学计算和可视化平台。

用法:
    python agent.py                    # 交互模式
    python agent.py --demo well        # 运行指定 demo
    python agent.py --demo all         # 运行所有 demo
    python agent.py --list             # 列出所有可用 demo

命令 (交互模式):
    evolve <potential> [params]  — 演化波函数
    matrix                       — 矩阵力学计算
    eigenstates <potential>      — 计算本征态
    plot <type>                  — 绘图
    animate <type>               — 动画
    demo <name>                  — 运行 demo
    help                         — 帮助
    quit / exit                  — 退出
"""

import sys
import os
import numpy as np
import argparse

# readline: 命令行历史和上下键选择
try:
    import readline
    HAS_READLINE = True
except ImportError:
    HAS_READLINE = False

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    yaml = None

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.core import (
    Grid, WaveFunction, create_potential, PotentialType,
    create_solver,
)
from src.matrix import MatrixMechanics
from src.viz import (
    plot_potential, plot_wavefunction, plot_eigenstates,
    plot_energy_levels, plot_matrix_element, plot_phase_space,
    animate_evolution,
)


def load_config(path: str = None) -> dict:
    """加载配置文件"""
    if path is None:
        path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(path) and HAS_YAML:
        with open(path, 'r') as f:
            return yaml.safe_load(f)
    return {}


class QuantumAgent:
    """量子力学智能体"""

    def __init__(self, config_path: str = None):
        self.config = load_config(config_path)

        # 默认参数
        self.hbar = self.config.get('hbar', 1.0)
        self.mass = self.config.get('mass', 1.0)

        numerics = self.config.get('numerics', {})
        self.default_x_range = numerics.get('x_range', [-10.0, 10.0])
        self.default_n_points = numerics.get('n_points', 1024)
        self.default_dt = numerics.get('dt', 0.001)
        self.default_t_max = numerics.get('t_max', 5.0)

        self.output_dir = os.path.join(os.path.dirname(__file__), 'output')
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, 'animations'), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, 'figures'), exist_ok=True)

        self.history = []

        # readline 历史文件
        self._history_file = os.path.join(os.path.expanduser('~'),
                                           '.quantum_agent_history')
        self._setup_readline()

        # 当前状态 (可被命令修改)
        self.current_grid = None
        self.current_potential = None
        self.current_wf = None
        self.current_result = None
        self.current_matrix = None

        print(self._welcome())

    def _welcome(self) -> str:
        msg = self.config.get('agent', {}).get('welcome_message',
                                                 'Welcome to Quantum Agent')
        return f"""
╔══════════════════════════════════════════════════════╗
║        {msg}         ║
║        v0.1.0 — 一维量子力学计算与可视化              ║
╚══════════════════════════════════════════════════════╝
Type 'help' for commands, 'demo all' to see examples.
"""

    # ============================================================
    # 命令处理
    # ============================================================

    def _parse_potential(self, pot_str: str):
        """解析势函数字符串

        示例:
            harmonic,omega=2.0
            infinite_well,width=3.0
            barrier,height=5.0,width=1.0
            double_well,separation=3.0,depth=8.0
        """
        parts = pot_str.split(',')
        ptype_str = parts[0].strip().lower()
        params = {}

        # 类型映射
        type_map = {
            'harmonic': PotentialType.HARMONIC,
            'ho': PotentialType.HARMONIC,
            'well': PotentialType.INFINITE_WELL,
            'infinite_well': PotentialType.INFINITE_WELL,
            'barrier': PotentialType.POTENTIAL_BARRIER,
            'potential_barrier': PotentialType.POTENTIAL_BARRIER,
            'finite_well': PotentialType.FINITE_WELL,
            'double_well': PotentialType.DOUBLE_WELL,
            'dw': PotentialType.DOUBLE_WELL,
            'morse': PotentialType.MORSE,
            'coulomb': PotentialType.COUlOMB_1D,
            'periodic': PotentialType.PERIODIC,
            'step': PotentialType.STEP,
            'zero': PotentialType.ZERO,
            'free': PotentialType.ZERO,
        }

        if ptype_str not in type_map:
            raise ValueError(f"未知势函数类型: {ptype_str}。支持: {list(type_map.keys())}")

        ptype = type_map[ptype_str]

        # 解析参数
        for part in parts[1:]:
            if '=' in part:
                key, val = part.split('=')
                key = key.strip()
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val.strip()

        return ptype, params

    def cmd_evolve(self, args: list):
        """演化波函数: evolve <potential> [options]

        选项:
            x0=0.0    初始位置
            p0=5.0    初始动量
            sigma=0.5 波包宽度
            dt=0.001  时间步长
            t_max=5.0 总时间
            n=1024    网格点数
            x_min=-10 空间下限
            x_max=10  空间上限
            method=ssfm  数值方法 (ssfm/cn)
            save=1    是否保存动画
        """
        if not args:
            print("用法: evolve <potential> [x0=X] [p0=X] [sigma=X] [dt=X] [t_max=X] [method=ssfm|cn]")
            print("示例: evolve harmonic,omega=2.0 x0=1.0 p0=3.0 t_max=10")
            return

        pot_str = args[0]
        kwargs = {'x0': 0.0, 'p0': 5.0, 'sigma': 0.5,
                  'dt': self.default_dt, 't_max': self.default_t_max,
                  'n': self.default_n_points, 'x_min': self.default_x_range[0],
                  'x_max': self.default_x_range[1], 'method': 'ssfm',
                  'save': True}

        for arg in args[1:]:
            if '=' in arg:
                k, v = arg.split('=', 1)
                try:
                    kwargs[k] = float(v) if '.' in v or v.lstrip('-').isdigit() else v
                except ValueError:
                    kwargs[k] = v

        ptype, params = self._parse_potential(pot_str)
        V = create_potential(ptype, **params)

        grid = Grid(kwargs['x_min'], kwargs['x_max'], int(kwargs['n']))
        wf = WaveFunction(grid)
        wf.set_gaussian(x0=kwargs['x0'], p0=kwargs['p0'], sigma=kwargs['sigma'])

        print(f"\n⚛  Initializing evolution...")
        print(f"   Potential: {V.label}")
        print(f"   Grid: {grid}")
        print(f"   Wavepacket: x₀={kwargs['x0']}, p₀={kwargs['p0']}, σ={kwargs['sigma']}")
        print(f"   Method: {kwargs['method']}, dt={kwargs['dt']}, t_max={kwargs['t_max']}")
        print(f"   Initial ⟨x⟩={wf.expectation_x():.4f}, ⟨p⟩={wf.expectation_p():.4f}")

        solver = create_solver(kwargs['method'], grid, V, self.hbar, self.mass)
        result = solver.evolve(wf, kwargs['t_max'], kwargs['dt'],
                               snapshot_interval=max(1, int(kwargs['t_max'] / kwargs['dt'] / 200)))

        # 输出统计
        print(f"\n✓ Evolution complete: {len(result.times)} snapshots")
        print(f"  Final ⟨x⟩ = {result.expectation_x[-1]:.4f}")
        print(f"  Final ⟨p⟩ = {result.expectation_p[-1]:.4f}")
        print(f"  Energy conservation: σ_E/Ē = {result.energy.std()/result.energy.mean():.2e}")
        print(f"  Norm conservation: max|1-norm| = {abs(result.norm_history - 1).max():.2e}")

        # 保存动画
        if kwargs.get('save', True):
            save_path = os.path.join(
                self.output_dir, 'animations',
                f"{ptype.value}_{kwargs['method']}.mp4"
            )
            print(f"  Generating animation → {save_path} ...")
            animate_evolution(result, V, save_path=save_path, fps=30, dark=True)
            print(f"  ✓ Saved to {save_path}")

        self.current_grid = grid
        self.current_potential = V
        self.current_wf = wf
        self.current_result = result
        self.history.append(('evolve', pot_str))

    def cmd_matrix(self, args: list):
        """矩阵力学计算: matrix [command]

        命令:
            report      — 打印算符关系报告
            x           — 坐标算符矩阵
            p           — 动量算符矩阵
            a / a_dag   — 产生/湮灭算符
            H [omega=N] — 哈密顿量
            eigen [k=N] — 本征值/本征态
            comm A B    — 对易子 [A, B]
        """
        if not self.current_matrix:
            self.current_matrix = MatrixMechanics(
                n_basis=50, hbar=self.hbar, mass=self.mass, omega=1.0
            )

        mm = self.current_matrix

        if not args or args[0] == 'report':
            print(mm.report())
            return

        cmd = args[0]

        if cmd == 'x':
            print("Coordinate operator x̂ (first 8×8):")
            print(mm.x[:8, :8].real)

        elif cmd == 'p':
            print("Momentum operator p̂ (first 8×8):")
            print(mm.p[:8, :8].real)

        elif cmd == 'a':
            print("Annihilation operator â (first 8×8):")
            print(mm.a[:8, :8])

        elif cmd == 'a_dag':
            print("Creation operator â† (first 8×8):")
            print(mm.a_dag[:8, :8])

        elif cmd == 'H':
            omega = float(args[1]) if len(args) > 1 else 1.0
            eigs, states = mm.eigensolve(k=10)
            print(f"Harmonic oscillator energies (ω={omega}):")
            for n, E in enumerate(eigs[:10]):
                expected = omega * (n + 0.5)
                print(f"  E[{n}] = {E:.6f}  (expected: {expected:.6f})")

        elif cmd == 'eigen':
            k = int(args[1]) if len(args) > 1 else 5
            eigs, states = mm.eigensolve(k=k)
            print(f"First {k} eigenvalues:")
            for i, E in enumerate(eigs):
                print(f"  E[{i}] = {E:.6f}")

        elif cmd == 'comm':
            if len(args) < 3:
                print("用法: matrix comm <A> <B>  (e.g., matrix comm x p)")
                return
            op_map = {'x': mm.x, 'p': mm.p, 'a': mm.a, 'a_dag': mm.a_dag}
            A = op_map.get(args[1])
            B = op_map.get(args[2])
            if A is None or B is None:
                print(f"未知算符。可用: {list(op_map.keys())}")
                return
            C = mm.commutator(A, B)
            result = mm.check_commutation(A, B)
            print(f"[{args[1]}, {args[2]}] Frobenius norm: {result['frobenius_norm']:.6e}")
            if args[2] == 'p' and args[1] == 'x':
                print(f"Expected: iℏI (should be near {self.hbar} * √N = {self.hbar*np.sqrt(mm._N):.4f})")

        else:
            print(f"未知矩阵命令: {cmd}")
            print("可用: report, x, p, a, a_dag, H, eigen, comm")

    def cmd_eigenstates(self, args: list):
        """计算本征态: eigenstates <potential> [n_states=5]"""
        if not args:
            print("用法: eigenstates <potential> [n_states=5]")
            print("示例: eigenstates harmonic,omega=2.0 5")
            return

        pot_str = args[0]
        n_states = int(args[1]) if len(args) > 1 else 5

        ptype, params = self._parse_potential(pot_str)
        V = create_potential(ptype, **params)

        grid = Grid(self.default_x_range[0], self.default_x_range[1],
                    self.default_n_points)

        save_path = os.path.join(self.output_dir, 'figures',
                                 f"eigenstates_{ptype.value}.png")
        fig, eigvals, eigvecs = plot_eigenstates(grid, V, n_states,
                                                  save_path=save_path)

        print(f"\n✓ Eigenstates computed: {n_states} states")
        print(f"  Energies:")
        for n, E in enumerate(eigvals):
            print(f"    E[{n}] = {E:.6f}")
        print(f"  Saved to {save_path}")

        self.current_grid = grid
        self.current_potential = V

    def cmd_plot(self, args: list):
        """绘图: plot <type> [options]

        类型: potential | wf | energy
        """
        if not args:
            print("用法: plot <potential|wf|energy>")
            return

        cmd = args[0]
        if cmd == 'potential' and self.current_potential:
            save_path = os.path.join(self.output_dir, 'figures', 'potential.png')
            plot_potential(self.current_potential, save_path=save_path)
            print(f"✓ Saved to {save_path}")

        elif cmd == 'wf' and self.current_wf:
            save_path = os.path.join(self.output_dir, 'figures', 'wavefunction.png')
            plot_wavefunction(self.current_wf, V=self.current_potential,
                              save_path=save_path)
            print(f"✓ Saved to {save_path}")

        else:
            print("请先运行 evolve 或 eigenstates。")

    def cmd_animate(self, args: list):
        """动画: animate <type>"""
        if not args:
            print("用法: animate <evolution>")
            return

        if args[0] == 'evolution' and self.current_result and self.current_potential:
            save_path = os.path.join(self.output_dir, 'animations', 'manual_animate.mp4')
            animate_evolution(self.current_result, self.current_potential,
                              save_path=save_path)
            print(f"✓ Saved to {save_path}")
        else:
            print("请先运行 evolve。")

    def cmd_demo(self, args: list):
        """运行 demo: demo <name|all>"""
        demos = {
            'harmonic': 'demos/harmonic_oscillator.py',
            'well': 'demos/infinite_well.py',
            'barrier': 'demos/potential_barrier.py',
            'matrix': 'demos/matrix_mechanics.py',
            'double_well': 'demos/double_well.py',
        }

        if not args or args[0] == 'all':
            for name, path in demos.items():
                print(f"\n{'='*50}")
                print(f"  Running demo: {name}")
                print(f"{'='*50}")
                self._run_demo(path)
            return

        name = args[0]
        if name in demos:
            self._run_demo(demos[name])
        else:
            print(f"未知 demo: {name}")
            print(f"可用: {list(demos.keys())}")

    def _run_demo(self, path: str):
        """执行 demo 脚本"""
        full_path = os.path.join(os.path.dirname(__file__), path)
        if os.path.exists(full_path):
            exec(open(full_path).read(), {
                '__name__': '__main__',
                'np': np,
                'Grid': Grid,
                'WaveFunction': WaveFunction,
                'create_potential': create_potential,
                'PotentialType': PotentialType,
                'create_solver': create_solver,
                'MatrixMechanics': MatrixMechanics,
                'animate_evolution': animate_evolution,
                'plot_potential': plot_potential,
                'plot_wavefunction': plot_wavefunction,
                'plot_eigenstates': plot_eigenstates,
                'output_dir': self.output_dir,
            })
        else:
            print(f"Demo 文件不存在: {full_path}")

    # ============================================================
    # readline 设置
    # ============================================================

    def _setup_readline(self):
        """设置 readline 历史"""
        if not HAS_READLINE:
            return
        # 设置历史文件
        readline.set_history_length(1000)
        # 加载历史
        try:
            readline.read_history_file(self._history_file)
        except (FileNotFoundError, PermissionError):
            pass

    def _save_history(self):
        """保存 readline 历史"""
        if not HAS_READLINE:
            return
        try:
            readline.write_history_file(self._history_file)
        except (IOError, PermissionError):
            pass

    # ============================================================
    # 交互循环
    # ============================================================

    def run_interactive(self):
        """交互模式主循环"""
        while True:
            try:
                cmd_line = input('\n⚛ > ').strip()
            except (EOFError, KeyboardInterrupt):
                print("\n👋 Goodbye!")
                self._save_history()
                break

            if not cmd_line:
                continue

            parts = cmd_line.split()
            cmd = parts[0].lower()
            args = parts[1:]

            if cmd in ('quit', 'exit', 'q'):
                print("👋 Goodbye!")
                self._save_history()
                break

            elif cmd == 'help':
                self._print_help()

            elif cmd == 'evolve':
                self.cmd_evolve(args)

            elif cmd == 'matrix':
                self.cmd_matrix(args)

            elif cmd == 'eigenstates':
                self.cmd_eigenstates(args)

            elif cmd == 'plot':
                self.cmd_plot(args)

            elif cmd == 'animate':
                self.cmd_animate(args)

            elif cmd == 'demo':
                self.cmd_demo(args)

            elif cmd == 'status':
                self._print_status()

            elif cmd in ('calc', 'eval', '='):
                self.cmd_calc(' '.join(args))

            else:
                print(f"未知命令: {cmd}。输入 'help' 查看帮助。")

    def cmd_calc(self, expr: str):
        """计算 Python 表达式: calc <expression>

        可直接使用:
            np          — numpy
            x, p        — 当前矩阵力学坐标/动量算符
            a, ad       — 湮灭/产生算符
            H           — 谐振子哈密顿量
            I           — 单位矩阵
            wf          — 当前波函数 psi 数组
            V           — 当前势函数
            grid        — 当前网格

        示例:
            calc np.dot(x, p) - np.dot(p, x)    # [x̂, p̂]
            calc np.trace(H)                     # Tr[Ĥ]
            calc np.linalg.eigvalsh(H)[:5]       # 前5个本征值
            calc np.abs(wf)**2                   # 概率密度
        """
        if not expr:
            print("用法: calc <expression>")
            print("示例: calc np.dot(x, p) - np.dot(p, x)")
            print("      calc np.linalg.eigvalsh(H)[:5]")
            print("      calc np.trace(a @ ad)")
            return

        # 准备命名空间 — 常用模块已预加载，无需 import
        ns = {
            'np': np,
            'numpy': np,
            '__builtins__': {},  # 安全沙箱，禁止 import / exec 等
        }

        # scipy 可选
        try:
            import scipy
            import scipy.linalg
            import scipy.integrate
            ns['sp'] = scipy
            ns['scipy'] = scipy
        except ImportError:
            pass

        # 注入当前状态
        if self.current_matrix is not None:
            ns['x'] = self.current_matrix.x
            ns['p'] = self.current_matrix.p
            ns['a'] = self.current_matrix.a
            ns['ad'] = self.current_matrix.a_dag
            ns['H'] = self.current_matrix.H_harmonic
            ns['I'] = np.eye(self.current_matrix._N)
            ns['N'] = self.current_matrix._N
        if self.current_wf is not None:
            ns['wf'] = self.current_wf.psi
            ns['wf_k'] = self.current_wf.psi_k
        if self.current_potential is not None:
            ns['V'] = self.current_potential
        if self.current_grid is not None:
            ns['grid'] = self.current_grid
            ns['x_arr'] = self.current_grid.x
            ns['k_arr'] = self.current_grid.k

        try:
            result = eval(expr, ns)
            self._display_result(result, expr)
        except Exception as e:
            print(f"Error: {type(e).__name__}: {e}")

    def _display_result(self, result, expr: str):
        """格式化显示计算结果"""
        if isinstance(result, np.ndarray):
            shape = result.shape
            if result.size == 1:
                val = complex(result.item())
                if abs(val.imag) < 1e-14:
                    print(f"= {val.real:.10g}")
                else:
                    print(f"= {val.real:.10g} + {val.imag:.10g}i")
            elif result.ndim == 2 and min(shape) <= 8:
                # 小矩阵：完整显示
                if np.allclose(result.imag, 0):
                    print(np.array2string(result.real, precision=4, suppress_small=True,
                                           max_line_width=100))
                else:
                    print("Real part:")
                    print(np.array2string(result.real, precision=4, suppress_small=True,
                                           max_line_width=100))
                    print("Imag part:")
                    print(np.array2string(result.imag, precision=4, suppress_small=True,
                                           max_line_width=100))
            else:
                print(f"array(shape={shape}, dtype={result.dtype})")
                if result.ndim <= 2:
                    # 显示角落
                    preview = result[:min(4, shape[0]), :min(4, shape[1])] if result.ndim == 2 else result[:8]
                    if np.allclose(preview.imag, 0):
                        print(f"  preview:\n{np.array2string(preview.real, precision=3, suppress_small=True)}")
                    else:
                        print(f"  preview (real):\n{np.array2string(preview.real, precision=3, suppress_small=True)}")
                print(f"  (use calc with slice to inspect: calc x[:5,:5])")
        elif isinstance(result, complex):
            if abs(result.imag) < 1e-14:
                print(f"= {result.real:.10g}")
            else:
                print(f"= {result.real:.10g} + {result.imag:.10g}i")
        elif isinstance(result, (int, float)):
            print(f"= {result:.10g}")
        elif result is None:
            pass
        else:
            print(repr(result))

    def _print_help(self):
        print("""
║  Quantum Agent Commands                                      ║
╠══════════════════════════════════════════════════════════════╣
║  evolve <potential> [...]  — 演化波函数                       ║
║    势函数: harmonic,infinite_well,barrier,double_well,morse  ║
║    示例: evolve harmonic,omega=2.0 x0=1.0 p0=5.0 t_max=10   ║
║                                                              ║
║  matrix [cmd]             — 矩阵力学                          ║
║    示例: matrix report / matrix eigen k=5 / matrix comm x p  ║
║                                                              ║
║  eigenstates <pot> [N]    — 计算本征态                        ║
║    示例: eigenstates harmonic,omega=2.0 5                    ║
║                                                              ║
║  plot <type>              — 绘制图形                          ║
║  animate <type>           — 生成动画                          ║
║  demo <name|all>          — 运行 demo                        ║
║  calc <expression>        — Python 表达式/矩阵运算             ║
║  status                   — 当前状态                          ║
║  help                     — 帮助                              ║
║  quit                     — 退出                              ║
╚══════════════════════════════════════════════════════════════╝
""")

    def _print_status(self):
        """打印当前状态"""
        print("\nCurrent state:")
        if self.current_potential:
            print(f"  Potential: {self.current_potential.label}")
        else:
            print("  Potential: (none)")
        if self.current_grid:
            print(f"  Grid: {self.current_grid}")
        else:
            print("  Grid: (none)")
        if self.current_wf:
            print(f"  Wavefunction: t={self.current_wf.t:.3f}, "
                  f"⟨x⟩={self.current_wf.expectation_x():.4f}, "
                  f"‖ψ‖={self.current_wf.norm:.6f}")
        else:
            print("  Wavefunction: (none)")
        if self.current_result:
            print(f"  Last evolution: {len(self.current_result.times)} snapshots, "
                  f"method={self.current_result.method}")
        print(f"  History: {len(self.history)} commands")


def main():
    parser = argparse.ArgumentParser(description='Quantum Agent - 量子力学智能体')
    parser.add_argument('--demo', type=str, help='Demo 名称 (或 "all")')
    parser.add_argument('--list', action='store_true', help='列出所有 demo')
    parser.add_argument('--config', type=str, help='配置文件路径')
    args = parser.parse_args()

    agent = QuantumAgent(config_path=args.config)

    if args.list:
        print("Available demos:")
        print("  harmonic       - 谐振子波包演化")
        print("  well           - 无限深势阱")
        print("  barrier        - 势垒隧穿")
        print("  matrix         - 矩阵力学")
        print("  double_well    - 双势阱")
        return

    if args.demo:
        agent.cmd_demo([args.demo])
        return

    agent.run_interactive()


if __name__ == '__main__':
    main()

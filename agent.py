#!/usr/bin/env python3
"""
Quantum Agent — QuTiP 风格量子力学智能体
==========================================

交互式量子力学计算和可视化。
所有量子态、算符、动力学函数都可直接在 calc 命令中调用。

用法:
    python agent.py                  # 交互模式
    python agent.py --demo           # 运行演示
    python agent.py --test           # 运行测试
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import readline
    HAS_READLINE = True
except ImportError:
    HAS_READLINE = False


class QuantumAgent:
    """量子力学智能体"""

    def __init__(self):
        self._calc_ns = {}   # calc 持久化变量
        self._init_readline()
        self._welcome()

    def _welcome(self):
        print(r"""
╔══════════════════════════════════════════════════╗
║     Quantum Agent  —  QuTiP-style QM toolkit     ║
║     Fock basis | States | Dynamics | Wigner      ║
╚══════════════════════════════════════════════════╝
Type 'help' for commands, 'demo' to see examples.
""")

    def _init_readline(self):
        if not HAS_READLINE:
            return
        self._hist_file = os.path.expanduser('~/.qm_agent_history')
        readline.set_history_length(1000)
        try:
            readline.read_history_file(self._hist_file)
        except (FileNotFoundError, PermissionError):
            pass
        # Tab 补全
        self._completions = [
            'calc', 'demo', 'test', 'help', 'quit', 'vars',
            'FockBasis', 'coherent', 'coherent_dm', 'squeezed', 'thermal_dm',
            'cat', 'fock', 'fock_dm', 'expect', 'variance', 'g2', 'mandel_q',
            'mean_photon', 'commutator', 'sesolve', 'mesolve', 'steadystate',
            'wigner', 'qfunc', 'plot_wigner', 'plot_photon_dist',
            'WaveGrid', 'gaussian_wavepacket', 'evolve_ssfm', 'animate_wave',
            'fidelity', 'purity', 'photon_dist', 'np',
            'animate', 'plot', 'wigner',
        ]
        readline.set_completer(self._completer)
        readline.parse_and_bind('tab: complete')

    def _completer(self, text, state):
        matches = [c for c in self._completions if c.startswith(text)]
        try:
            return matches[state]
        except IndexError:
            return None

    def _save_hist(self):
        if HAS_READLINE:
            try:
                readline.write_history_file(self._hist_file)
            except (IOError, PermissionError):
                pass

    # ================================================================
    # 命令分发
    # ================================================================

    def run(self):
        while True:
            try:
                line = input('\n⚛ > ').strip()
            except (EOFError, KeyboardInterrupt):
                print("\nGoodbye!")
                self._save_hist()
                break

            if not line:
                continue

            parts = line.split()
            cmd = parts[0].lower()
            args = parts[1:]

            if cmd in ('q', 'quit', 'exit'):
                self._save_hist()
                print("Goodbye!")
                break
            elif cmd == 'help':
                self._help()
            elif cmd == 'demo':
                self._demo()
            elif cmd in ('calc', '=', 'eval'):
                self.calc(' '.join(args))
            elif cmd == 'test':
                self._run_tests()
            elif cmd == 'animate':
                self.calc(f"animate_wave({args[0]}, save_path='{args[1] if len(args)>1 else 'output/animation.mp4'}')" if args else "print('Usage: animate <result_var> [save_path]')")
            elif cmd == 'plot' and args and args[0] == 'wigner':
                self.calc(f"plot_wigner(x, p, W, save='output/wigner.png')" if len(args) < 2 else f"plot_wigner({args[1]}, {args[2] if len(args)>2 else 'p'}, {args[3] if len(args)>3 else 'W'}, save='output/wigner.png')")
            elif cmd == 'wigner':
                self.calc("x, p, W = wigner(psi) if 'psi' in dir() else print('Set psi first: calc psi = coherent(20, 2.0)')")
            else:
                print(f"Unknown: {cmd}.  Type 'help'.")

    # ================================================================
    # calc — Python 表达式求值
    # ================================================================

    def calc(self, expr: str):
        """计算 Python 表达式

        预加载模块和函数可直接使用:
            np, fb, FockBasis
            fock, coherent, squeezed, thermal_dm, cat
            expect, variance, g2, mandel_q, mean_photon
            commutator, sesolve, mesolve, steadystate
            wigner, qfunc, plot_wigner, plot_photon_dist
        """
        if not expr:
            print("Usage: calc <expression>")
            print("       calc <var> = <expression>")
            return

        # 准备命名空间
        import src.qm as qm
        import src.viz as viz

        ns = {
            'np': np, 'numpy': np,
            '__builtins__': {},
            # qm 模块
            'qm': qm, 'FockBasis': qm.FockBasis,
            'fock': qm.fock, 'fock_dm': qm.fock_dm,
            'coherent': qm.coherent, 'coherent_dm': qm.coherent_dm,
            'squeezed': qm.squeezed, 'thermal_dm': qm.thermal_dm,
            'cat': qm.cat,
            'expect': qm.expect, 'variance': qm.variance,
            'g2': qm.g2, 'mandel_q': qm.mandel_q,
            'mean_photon': qm.mean_photon,
            'commutator': qm.commutator,
            'sesolve': qm.sesolve, 'mesolve': qm.mesolve,
            'steadystate': qm.steadystate,
            'photon_dist': qm.photon_dist,
            'fidelity': qm.fidelity, 'purity': qm.purity,
            # wave
            'WaveGrid': qm.WaveGrid,
            'gaussian_wavepacket': qm.gaussian_wavepacket,
            'evolve_ssfm': qm.evolve_ssfm,
            'animate_wave': qm.animate_wave,
            # viz
            'wigner': viz.wigner, 'qfunc': viz.qfunc,
            'plot_wigner': viz.plot_wigner,
            'plot_photon_dist': viz.plot_photon_dist,
            # 全局 fb
            'fb': qm.FockBasis(50),
        }

        # qft 量子场论模块 (懒加载)
        try:
            import src.qft as qft
            ns['qft'] = qft
            ns['ScalarField'] = qft.ScalarField
        except Exception:
            pass
        ns.update(self._calc_ns)

        # 特殊命令
        if expr.strip() == 'vars':
            if self._calc_ns:
                print("Variables:")
                for k, v in sorted(self._calc_ns.items()):
                    s = f"array({v.shape}, {v.dtype})" if isinstance(v, np.ndarray) else repr(v)[:60]
                    print(f"  {k}: {s}")
            else:
                print("(no variables)")
            return

        # 赋值或求值
        if '=' in expr and not expr.startswith('='):
            try:
                exec(expr, ns)
                for k, v in ns.items():
                    if k not in ('np', 'numpy', 'qm', 'fb', '__builtins__') and not k.startswith('_'):
                        if k not in self._calc_ns or self._calc_ns.get(k) is not v:
                            self._calc_ns[k] = v
                vname = expr.split('=')[0].strip()
                if vname in self._calc_ns:
                    self._show(self._calc_ns[vname], prefix=f"  {vname} = ")
            except Exception as e:
                print(f"Error: {e}")
        else:
            try:
                result = eval(expr, ns)
                self._show(result)
            except Exception as e:
                print(f"Error: {e}")

    def _show(self, val, prefix=""):
        if isinstance(val, np.ndarray):
            if val.size == 1:
                c = complex(val.item())
                print(f"{prefix}{c.real:.6g}" if abs(c.imag) < 1e-12
                      else f"{prefix}{c.real:.6g}+{c.imag:.6g}i")
            elif val.ndim == 1 and val.size <= 12:
                print(f"{prefix}{np.array2string(val, precision=4, suppress_small=True)}")
            elif val.ndim == 2 and min(val.shape) <= 6:
                print(np.array2string(val.real if np.allclose(val.imag, 0) else val,
                                      precision=3, suppress_small=True))
            else:
                print(f"{prefix}array(shape={val.shape}, dtype={val.dtype})")
                if val.ndim <= 2:
                    p = val[:min(3, val.shape[0]), :min(3, val.shape[1])] if val.ndim == 2 else val[:6]
                    print(f"  preview: {np.array2string(p.real, precision=2, suppress_small=True)}")
        elif isinstance(val, (int, float, complex)):
            print(f"{prefix}{val:.6g}")
        elif isinstance(val, dict):
            print(f"{prefix}{{k: type(v).__name__ for k, v in val.items()}}")
        elif val is not None:
            print(f"{prefix}{repr(val)[:120]}")

    # ================================================================
    # demo / help / test
    # ================================================================

    def _demo(self):
        import src.qm as qm
        fb = qm.FockBasis(30)

        print("\n" + "=" * 55)
        print("  Quantum Agent Demo")
        print("=" * 55)

        # 1. 相干态
        alpha = 2.0 + 0j
        psi = qm.coherent(30, alpha)
        print(f"\n1. Coherent |α={alpha}⟩")
        print(f"   ⟨n⟩ = {qm.mean_photon(psi, fb):.3f}  (|α|² = {abs(alpha)**2:.1f})")
        print(f"   g²(0) = {qm.g2(psi, fb):.4f}  (Poisson: 1.0)")

        # 2. 热态
        rho_th = qm.thermal_dm(30, 1.0)
        print(f"\n2. Thermal ⟨n⟩=1.0")
        print(f"   ⟨n⟩ = {qm.mean_photon(rho_th, fb):.3f}")
        print(f"   g²(0) = {qm.g2(rho_th, fb):.4f}  (Bunched: 2.0)")

        # 3. 压缩态
        r = 0.5
        psi_sq = qm.squeezed(30, r)
        print(f"\n3. Squeezed vacuum r={r}")
        print(f"   ⟨n⟩ = {qm.mean_photon(psi_sq, fb):.4f}  (sinh²(r) = {np.sinh(r)**2:.4f})")
        print(f"   g²(0) = {qm.g2(psi_sq, fb):.4f}")

        # 4. 猫态
        psi_cat = qm.cat(30, 2.0, 0.0)
        print(f"\n4. Even cat |α=2⟩+|−α=2⟩")
        print(f"   ⟨n⟩ = {qm.mean_photon(psi_cat, fb):.3f}")

        # 5. 对易子
        xp = qm.commutator(fb.x, fb.p)
        print(f"\n5. [x̂, p̂] submatrix norm: {np.linalg.norm(xp[:5,:5] - 1j*np.eye(5), 'fro'):.2e}")

        # 6. 衰减
        print(f"\n6. Lindblad decay demo:")
        H = fb.hamiltonian()
        rho0 = qm.coherent_dm(30, 2.0)
        tlist = np.linspace(0, 5, 30)
        gamma = 0.2
        res = qm.mesolve(H, rho0, tlist, c_ops=[np.sqrt(gamma)*fb.a],
                         e_ops=[fb.n_op])
        n_t = np.real(res['expect'][0])
        print(f"   ⟨n⟩(0) = {n_t[0]:.2f}  →  ⟨n⟩(5) = {n_t[-1]:.3f}")

        # 7. Wigner
        print(f"\n7. Try in calc:")
        print(f"   calc psi = coherent(20, 1+0.5j)")
        print(f"   calc x, p, W = wigner(psi)")
        print(f"   calc plot_wigner(x, p, W)")

        print("\n" + "=" * 55)

    def _help(self):
        print("""
Commands:
  calc <expr>        Evaluate Python expression
  calc <var> = <expr> Assign variable
  calc vars           List variables
  animate <var> [path] Animate wavefunction result
  plot wigner [x] [p] [W]  Plot Wigner function
  wigner              Quick Wigner of current psi
  demo                Run Fock-basis demonstration
  test                Run self-tests
  help                This help
  quit                Exit

Preloaded in calc:
  FockBasis, fock, coherent, squeezed, thermal_dm, cat
  expect, variance, g2, mandel_q, mean_photon
  commutator, sesolve, mesolve, steadystate
  wigner, qfunc, plot_wigner, plot_photon_dist
  fb (default FockBasis(50)), np (numpy)
""")

    def _run_tests(self):
        import src.qm as qm
        passed = 0
        total = 0

        def check(name, cond):
            nonlocal passed, total
            total += 1
            if cond:
                passed += 1
                print(f"  ✓ {name}")
            else:
                print(f"  ✗ {name}")

        fb = qm.FockBasis(20)

        check("fock norm", abs(np.linalg.norm(qm.fock(20, 5)) - 1) < 1e-12)
        check("coherent norm", abs(np.linalg.norm(qm.coherent(20, 2+1j)) - 1) < 1e-12)

        psi = qm.coherent(20, 2.0)
        check("mean_photon ≈ |α|²", abs(qm.mean_photon(psi, fb) - 4.0) < 0.1)
        check("g2 coherent ≈ 1", abs(qm.g2(psi, fb) - 1.0) < 0.1)

        rho_th = qm.thermal_dm(20, 1.0)
        check("mean_photon thermal", abs(qm.mean_photon(rho_th, fb) - 1.0) < 0.1)
        check("g2 thermal ≈ 2", abs(qm.g2(rho_th, fb) - 2.0) < 0.3)

        xp = qm.commutator(fb.x, fb.p)
        k = 15
        check("[x,p] ≈ iI",
              np.linalg.norm(xp[:k, :k] - 1j * np.eye(k), 'fro') < 0.5)

        check("purity fock=1", abs(qm.purity(qm.fock_dm(20, 3)) - 1.0) < 1e-12)
        check("purity thermal<1", qm.purity(rho_th) < 0.99)

        print(f"\n  {passed}/{total} passed")


def main():
    import argparse
    p = argparse.ArgumentParser(description='Quantum Agent')
    p.add_argument('--demo', action='store_true')
    p.add_argument('--test', action='store_true')
    p.add_argument('--list', action='store_true')
    args = p.parse_args()

    agent = QuantumAgent()
    if args.list:
        print("Available demos:")
        demos = [
            ('heisenberg_uncertainty', 'Δx·Δp ≥ ℏ/2'),
            ('free_particle', 'Free particle spreading'),
            ('measurement_collapse', 'Position measurement collapse'),
            ('momentum_collapse', 'Momentum measurement collapse'),
            ('energy_collapse', 'Energy measurement collapse'),
            ('double_slit', 'Double-slit interference (2D TDSE)'),
            ('quantum_eraser', 'Quantum eraser experiment'),
        ]
        for name, desc in demos:
            print(f"  {name:<28s} {desc}")
        print(f"\nRun: python demos/<name>.py")
        return
    if args.demo:
        agent._demo()
    elif args.test:
        agent._run_tests()
    else:
        agent.run()


if __name__ == '__main__':
    main()

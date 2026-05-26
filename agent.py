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
        # Tab 补全: 命令 + 智能文件路径
        self._completions = [
            'calc', 'demo', 'test', 'help', 'quit', 'vars',
            'cd', 'pwd', 'ls', 'run', 'animate', 'plot', 'wigner', 'formula',
            'FockBasis', 'coherent', 'coherent_dm', 'squeezed', 'thermal_dm',
            'cat', 'fock', 'fock_dm', 'expect', 'variance', 'g2', 'mandel_q',
            'mean_photon', 'commutator', 'sesolve', 'mesolve', 'steadystate',
            'wigner', 'qfunc', 'plot_wigner', 'plot_photon_dist',
            'WaveGrid', 'gaussian_wavepacket', 'evolve_ssfm', 'animate_wave', 'double_well', 'periodic_potential', 'delta_barrier', 'finite_well',
            'fidelity', 'purity', 'photon_dist', 'np',
            'ScalarField', 'LatticePhi4',
        ]
        readline.set_completer_delims(' \t\n`~!@#$%^&*()-=+[{]}\\|;:\'",<>?')
        readline.set_completer(self._completer)
        readline.parse_and_bind('tab: complete')

    def _completer(self, text, state):
        """Smart completion: command names + file paths after run/animate"""
        import glob as _glob

        line = readline.get_line_buffer()
        stripped = line.lstrip()
        words = stripped.split()
        cmd = words[0].lower() if words else ''

        # File-path completion after 'run', 'animate', 'cd'
        if cmd in ('run', 'animate', 'cd') and len(words) >= 1:
            if not line.rstrip().endswith(cmd) or len(words) > 1:
                prefix = text if text else '.'
                agent_dir = os.path.dirname(__file__)
                matches = []
                for base in ['scripts', 'demos', '.', 'output']:
                    d = os.path.join(agent_dir, base)
                    if os.path.isdir(d):
                        for p in _glob.glob(os.path.join(d, prefix + '*')):
                            rel = os.path.relpath(p, agent_dir)
                            if os.path.isdir(p):
                                rel += os.sep
                            if rel.startswith(prefix):
                                matches.append(rel)
                # Also search user-provided paths (./ ../ ~/)
                if prefix.startswith(('.', os.sep, '~')):
                    pat = os.path.expanduser(prefix + '*')
                    for p in _glob.glob(pat):
                        s = p + os.sep if os.path.isdir(p) else p
                        if s not in matches:
                            matches.append(s)
                # cd/run/animate: also search from cwd for any prefix
                if cmd in ('run', 'animate', 'cd'):
                    cwd = os.getcwd()
                    for p in _glob.glob(os.path.join(cwd, prefix + '*')):
                        if os.path.isdir(p):
                            s = os.path.basename(p) + os.sep
                        else:
                            s = os.path.basename(p)
                        if s.startswith(prefix) and s not in matches:
                            matches.append(s)
                matches = sorted(set(matches))
                if state < len(matches):
                    return matches[state]
                return None

        # Default: command/variable-name completion
        matches = [c for c in self._completions if c.startswith(text)]
        for v in self._calc_ns:
            if v.startswith(text) and v not in matches:
                matches.append(v)
        matches = sorted(matches)
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
    # LaTeX → Unicode 转换表
    # ================================================================

    @staticmethod
    def _latex_to_unicode(latex_str: str) -> str:
        """将 LaTeX 数学公式转换为 Unicode 字符，直接在终端显示"""
        import re

        _LATEX_MAP = {
            # Greek
            r'\hbar': 'ℏ', r'\Psi': 'Ψ', r'\psi': 'ψ',
            r'\alpha': 'α', r'\beta': 'β', r'\gamma': 'γ', r'\delta': 'δ',
            r'\epsilon': 'ε', r'\zeta': 'ζ', r'\eta': 'η', r'\theta': 'θ',
            r'\iota': 'ι', r'\kappa': 'κ', r'\lambda': 'λ', r'\mu': 'μ',
            r'\nu': 'ν', r'\xi': 'ξ', r'\pi': 'π', r'\rho': 'ρ',
            r'\sigma': 'σ', r'\tau': 'τ', r'\upsilon': 'υ', r'\phi': 'φ',
            r'\chi': 'χ', r'\omega': 'ω',
            r'\Gamma': 'Γ', r'\Delta': 'Δ', r'\Theta': 'Θ', r'\Lambda': 'Λ',
            r'\Xi': 'Ξ', r'\Pi': 'Π', r'\Sigma': 'Σ', r'\Upsilon': 'Υ',
            r'\Phi': 'Φ', r'\Omega': 'Ω',
            # Math operators
            r'\partial': '∂', r'\nabla': '∇', r'\int': '∫', r'\sum': '∑',
            r'\prod': '∏', r'\infty': '∞', r'\approx': '≈', r'\propto': '∝',
            r'\sim': '~', r'\simeq': '≃', r'\times': '×', r'\cdot': '·',
            r'\otimes': '⊗', r'\oplus': '⊕', r'\ominus': '⊖',
            r'\langle': '⟨', r'\rangle': '⟩', r'\bra': '⟨', r'\ket': '|',
            r'\equiv': '≡', r'\neq': '≠', r'\pm': '±', r'\mp': '∓',
            r'\leq': '≤', r'\geq': '≥', r'\ll': '≪', r'\gg': '≫',
            r'\rightarrow': '→', r'\Rightarrow': '⇒', r'\leftarrow': '←',
            r'\uparrow': '↑', r'\downarrow': '↓', r'\leftrightarrow': '↔',
            r'\mapsto': '↦', r'\longrightarrow': '→',
            r'\dagger': '†', r'\ddagger': '‡',
            r'\dots': '…', r'\cdots': '⋯', r'\vdots': '⋮',            r'\ddots': '⋱',
            r'\xrightarrow': '→', r'\xleftarrow': '←',
            r'\mapsto': '↦', r'\to': '→', r'\implies': '⇒',
            r'\gg': '≫', r'\ll': '≪', r'\ggg': '⋙', r'\lll': '⋘',
            r'\mid': '∣', r'\nmid': '∤',
            r'\sinh': 'sinh', r'\cosh': 'cosh', r'\tanh': 'tanh',
            # Hatted operators (precomposed chars for better terminal support)
            r'\hat{H}': 'Ĥ', r'\hat{x}': 'x̂', r'\hat{p}': 'p̂',
            r'\hat{a}': 'â', r'\hat{N}': 'N̂', r'\hat{\rho}': 'ρ̂',
            # Superscripts (>2 chars first)
            r'^{(0)}': '⁽⁰⁾', r'^0': '⁰', r'^1': '¹', r'^2': '²',
            r'^3': '³', r'^4': '⁴', r'^5': '⁵', r'^6': '⁶',
            r'^7': '⁷', r'^8': '⁸', r'^9': '⁹', r'^+': '⁺', r'^-': '⁻',
            r'^*': '*',  # complex conjugate
            # Subscripts
            r'_0': '₀', r'_1': '₁', r'_2': '₂', r'_3': '₃', r'_4': '₄',
            r'_5': '₅', r'_6': '₆', r'_7': '₇', r'_8': '₈', r'_9': '₉',
        }

        s = latex_str
        # Longest patterns first
        for k, v in sorted(_LATEX_MAP.items(), key=lambda x: -len(x[0])):
            s = s.replace(k, v)

        # \\frac{num}{den}
        def frac_repl(m):
            num = QuantumAgent._latex_to_unicode(m.group(1))
            den = QuantumAgent._latex_to_unicode(m.group(2))
            # Single char or number: use compact a/b
            if len(num) <= 2 and not (' ' in num or '(' in num):
                if len(den) <= 2 and not (' ' in den or '(' in den):
                    return f'{num}/{den}'
            return f'({num})/({den})'
        s = re.sub(
            r'\\frac\{([^{}]*(?:\{[^{}]*\}[^{}]*)*?)\}\{([^{}]*(?:\{[^{}]*\}[^{}]*)*?)\}',
            frac_repl, s)

        # \\sqrt{...}
        s = re.sub(
            r'\\sqrt\{([^{}]*(?:\{[^{}]*\}[^{}]*)*?)\}',
            lambda m: f'√({QuantumAgent._latex_to_unicode(m.group(1))})', s)

        # Cleanup
        s = s.replace('\\left', '').replace('\\right', '')
        s = s.replace('{', '').replace('}', '')
        s = s.replace('\\,', ' ').replace('\\;', '  ').replace('\\ ', ' ')
        return s

    # ================================================================
    # formula — LaTeX 公式终端显示
    # ================================================================

    def _render_formula(self, latex: str):
        """在终端显示 LaTeX 公式 (Unicode) + 保存 PNG 文件

        用法:
            formula i\\hbar\\frac{\\partial}{\\partial t}\\Psi = \\hat{H}\\Psi
            formula \\Delta x \\cdot \\Delta p \\geq \\frac{\\hbar}{2}

        Unicode 显示在终端，PNG 保存到 output/formulas/
        """
        if not latex:
            print("Usage: formula <LaTeX expression>")
            print(r"Example: formula i\hbar\frac{\partial}{\partial t}\Psi = \hat{H}\Psi")
            return

        # 先显示 Unicode 版本
        unicode_version = self._latex_to_unicode(latex)
        print(f'\n  {unicode_version}\n')

        # 渲染 PNG 保存
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import time

            fig, ax = plt.subplots(figsize=(10, 1.0), facecolor='#ffffff')
            ax.axis('off')
            try:
                ax.text(0.5, 0.5, f'${latex}$', transform=ax.transAxes,
                        fontsize=30, ha='center', va='center', color='#1f2328')
            except Exception:
                ax.text(0.5, 0.5, latex, transform=ax.transAxes,
                        fontsize=16, ha='center', va='center', color='#1f2328',
                        fontfamily='monospace')

            save_dir = os.path.join(os.path.dirname(__file__), 'output', 'formulas')
            os.makedirs(save_dir, exist_ok=True)
            ts = time.strftime('%H%M%S')
            save_path = os.path.join(save_dir, f'formula_{ts}.png')
            fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='#ffffff')
            plt.close(fig)
            print(f'  PNG: {save_path}')
        except Exception:
            pass  # PNG save is best-effort

    # ================================================================
    # run — 执行脚本 (.qms 文件)
    # ================================================================

    def _run_script(self, script_path: str):
        """执行 .qms 量子脚本文件

        类似 MATLAB .m 文件: 每行一条命令，支持:
          - Python 表达式 (自动求值)
          - formula / animate / plot wigner / demo / test / help
          - run <another_script> (嵌套调用)
          - # 注释
          - 变量跨行共享 (同一命名空间)
        """
        if not script_path:
            print("Usage: run <script.qms>")
            print("       python agent.py --run <script.qms>")
            return

        import os as _os
        # Resolve relative path: cwd first, then agent_dir
        if not _os.path.isabs(script_path):
            cwd_path = _os.path.join(_os.getcwd(), script_path)
            if _os.path.exists(cwd_path):
                script_path = cwd_path
            elif not _os.path.exists(script_path):
                script_path = _os.path.join(_os.path.dirname(__file__), script_path)

        if not _os.path.exists(script_path):
            print(f"Script not found: {script_path}")
            return

        script_dir = _os.path.dirname(_os.path.abspath(script_path))
        print(f"  [Running: {script_path}]")
        try:
            with open(script_path, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"  Error reading script: {e}")
            return

        # 切换工作目录到脚本所在目录 (便于 output/animation.mp4 等相对路径)
        old_cwd = _os.getcwd()
        _os.chdir(script_dir)

        import ast as _ast

        try:
            buf_lines = []  # 多行语句缓冲区
            lineno = 1
            raw_idx = 0
            while raw_idx < len(lines):
                raw_line = lines[raw_idx]
                raw_idx += 1
                line = raw_line.strip()

                # 空行和注释
                if not line or line.startswith('#'):
                    if line.startswith('#'):
                        print(f"  #{line[1:].strip()}")
                    lineno += 1
                    continue

                # Known commands: dispatch immediately (flush buffer first)
                _KNOWN = ('formula', 'run', 'cd', 'ls', 'pwd', 'help', 'demo',
                          'test', 'animate', 'plot', 'wigner', 'q', 'quit', 'exit')
                first_word = line.split()[0].lower() if line else ''
                if first_word in _KNOWN:
                    # Flush pending buffer
                    if buf_lines:
                        buf = chr(10).join(buf_lines)
                        disp = buf[:70] + ('...' if len(buf) > 70 else '')
                        print(f"  [{lineno - len(buf_lines) + 1}] {disp}")
                        try:
                            self._dispatch(buf)
                        except SystemExit:
                            raise
                        except Exception as e:
                            print(f"  Error: {e}")
                        buf_lines = []
                    # Dispatch this command
                    disp = line[:70] + ('...' if len(line) > 70 else '')
                    print(f"  [{lineno}] {disp}")
                    try:
                        self._dispatch(line)
                    except SystemExit:
                        raise
                    except Exception as e:
                        print(f"  Error (line {lineno}): {e}")
                    lineno += 1
                    continue

                # 多行语句: 首行 strip, 后续行保留原始缩进
                if not buf_lines:
                    buf_lines.append(line)
                else:
                    buf_lines.append(raw_line.rstrip())
                buf = chr(10).join(buf_lines)

                # 检查语句是否完整 (用 ast.parse 试解析)
                complete = False
                try:
                    _ast.parse(buf, mode='exec')
                    complete = True
                    # 复合语句: 缩进未退回时继续累积
                    _first = buf_lines[0].lstrip()
                    if _first.startswith(('def ', 'class ', 'if ', 'elif ',
                        'else:', 'for ', 'while ', 'with ', 'try:', 'except ',
                        'except:', 'finally:')):
                        # 从当前位置向前找到下一个非空非注释行
                        peek_idx = raw_idx
                        while peek_idx < len(lines):
                            peek_line = lines[peek_idx].rstrip(chr(10))
                            if peek_line and not peek_line.lstrip().startswith('#'):
                                break
                            peek_idx += 1
                        if peek_idx < len(lines):
                            next_line = lines[peek_idx]
                            def_indent = len(buf_lines[0]) - len(buf_lines[0].lstrip())
                            next_indent = len(next_line) - len(next_line.lstrip())
                            if next_indent > def_indent:
                                complete = False
                except SyntaxError:
                    # 不完整语句 (unclosed brackets / trailing colon), 继续累积
                    if raw_idx < len(lines):
                        lineno += 1
                        continue

                if complete:
                    disp = buf[:70] + ('...' if len(buf) > 70 else '')
                    print(f"  [{lineno - len(buf_lines) + 1}] {disp}")
                    try:
                        self._dispatch(buf)
                    except SystemExit:
                        raise
                    except Exception as e:
                        print(f"  Error: {e}")
                    buf_lines = []
                lineno += 1

            # 尾残 (如文件以不完整语句结束)
            if buf_lines:
                buf = '\n'.join(buf_lines)
                print(f"  [{lineno - len(buf_lines) + 1}] {buf[:70]}...")
                try:
                    self._dispatch(buf)
                except SystemExit:
                    raise
                except Exception as e:
                    print(f"  Error: {e}")
        finally:
            _os.chdir(old_cwd)

        print(f"  [Done: {script_path}]")

    # ================================================================
    # cd / pwd — 路径导航
    # ================================================================

    def _cd(self, args):
        """cd [path] — change working directory, tab-completes paths"""
        path = args[0] if args else os.path.expanduser('~')
        path = os.path.expanduser(path)
        try:
            os.chdir(path)
            print(f'  {os.getcwd()}')
        except (FileNotFoundError, NotADirectoryError, PermissionError) as e:
            print(f'  cd: {e}')

    # ================================================================
    # ls — 列出当前目录文件
    # ================================================================

    def _ls(self):
        """列出当前目录的文件和子目录"""
        cwd = os.getcwd()
        try:
            entries = sorted(os.listdir(cwd))
        except PermissionError as e:
            print(f'  ls: {e}')
            return
        if not entries:
            print('  (empty)')
            return
        # 列排版
        max_len = max(len(e) for e in entries) + 2
        cols = max(1, 80 // max_len)
        rows = (len(entries) + cols - 1) // cols
        for r in range(rows):
            line_parts = []
            for c in range(cols):
                idx = r + c * rows
                if idx < len(entries):
                    e = entries[idx]
                    if os.path.isdir(os.path.join(cwd, e)):
                        e = e + '/'
                    line_parts.append(e.ljust(max_len))
            print('  ' + ''.join(line_parts).rstrip())

    # ================================================================
    # 命令分发
    # ================================================================

    def _dispatch(self, line: str) -> bool:
        """分发单条命令。返回 False 表示退出。"""
        parts = line.split()
        cmd = parts[0].lower()
        args = parts[1:]

        if cmd in ('q', 'quit', 'exit'):
            return False
        elif cmd == 'help':
            if args:
                self._function_help(args[0].lower())
            else:
                self._help()
        elif cmd == 'demo':
            self._demo()
        elif cmd == 'cd':
            self._cd(args)
        elif cmd == 'pwd':
            print(os.getcwd())
        elif cmd == 'ls':
            self._ls()
        elif cmd in ('calc', '=', 'eval'):
            self.calc(' '.join(args))
        elif cmd == 'test':
            self._run_tests()
        elif cmd == 'animate':
            if args:
                save = args[1] if len(args) > 1 else 'output/animation.mp4'
                self.calc(f"animate_wave({args[0]}, save_path='{save}')")
            else:
                print('Usage: animate <result_var> [save_path]')
        elif cmd == 'plot' and args and args[0] == 'wigner':
            xv = args[1] if len(args) > 1 else 'x'
            pv = args[2] if len(args) > 2 else 'p'
            wv = args[3] if len(args) > 3 else 'W'
            self.calc(f"plot_wigner({xv}, {pv}, {wv}, save='output/wigner.png')")
        elif cmd == 'wigner':
            self.calc("x, p, W = wigner(psi) if 'psi' in dir() else print('Set psi first')")
        elif cmd == 'formula':
            formula_args = line[len('formula'):].strip()
            self._render_formula(formula_args)
        elif cmd == 'run':
            self._run_script(' '.join(args))
        else:
            # 不是已知命令, 尝试作为 Python 表达式求值
            self.calc(line)
        return True

    def run(self):
        while True:
            try:
                cwd = os.getcwd()
                home = os.path.expanduser('~')
                if cwd.startswith(home):
                    cwd = '~' + cwd[len(home):]
                if len(cwd) > 35:
                    cwd = '...' + cwd[-32:]
                prompt = f'\n⚛ {cwd} > '
                line = input(prompt).strip()
            except (EOFError, KeyboardInterrupt):
                print("\nGoodbye!")
                self._save_hist()
                break

            if not line:
                continue

            if not self._dispatch(line):
                self._save_hist()
                print("Goodbye!")
                break

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
            '__builtins__': {'__import__': __import__, 'print': print, 'abs': abs, 'len': len, 'range': range,
                             'int': int, 'float': float, 'complex': complex,
                             'bool': bool, 'str': str, 'list': list, 'dict': dict,
                             'tuple': tuple, 'set': set, 'min': min, 'max': max,
                             'sum': sum, 'round': round, 'zip': zip, 'enumerate': enumerate,
                             'sorted': sorted, 'reversed': reversed, 'isinstance': isinstance},
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
            'double_well': qm.double_well,
            'periodic_potential': qm.periodic_potential,
            'delta_barrier': qm.delta_barrier,
            'finite_well': qm.finite_well,
            'harmonic_oscillator_potential': qm.harmonic_oscillator_potential,
            'step_potential': qm.step_potential,
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
            ns['LatticePhi4'] = qft.LatticePhi4
            ns['wick_expand'] = qft.wick_expand
            ns['feynman_amplitude_phi4_2to2'] = qft.feynman_amplitude_phi4_2to2
            ns['differential_cross_section'] = qft.differential_cross_section
            ns['transition_probability'] = qft.transition_probability
            ns['PathIntegralMC'] = qft.PathIntegralMC
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

        # 支持 import 语句 (脚本中常用)
        if expr.strip().startswith('import ') or expr.strip().startswith('from '):
            try:
                exec(expr, ns)
                # 追踪新增的模块/变量
                for k, v in ns.items():
                    if k not in ('np', 'numpy', 'qm', 'fb', '__builtins__') and not k.startswith('_'):
                        if k not in self._calc_ns or self._calc_ns.get(k) is not v:
                            self._calc_ns[k] = v
            except Exception as e:
                print(f"Error: {e}")
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
            # Try eval first, fall back to exec for statements (for, if, etc.)
            try:
                result = eval(expr, ns)
                self._show(result)
            except SyntaxError:
                try:
                    exec(expr, ns)
                except Exception as e2:
                    print(f"Error: {e2}")
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

    def _function_help(self, name):
        '''显示函数的物理公式'''
        formulas = {
            'coherent': '|α⟩ = e^{-|α|²/2} Σ αⁿ/√(n!) |n⟩    ⟨n⟩ = |α|²,  g²=1',
            'squeezed': '|ζ⟩ = exp[½(ζ*a² - ζ a†²)]|0⟩    ⟨n⟩ = sinh²(r),  g²=3+1/sinh²(r)',
            'thermal_dm': 'ρ_th = Σ n_thⁿ/(1+n_th)^{n+1} |n⟩⟨n|    g²=2',
            'cat': '|cat⟩ ∝ |α⟩ + e^{iφ}|-α⟩    W<0 (non-classical)',
            'fock': '|n⟩ = (a†)^n/√(n!) |0⟩    g² = 1 - 1/n (sub-Poisson)',
            'g2': 'g²(0) = ⟨a†a†aa⟩/⟨a†a⟩²     classical: g²≥1, quantum: g²<1 possible',
            'mandel_q': 'Q = ⟨n⟩(g²-1)    Q=0 Poisson, Q<0 sub-Poisson (non-classical)',
            'commutator': '[A,B] = AB - BA    [x̂,p̂] = iħ    [a,a†] = 1',
            'expect': '⟨O⟩ = Tr[ρ O] (DM) or ⟨ψ|O|ψ⟩ (pure)',
            'variance': 'ΔO² = ⟨O²⟩ - ⟨O⟩²    Heisenberg: Δx·Δp ≥ ħ/2',
            'sesolve': 'iħ ∂|ψ⟩/∂t = H|ψ⟩    exact diagonalization',
            'mesolve': 'dρ/dt = -i[H,ρ] + Σ γ(LρL† - ½{L†L,ρ})    RK4 integrator',
            'steadystate': 'solve L[ρ_ss] = 0    Liouvillian superoperator',
            'wigner': 'W(x,p) = 1/πħ ∫ ⟨x+y|ρ|x-y⟩ e^{-2ipy/ħ} dy    W<0 = non-classical',
            'wignerg': 'W negative values = quantumness signature',
            'gaussian_wavepacket': 'ψ(x,0) = (πσ²)^{-1/4} exp[-(x-x₀)²/2σ² + ip₀x/ħ]',
            'evolve_ssfm': 'iħ ∂ψ/∂t = -ħ²/2m ∂²ψ/∂x² + V(x)ψ    Split-Step Fourier',
            'wavegrid': '1D grid: x ∈ [x_min, x_max], N points, dx, k-space',
            'scalarfield': 'φ̂(x) = ∫ d³k (â_k e^{-ikx} + â†_k e^{ikx}) / √(2ω_k)',
            'latticephi4': 'H = Σ [½π² + ½(∇φ)² + ½m²φ² + λφ⁴/4!]    exact diag',
            'pathintegralmc': 'Z = ∫ Dx e^{-S_E[x]/ħ}    Metropolis sampling, β → ∞ = ground state',
            'double_well': 'V(x)=V₀[(x/a)²-1]²/2    quantum tunneling between wells',
            'delta_barrier': 'V(x)=g/√(πσ) e^{-(x-x₀)²/σ²}    narrow Gaussian ≈ δ(x-x₀)',
            'periodic_potential': 'V(x)=A cos(2πx/λ)    optical lattice / Bloch waves',
            'step_potential': 'V(x)=V₀·θ(x-x₀)    scattering: transmission/reflection',
            'finite_well': 'V(x)=-V₀ (|x-x₀|<w/2)    bound states, discrete spectrum',
            'fidelity': 'F=|⟨ψ₁|ψ₂⟩|²    state overlap, F∈[0,1]',
            'purity': 'P=Tr[ρ²]    P=1 pure, P<1 mixed',
        }
        if name in formulas:
            print(f'\n  {name}:\n    {formulas[name]}\n')
        else:
            print(f'  No formula reference for: {name}')
            print(f'  Try: {", ".join(sorted(formulas.keys())[:10])}...')

    def _help(self):
        print("""
Commands:
  <expression>        Direct Python expression (auto-evaluated)
  calc <var> = <expr> Assign variable (prefix optional)
  calc vars           List variables
  cd [path]           Change working directory (tab-completes)
  pwd                 Print working directory
  ls                  List files in current directory
  formula <latex>     Render LaTeX to Unicode in terminal
  animate <var> [path] Animate wavefunction result
  plot wigner [x] [p] [W]  Plot Wigner function
  wigner              Quick Wigner of current psi
  run <script.qms>    Execute quantum script file
  demo                Run Fock-basis demonstration
  test                Run self-tests
  help                This help
  quit                Exit

Scripts (.qms files):
  python agent.py --run scripts/harmonic.qms
  (or inside agent:  run scripts/harmonic.qms)

Physics reference (type 'help <function>'):
  help coherent     g2        wigner    sesolve
  help mesolve      commutator     expect

Preloaded:
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
    p.add_argument('--run', metavar='SCRIPT', help='Execute a .qms script file', dest='script')
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
            ('qft_scalar_field', 'QFT: free scalar field'),
            ('qft_lattice', 'QFT: lattice φ⁴ theory'),
            ('qft_scattering', 'QFT: Feynman diagrams'),
        ]
        for name, desc in demos:
            print(f"  {name:<28s} {desc}")
        print(f"\nRun: python demos/<name>.py")
        print(f"\nScript: python agent.py --run scripts/<name>.qms")
        return
    if args.script:
        agent._run_script(args.script)
    elif args.demo:
        agent._demo()
    elif args.test:
        agent._run_tests()
    else:
        agent.run()


if __name__ == '__main__':
    main()

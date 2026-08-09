#!/usr/bin/env python3
"""Quick validation of gauge.py, dirac.py, qed.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np

passed = 0
failed = 0

def check(name, cond):
    global passed, failed
    if cond:
        passed += 1
        print(f'  ✓ {name}')
    else:
        failed += 1
        print(f'  ✗ {name}')

# ================================================================
# gauge.py
# ================================================================
print("── gauge.py ──")
from src.qft.gauge import (
    minkowski_metric, lorentz_dot, k_squared,
    polarization_vectors, check_polarization,
    photon_propagator, photon_propagator_landau,
    GaugeField, field_strength_tensor, ward_identity_check, G_MUNU,
)

g = minkowski_metric()
check('minkowski_metric diag(1,-1,-1,-1)', g[0,0]==1 and g[1,1]==-1)

k_light = np.array([5.0, 3.0, 0.0, 4.0])
check('k²=0 for lightlike', abs(k_squared(k_light)) < 1e-12)

k_sp = np.array([3.0, 1.0, 1.0, 1.0])
check('k²=6 for spacelike', abs(k_squared(k_sp) - 6.0) < 1e-12)

kz = np.array([2.0, 0.0, 0.0, 2.0])
c = check_polarization(kz)
check('k·ε₁=0', c['k·ε₁'] < 1e-10)
check('k·ε₂=0', c['k·ε₂'] < 1e-10)
check('ε₁²=-1', abs(c['ε₁²'] + 1) < 1e-10)
check('ε₂²=-1', abs(c['ε₂²'] + 1) < 1e-10)

D = photon_propagator(1.0)
check('Feynman propagator shape', D.shape == (4, 4))

Dl = photon_propagator_landau(kz)
check('Landau propagator shape', Dl.shape == (4, 4))

gf = GaugeField(volume=1.0)
A = gf.field_at(np.array([1.0, 0.5, 0.0, 0.0]))
check('GaugeField A_μ returns (4,)', A.shape == (4,))

ward = ward_identity_check(lambda kk: polarization_vectors(kk, 1), kz)
check('Ward identity passes', ward['passes'])

# Field strength placeholder
def A_fn(x): return np.array([np.sin(x[0]), 0., 0., 0.])
F = field_strength_tensor(A_fn, np.array([0., 0., 0., 0.]))
check('field_strength_tensor shape', F.shape == (4, 4))

# ================================================================
# dirac.py
# ================================================================
print("── dirac.py ──")
from src.qft.dirac import (
    GammaMatrices, DiracSpinor, dirac_slash,
    spin_sum_u, spin_sum_v,
    spin_sum_u_from_spinors, spin_sum_v_from_spinors,
    dirac_equation_check, I4,
)

gm = GammaMatrices(representation='dirac')
cliff = gm.check_clifford_algebra()
check('Clifford algebra ({γμ,γν}=2gμν)', max(cliff.values()) < 1e-12)

g5 = gm.check_gamma5_properties()
check('(γ⁵)² = I', g5['square_is_I'])
check('{γ⁵,γμ} = 0', g5['anticommutes'])
check('(γ⁵)† = γ⁵', g5['hermitian'])

sigma01 = gm.get_sigma_munu(0, 1)
check('σμν shape (4,4)', sigma01.shape == (4, 4))
check('σνμ = -σμν', np.allclose(gm.get_sigma_munu(1, 0), -sigma01))

ds = DiracSpinor(gm)
p = np.array([0.0, 0.0, 1.0])
mass = 0.511
u = ds.u_spinor(p, mass, spin=1)
v = ds.v_spinor(p, mass, spin=1)
check('u_spinor shape (4,)', u.shape == (4,))
check('v_spinor shape (4,)', v.shape == (4,))

u_bar = ds.adjoint(u)
v_bar = ds.adjoint(v)
ubaru = float((u_bar @ u).real)
vbarv = float((v_bar @ v).real)
check('ūu = 2m', abs(ubaru - 2*mass) < 1e-10)
check('v̄v = -2m', abs(vbarv + 2*mass) < 1e-10)

E = np.sqrt(1.0 + mass**2)
p4 = np.array([E, 0.0, 0.0, 1.0])
sum_u_f = spin_sum_u(p4, mass, gm)
sum_u_ex = spin_sum_u_from_spinors(p, mass, gm)
check('Σ uū formula=explicit', np.max(np.abs(sum_u_f - sum_u_ex)) < 1e-10)

sum_v_f = spin_sum_v(p4, mass, gm)
sum_v_ex = spin_sum_v_from_spinors(p, mass, gm)
check('Σ vv̄ formula=explicit', np.max(np.abs(sum_v_f - sum_v_ex)) < 1e-10)

chk_u = dirac_equation_check(u, p4, mass, gm)
chk_v = dirac_equation_check(v, p4, mass, gm)
check('(p̸-m)u = 0', chk_u['is_u_solution'])
check('(p̸+m)v = 0', chk_v['is_v_solution'])

# Chiral representation
gm_c = GammaMatrices(representation='chiral')
g5c = gm_c.check_gamma5_properties()
check('Chiral (γ⁵)²=I', g5c['square_is_I'])

# ================================================================
# qed.py
# ================================================================
print("── qed.py ──")
from src.qft.qed import (
    ALPHA_QED, M_E, M_MU,
    mandelstam, mandelstam_from_cm,
    compton_cross_section, compton_amplitude_squared,
    pair_annihilation_cross_section, pair_annihilation_amplitude_squared,
    moller_cross_section, moller_amplitude_squared,
    phase_space_factor,
)

p1 = np.array([1.,0.,0.,1.])
p2 = np.array([1.,0.,0.,-1.])
p3 = np.array([1.,0.,1.,0.])
p4 = np.array([1.,0.,-1.,0.])
s, t, u = mandelstam(p1, p2, p3, p4)
check('Mandelstam s=4 (2→2 forward)', abs(s-4.)<1e-10)
check('Mandelstam t=-2', abs(t+2.)<1e-10)
check('Mandelstam u=-2', abs(u+2.)<1e-10)

s_cm = 100.0
s2, t2, u2 = mandelstam_from_cm(s_cm, np.pi/3)
check('CM Mandelstam s invariant', abs(s2-s_cm)<1e-10)
check('s+t+u=0 (m=0)', abs(s2+t2+u2)<1e-10)

cs = compton_cross_section(1.0, np.pi/2)
check('Compton dσ/dΩ > 0', cs > 0)

cs_low = compton_cross_section(1e-6, np.pi/2)
thomson = ALPHA_QED**2 / (2*M_E**2)  # low-energy → Thomson
check('Low-E Compton ~ Thomson', abs(cs_low/thomson - 1) < 0.01)

amp2_c = compton_amplitude_squared(10.0, -5.0)
check('Compton |M|² > 0', amp2_c > 0)

amp2_ann = pair_annihilation_amplitude_squared(50000.0, -12500.0, -12500.0)
check('Pair ann |M|² > 0', amp2_ann > 0)

ds_ann = pair_annihilation_cross_section(50000.0, np.pi/3)
check('e⁺e⁻→μ⁺μ⁻ dσ/dΩ > 0 (above thr)', ds_ann > 0)

amp2_mol = moller_amplitude_squared(100.0, -25.0, -25.0)
check('Møller |M|² > 0', amp2_mol > 0)

ds_mol = moller_cross_section(10000.0, np.pi/3)
check('Møller dσ/dΩ > 0', ds_mol > 0)

ps = phase_space_factor(100.0, 0.5, 0.5)
check('Phase space factor > 0', ps > 0)

# ================================================================
print()
print('='*55)
total = passed + failed
print(f'  {passed}/{total} passed' + ('  ✓ All clear!' if not failed else f'  ({failed} FAILED)'))
print('='*55)
sys.exit(1 if failed else 0)

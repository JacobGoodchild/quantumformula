"""Quick re-verification of every formula in formulas.txt (formula vs an independent numerical
computation at one parameter point each). Run from the repo root:  python3 code/verify_all.py
Each line prints the formula value, the independent value and their difference."""
import sys, os
here=os.path.dirname(os.path.abspath(__file__))
for d in ['lieb_metric','twoband_metric','diamond_metric','lieb3d_metric']:
    sys.path.insert(0,os.path.join(here,d))
import mpmath as mp
mp.mp.dps=30
def show(name,a,b):
    print(f"{name:58s} formula {mp.nstr(a,18):>24s}  numeric {mp.nstr(b,18):>24s}  diff {mp.nstr(a-b,2)}")

# F1(a,b): Lieb metric
from aniso import Ig
d1,d2=mp.mpf('0.37'),mp.mpf('0.81'); k1=(1-d1**2)/(1+d1**2); k2=(1-d2**2)/(1+d2**2); m=k1*k2
f1=(2+k1+k2)/(8*mp.sqrt((1+k1)*(1+k2)))*(2*mp.ellipk(m)-mp.ellipe(m))/mp.pi
show("F1(b) Lieb metric, (0.37,0.81)",f1,Ig(d1,d2))
# F1(d) component
from gxx_formula import gxx
from aniso3 import parts
show("F1(d) Lieb g_xx, (0.37,0.81)",gxx(d1,d2),parts(d1,d2)[0])
from lieb_general_formula import Mgen
from lieb4 import M4
show("F1(g) general Lieb metric, (1.37,0.61,2.2,0.93)",Mgen('1.37','0.61','2.2','0.93'),M4(*[mp.mpf(x) for x in ('1.37','0.61','2.2','0.93')]))
# F2 QWZ metric
from qwz_formula import Mformula
from qwz_fast import M as Mq
show("F2 QWZ metric, m = 0.7",Mformula('0.7'),Mq('0.7',30)); mp.mp.dps=30
# F3 gapped graphene
from honey_formula import Mhoney
from honey_trap import M as Mh
show("F3 gapped graphene metric, Delta = 1.3",Mhoney('1.3'),Mh('1.3',110,30)); mp.mp.dps=30
# F4 diamond (vs Joyce form; direct 3D check is slow - see hp_direct logs)
from diamond_explicit import M_explicit
from diamond import M_formula
show("F4 diamond metric, Delta = 1.3 (explicit vs Joyce form)",M_explicit('1.3')[0],M_formula('1.3')[0])
# F5 3D Lieb: delta -> 0 limit vs Watson/4
from lieb3d import M3, watson
show("F5 3D Lieb metric at delta=1e-6 vs Watson/4 (~1e-6 off)",M3('1e-6'),watson()/4)
# F6 QWZ <Omega^2>
from qwz_omega2_formula import Om2formula as O6
from qwz_omega2 import Om2 as O6n
show("F6 QWZ <Omega^2>, m = 0.7",O6('0.7'),O6n('0.7',30)); mp.mp.dps=30
# F7 gapped graphene <Omega^2>
from honey_omega2_formula import Om2formula as O7
from honey_omega2_pslq import Om2s
show("F7 graphene <Omega^2>, Delta^2 = 1.7",O7('1.7'),Om2s('1.7',120,30)); mp.mp.dps=30
# F8 diamond <|Omega|^2>: compare with stored 32-digit direct result at s = 9
from diamond_omega2_formula import Om2 as O8
show("F8 diamond <|Omega|^2>, s = 9 (vs stored direct value)",O8(9),mp.mpf('0.000177802179088358877213937101921428'))
# F9 Lieb trapping
from trap_formula import wB_formula
from trap import wB
show("F9 Lieb CTQW trapped fraction, (0.37,0.81)",wB_formula(d1,d2),wB(d1,d2))
from trap_general_formula import wBgen
from trap_general_fit import wB4
show("F9(general) trapped fraction, (1.37,0.61,2.2,0.93)",wBgen('1.37','0.61','2.2','0.93'),wB4(*[mp.mpf(x) for x in ('1.37','0.61','2.2','0.93')]))
# F10 strained graphene magic mass
from strained_hp import M as Ms
mp.mp.dps=30; t=['1.3','0.8','1']; T=sum(mp.mpf(x)**2 for x in t)
show("F10 strained graphene at magic mass (= 1/48)",mp.mpf(1)/48,Ms(t,3*T,90))

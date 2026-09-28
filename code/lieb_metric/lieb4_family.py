from lieb4 import M4
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp
mp.dps=60
def basis(a1,a2,b1,b2):
    Px=a1**2+a2**2; Py=b1**2+b2**2; B1=2*a1*a2; B2=2*b1*b2
    c1=(B2-Px-Py)/B1; c2=(-B2-Px-Py)/B1
    k2=2*(c1-c2)/((1-c1)*(-1-c2)); pre=sqrt((1-c1)*(-1-c2))
    return k2,pre
rows=[]
for a in [F(3,2),F(2),F(5,2),F(3),F(4),F(5),F(6),F(7,2),F(9,2)]:
    A=mpf(a.numerator)/a.denominator
    v=M4(A,mpf(1),mpf(1),mpf(1)); k2,pre=basis(A,1,1,1)
    rel=pslq([v,ellipk(k2)/pi/pre,ellipe(k2)/pi/pre],maxcoeff=10**9,maxsteps=10**6)
    al,be=F(-rel[1],rel[0]),F(-rel[2],rel[0]); rows.append((a,al,be))
    print(a,'k2=',mp.nstr(k2,12),'pre^2=',F(str(mp.nstr(pre**2,40))).limit_denominator(10**6),'alpha',al,'beta',be)
x=sp.Symbol('a')
for idx,name in ((1,'alpha'),(2,'beta')):
    pts=[(sp.Rational(r[0].numerator,r[0].denominator),sp.Rational(r[idx].numerator,r[idx].denominator)) for r in rows]
    for deg in range(1,5):
        try: rr=sp.polys.polyfuncs.rational_interpolate(pts[:2*deg+1],deg,X=x)
        except Exception: continue
        if all(rr.subs(x,p)==q for p,q in pts): print(name,'=',sp.factor(rr)); break

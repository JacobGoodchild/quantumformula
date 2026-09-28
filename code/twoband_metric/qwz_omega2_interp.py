from qwz_omega2 import Om2
from mpmath import mp, mpf, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp
xs=[];As=[];Bs=[]
for M in range(3,25):
    mm=mpf(M); v=Om2(mm,70); mp.dps=70
    k2=8*(mm**2-2)/mm**4; K=ellipk(k2);E=ellipe(k2)
    rel=pslq([v,K/pi,E/pi],maxcoeff=10**15,maxsteps=10**7)
    a=F(-rel[1],rel[0]); b=F(-rel[2],rel[0])
    xs.append(M);As.append(a);Bs.append(b)
x=sp.Symbol('m')
for name,ys in (('a',As),('b',Bs)):
    for deg in range(2,10):
        try:
            r=sp.polys.polyfuncs.rational_interpolate(list(zip(xs[:2*deg+2],[sp.Rational(y.numerator,y.denominator) for y in ys[:2*deg+2]])),deg,X=x)
        except Exception as e: continue
        ok=all(sp.Rational(y.numerator,y.denominator)==r.subs(x,xx) for xx,y in zip(xs,ys))
        if ok:
            print(name,'=',sp.factor(r)); break

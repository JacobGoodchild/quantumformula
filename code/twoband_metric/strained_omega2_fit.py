from strained_omega2 import Om2, Gw
from mpmath import mp, mpf, diff, pslq
from fractions import Fraction as F
import sympy as sp
W=sp.Symbol('w')
for a in [F(3,2),F(1,2),F(2)]:
    ts=[mpf(a.numerator)/a.denominator,mpf(1),mpf(1)]; T=a*a+2
    A=[];Bc=[]
    for D2 in [1,2,3,4,5,6]:
        mp.dps=45; v=Om2(ts,D2,150); mp.dps=45
        w=mpf(T.numerator)/T.denominator+D2
        g=Gw(ts,w); g1=diff(lambda x: Gw(ts,x), w)
        rel=pslq([v,g,g1],maxcoeff=10**10,maxsteps=10**6,tol=mpf(10)**-35)
        al=F(-rel[1],rel[0])*384; be=F(-rel[2],rel[0])*384
        wr=T+D2; A.append((wr,al)); Bc.append((wr,be))
    pa=sp.interpolate([(sp.Rational(x.numerator,x.denominator),sp.Rational(y.numerator,y.denominator)) for x,y in A],W)
    pb=sp.interpolate([(sp.Rational(x.numerator,x.denominator),sp.Rational(y.numerator,y.denominator)) for x,y in Bc],W)
    print('a=',a,' -384*coef(G) =',sp.expand(-pa),'   -384*coef(dG/dw) =',sp.expand(-pb))

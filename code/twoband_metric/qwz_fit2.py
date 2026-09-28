from qwz_fit import data
from mpmath import mp, mpf, pi, pslq
from fractions import Fraction as F
mp.dps=60
for q in [F(3,2),F(5,4),F(7,4),F(1,2),F(1,3),F(2,3),F(4,3)]:
    mm=mpf(q.numerator)/q.denominator
    v,K,E,P,k2,n=data(mm)
    rel=pslq([v,1,K/pi,E/pi,P/pi],maxcoeff=10**7,maxsteps=10**6)
    out=[F(-rel[i],rel[0]) for i in (1,2,3,4)] if rel and rel[0] else rel
    print(q,'v=',mp.nstr(v,20),'k2=',F(str(k2)).limit_denominator(10**6),'n=',F(str(n)).limit_denominator(10**6),'c,a,e,b=',out)

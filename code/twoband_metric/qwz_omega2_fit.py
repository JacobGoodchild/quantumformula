from qwz_omega2 import Om2
from mpmath import mp, mpf, pi, ellipk, ellipe, ellippi, pslq
from fractions import Fraction as F
for q in [F(3),F(4),F(5),F(5,2),F(6),F(3,2),F(1,2),F(8)]:
    mm=mpf(q.numerator)/q.denominator
    v=Om2(mm,60); mp.dps=60
    k2=8*(mm**2-2)/mm**4; n=4*(mm-1)/mm**2
    K=ellipk(k2);E=ellipe(k2);P=ellippi(n,k2)
    rel=pslq([v,1,K/pi,E/pi,P/pi],maxcoeff=10**9,maxsteps=10**6)
    print(q, [F(-rel[i],rel[0]) for i in (1,2,3,4)] if rel and rel[0] else rel)

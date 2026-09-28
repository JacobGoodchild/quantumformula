from aniso3 import parts
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import math
mp.dps=60
def issq(x): return math.isqrt(x.numerator)**2==x.numerator and math.isqrt(x.denominator)**2==x.denominator
for k1 in [F(1,3),F(2,3),F(1,4),F(3,4),F(1,5),F(2,5)]:
    for k2 in [F(1,3),F(1,7),F(2,7),F(3,8),F(5,9)]:
        if k1==k2 or issq((1+k1)*(1+k2)): continue
        K1=mpf(k1.numerator)/k1.denominator; K2=mpf(k2.numerator)/k2.denominator
        d1=sqrt((1-K1)/(1+K1)); d2=sqrt((1-K2)/(1+K2))
        gx,gy=parts(d1,d2); m=K1*K2; K=ellipk(m);E=ellipe(m); r=sqrt((1+K1)*(1+K2))
        rel=pslq([gx,K*r/pi,E*r/pi],maxcoeff=10**9,maxsteps=10**6)
        a,b=F(-rel[1],rel[0]),F(-rel[2],rel[0])
        print(k1,k2,'alpha',a,'beta*(1-m)',b*(1-k1*k2))

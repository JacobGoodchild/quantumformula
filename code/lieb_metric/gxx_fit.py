from aniso3 import parts
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
mp.dps=60
out=[]
for k1,k2 in [(F(1,2),F(1,3)),(F(1,2),F(1,5)),(F(2,3),F(1,4)),(F(1,3),F(1,7)),(F(3,4),F(1,2)),(F(1,5),F(2,7)),(F(3,5),F(1,6))]:
    K1=mpf(k1.numerator)/k1.denominator; K2=mpf(k2.numerator)/k2.denominator
    d1=sqrt((1-K1)/(1+K1)); d2=sqrt((1-K2)/(1+K2))
    gx,gy=parts(d1,d2)
    m=K1*K2; K=ellipk(m);E=ellipe(m); r=sqrt((1+K1)*(1+K2))
    rel=pslq([gx,K/pi/r,E/pi/r,K*r/pi,E*r/pi],maxcoeff=10**9,maxsteps=10**6)
    co=[F(-rel[i],rel[0]) for i in range(1,5)] if rel and rel[0] else rel
    print(k1,k2,co)

from lieb4 import M4
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, ellippi, pslq
from fractions import Fraction as F
mp.dps=60
for h in [(2,1,1,F(1,3)),(F(3,2),F(1,2),F(2,3),F(1,3)),(3,1,1,1),(2,1,1,1)]:
    a1,a2,b1,b2=[mpf(x.numerator)/x.denominator if isinstance(x,F) else mpf(x) for x in h]
    v=M4(a1,a2,b1,b2)
    Px=a1**2+a2**2; Py=b1**2+b2**2; B1=2*a1*a2; B2=2*b1*b2
    c1=(B2-Px-Py)/B1; c2=(-B2-Px-Py)/B1
    k2=2*(c1-c2)/((1-c1)*(-1-c2)); n=2/(1-c1); pre=sqrt((1-c1)*(-1-c2))
    K=ellipk(k2);E=ellipe(k2);P=ellippi(n,k2)
    print([str(x) for x in h], pslq([v,1,K/pi/pre,E/pi/pre,P/pi/pre,K/pi*pre,E/pi*pre],maxcoeff=10**8,maxsteps=10**6))

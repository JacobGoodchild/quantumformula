from aniso import Ig
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, identify
mp.dps=40
from fractions import Fraction as F
for d1,d2 in [(F(1,2),F(1,3)),(F(1,2),F(1,2)),(F(1,2),F(1,5)),(F(1,3),F(1,5)),(F(2,3),F(1,4)),(F(1,2),F(1,1)), (F(1,2),F(0,1))]:
    D1=mpf(d1.numerator)/d1.denominator; D2=mpf(d2.numerator)/d2.denominator
    v=Ig(D1,D2)
    A0=2*(1+D1**2)+2*(1+D2**2); B1=2*(1-D1**2); B2=2*(1-D2**2)
    s=sqrt(A0**2-(B1-B2)**2); m=4*B1*B2/s**2
    c=v*pi/(2*ellipk(m)-ellipe(m))
    print(d1,d2, c*s, identify(c*s), identify(c), s**2)

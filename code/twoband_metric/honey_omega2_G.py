from honey_omega2_pslq import Om2s
from mpmath import mp, mpf, pi, sqrt, ellipk, diff, pslq
from fractions import Fraction as F
def Gs(s):   # G_tri(3+s) as function of s = Delta^2
    P=sqrt((1+s)*(9+s)); k2=mpf(1)/2+(3-6*s-s*s)/(2*(1+s)*P)
    return 2*ellipk(k2)/(pi*(1+s)**mpf(0.75)*(9+s)**mpf(0.25))
res={}
for s,N in [(2,130),(3,120),(5,110),(7,110),(1,150),(4,115),(6,110),(9,100),(11,100)]:
    v=Om2s(s,N,70); mp.dps=70; s=mpf(s)
    g0=Gs(s); g1=diff(Gs,s); g2=diff(Gs,s,2)
    rel=pslq([v,1,g0,g1],maxcoeff=10**12,maxsteps=10**7)
    print(int(s), [F(-rel[i],rel[0]) for i in range(1,4)] if rel and rel[0] else rel)

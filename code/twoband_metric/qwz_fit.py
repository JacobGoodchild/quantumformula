from qwz_fast import M
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, ellippi, pslq
from fractions import Fraction as F
mp.dps=60
def data(mm):
    v=M(mm,60)
    r1=-(mm*mm-2*mm+2)/(2*(mm-1)); r2=-(mm*mm+2*mm+2)/(2*(mm+1))
    e1,e2,e3,e4=sorted([mpf(1),mpf(-1),r1,r2],reverse=True)
    k2=(e1-e2)*(e3-e4)/((e1-e3)*(e2-e4)); n=(e1-e2)/(e1-e3)
    K=ellipk(k2);E=ellipe(k2);P=ellippi(n,k2)
    return v,K,E,P,k2,n
for q in [F(5,2),F(3),F(7,2),F(4),F(5),F(6),F(9,4),F(8)]:
    mm=mpf(q.numerator)/q.denominator
    v,K,E,P,k2,n=data(mm)
    rel=pslq([v-mpf(1)/8,K/pi,E/pi,P/pi],maxcoeff=10**7,maxsteps=10**6)
    out=None
    if rel and rel[0]!=0:
        out=[F(-rel[i],rel[0]) for i in (1,2,3)]
    print(q, 'k2=',F(str(k2)).limit_denominator(10**6),'n=',F(str(n)).limit_denominator(10**6), 'alpha,eps,beta=',out)

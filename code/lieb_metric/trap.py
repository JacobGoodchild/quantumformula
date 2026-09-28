# Flat-band weight on an edge site B (x-bond site) of the anisotropic staggered Lieb lattice:
# w_B = < |g|^2 / (|f|^2 + |g|^2) >  = probability a CTQW started on B stays trapped forever;
# long-time-average return probability = w_B^2.
from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, ellippi, pslq
def wB(d1,d2):
    p1=1+d1*d1;q1=1-d1*d1;p2=1+d2*d2;q2=1-d2*d2
    def fx(x):
        a=2*p1+2*q1*cos(x); return a/sqrt((a+2*p2)**2-4*q2**2)
    return 1-quad(fx,[0,pi])/pi
if __name__=='__main__':
    from fractions import Fraction as F
    mp.dps=50
    for d1,d2 in [(mpf(1)/2,mpf(1)/3),(mpf(1)/3,mpf(1)/2),(mpf(1)/2,mpf(1)/5),(mpf(2)/3,mpf(1)/4)]:
        w=wB(d1,d2)
        p1=1+d1*d1;q1=1-d1*d1;p2=1+d2*d2;q2=1-d2*d2
        c1=-(p1+p2-q2)/q1; c2=-(p1+p2+q2)/q1
        e1,e2,e3,e4=1,-1,c1,c2
        k2=(e1-e2)*(e3-e4)/((e1-e3)*(e2-e4)); n=(e1-e2)/(e1-e3); pre=sqrt((e1-e3)*(e2-e4))
        K=ellipk(k2);P=ellippi(n,k2);E=ellipe(k2)
        print(mp.nstr(d1,5),mp.nstr(d2,5),'w=',mp.nstr(w,25),'k2=',mp.nstr(k2,15),'n=',mp.nstr(n,15),'pre^2=',mp.nstr(pre**2,10))
        print('  ',pslq([w,1,K/pi/pre,P/pi/pre,K/pi,P/pi],maxcoeff=10**7,maxsteps=10**6))

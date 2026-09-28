from mpmath import mp, mpf, mpc, exp, pi, sqrt, quad, cos, diff, pslq
from fractions import Fraction as F
I=mpc(0,1)
def om2(k1,k2,t,D):
    ph=(k1+k2)/3; e0=exp(I*ph); e1=e0*exp(-I*k1); e2=e0*exp(-I*k2)
    f=t[0]*e0+t[1]*e1+t[2]*e2
    f1=I/3*f - I*t[1]*e1; f2=I/3*f - I*t[2]*e2
    fx=f1+f2/2; fy=sqrt(3)/2*f2
    J=(fx.conjugate()*fy).imag
    return D*D*J*J/(4*(abs(f)**2+D*D)**3)
def Om2(t,D2,N):
    t=[mpf(x) for x in t]; D=sqrt(mpf(D2)); h=2*pi/N; s=0
    for i in range(N):
        for j in range(N): s+=om2(i*h,j*h,t,D)
    return s/N**2
def Gw(t,w):   # <1/(|f|^2 - T + w)>  as function of w
    t1,t2,t3=[mpf(x) for x in t]; T=t1**2+t2**2+t3**2
    def fk(k1):
        A=w+2*t1*t2*cos(k1); R2=4*t3**2*(t1**2+t2**2+2*t1*t2*cos(k1))
        return 1/sqrt(A*A-R2)
    return quad(fk,[0,pi])/pi
if __name__=='__main__':
    mp.dps=40
    for t,D2 in [((F(3,2),F(1),F(1)),2),((F(3,2),F(1),F(1)),3),((F(1,2),F(1),F(1)),2),((F(3,2),F(1,2),F(1)),2)]:
        ts=[mpf(x.numerator)/x.denominator for x in t]; T=sum(x*x for x in ts); w=T+D2
        v=Om2(ts,D2,140); mp.dps=40
        g=Gw(ts,w); g1=diff(lambda x: Gw(ts,x), w); g2=diff(lambda x: Gw(ts,x), w, 2)
        rel=pslq([v,1,g,g1],maxcoeff=10**8,maxsteps=10**6,tol=mpf(10)**-30)
        rel2=pslq([v,1,g,g1,g2],maxcoeff=10**8,maxsteps=10**6,tol=mpf(10)**-30)
        print([str(x) for x in t],D2,'T=',T, rel, rel2)

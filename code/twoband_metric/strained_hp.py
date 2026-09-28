# 50-digit spectral check of the universal value at the magic mass for arbitrary hoppings
from mpmath import mp, mpf, mpc, exp, pi, sqrt
I=mpc(0,1)
def trg(k1,k2,t,D):
    ph=(k1+k2)/3; e0=exp(I*ph); e1=e0*exp(-I*k1); e2=e0*exp(-I*k2)
    f=t[0]*e0+t[1]*e1+t[2]*e2
    f1=I/3*f - I*t[1]*e1; f2=I/3*f - I*t[2]*e2
    fx=f1+f2/2; fy=sqrt(3)/2*f2
    d2=abs(f)**2+D*D; s=0
    for g in (fx,fy):
        s+=(abs(g)**2*d2-(mp.re(f.conjugate()*g))**2)/(4*d2*d2)
    return s
def M(t,D2,N):
    t=[mpf(x) for x in t]; D=sqrt(mpf(D2)); h=2*pi/N; s=0
    for i in range(N):
        for j in range(N): s+=trg(i*h,j*h,t,D)
    return s/N**2
if __name__=='__main__':
    mp.dps=50
    t=['1.3','0.8','1']; T=sum(mpf(x)**2 for x in t)
    for N in (80,110): print(N, M(t,3*T,N)-mpf(1)/48)

# <Omega^2> for gapped graphene (Bravais a=1, physical positions) by spectral trapezoid rule
from mpmath import mp, mpf, mpc, exp, pi, sqrt
import sys
I=mpc(0,1)
def om2(k1,k2,D):
    ph=(k1+k2)/3
    f = exp(I*ph)*(1+exp(-I*k1)+exp(-I*k2))
    f1= I/3*f - I*exp(I*ph)*exp(-I*k1); f2= I/3*f - I*exp(I*ph)*exp(-I*k2)
    fx=f1+f2/2; fy=sqrt(3)/2*f2
    J=(fx.conjugate()*fy).imag
    return D*D*J*J/(4*(abs(f)**2+D*D)**3)
def Om2(D,N,dps):
    mp.dps=dps; D=mpf(D); h=2*pi/N; s=0
    for i in range(N):
        for j in range(N): s+=om2(i*h,j*h,D)
    return s/N**2
if __name__=='__main__':
    for N in map(int,sys.argv[2:]): print(N, Om2(sys.argv[1],N,50))

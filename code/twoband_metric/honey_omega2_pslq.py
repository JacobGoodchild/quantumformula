from honey_omega2 import om2
from mpmath import mp, mpf, pi, sqrt, ellipk, ellipe, pslq
def Om2s(s,N,dps):
    mp.dps=dps; D=sqrt(mpf(s)); h=2*pi/N; t=0
    for i in range(N):
        for j in range(N): t+=om2(i*h,j*h,D)
    return t/N**2
if __name__=='__main__':
    for s,N in [(7,110),(3,120),(15,100)]:
        v=Om2s(s,N,60); mp.dps=60; s=mpf(s)
        P=sqrt((1+s)*(9+s)); k2=mpf(1)/2+(3-6*s-s*s)/(2*(1+s)*P)
        K=ellipk(k2);E=ellipe(k2); r=(1+s)**mpf(0.25); q=(9+s)**mpf(0.25)
        basis=[v,1,K/pi/(r**3*q),E/pi/(r**3*q),K/pi*q**3/r**3, E/pi*q**3/r**3, K/pi/(r*q**3),E/pi/(r*q**3), E/pi*r/q**3, E/pi*q/r]
        print(int(s), v, pslq(basis,maxcoeff=10**8,maxsteps=10**6))

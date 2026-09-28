# independent: sum over l=1 excited states, alpha = 2 sum |<n|z|0>|^2/(E_n-E0)
from mpmath import mp, mpf, pi, sin, cos, sqrt, quad, findroot, tan
mp.dps=30
def j1(x): return sin(x)/x**2-cos(x)/x
def zeros(N):
    out=[]
    for n in range(1,N+1):
        x0=(n+mpf(1)/2)*pi-mpf(1)/((n+mpf(1)/2)*pi)
        out.append(findroot(lambda x: tan(x)-x, x0))
    return out
norm0=sqrt(2*pi)   # int |sin(pi r)/r|^2 d^3r = 2 pi
tot=0; terms=[]
for k in zeros(400):
    # l=1 state: j1(k r) Y10 ; radial norm int_0^1 j1(kr)^2 r^2 dr
    nrm=quad(lambda r: (j1(k*r))**2*r**2,[0,1])
    # <n|z|0> = (1/sqrt3) * int j1(kr) r (sin(pi r)/r) r^2 dr * sqrt(4pi)/ ... do it with Y00=1/sqrt(4pi):
    # psi0 = sin(pi r)/r /sqrt(2pi) = R0(r) Y00 with R0 = sin(pi r)/r*sqrt(2/1)... handle via radial norms:
    R0n=quad(lambda r: (sin(pi*r)/r)**2*r**2,[0,1])
    rad=quad(lambda r: j1(k*r)*r*(sin(pi*r)/r)*r**2,[0,1])
    me2=rad**2/(nrm*R0n)/3          # |<n|z|0>|^2 with <Y10|cos|Y00> = 1/sqrt3
    t=2*me2/((k**2-pi**2)/2); tot+=t; terms.append(t)
print('partial sum 400 states:',tot,' last term',terms[-1])

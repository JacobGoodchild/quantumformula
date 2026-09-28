import sympy as sp
from mpmath import mp, mpf, pi as mpi, findroot, tan, sin, cos, nsum, inf
r,k=sp.symbols('r k',positive=True); pi=sp.pi
j1=sp.sin(k*r)/(k*r)**2-sp.cos(k*r)/(k*r)
rad=sp.simplify(sp.integrate(j1*sp.sin(pi*r)*r**2,(r,0,1)))
nrm=sp.simplify(sp.integrate(sp.expand_trig(j1**2*r**2),(r,0,1)))
print('rad=',rad); print('nrm=',nrm)
radf=sp.lambdify(k,rad,'mpmath'); nrmf=sp.lambdify(k,nrm,'mpmath')
mp.dps=30
R0n=mpf(1)/2
def term(n):
    x0=(n+mpf(1)/2)*mpi-1/((n+mpf(1)/2)*mpi)
    kk=findroot(lambda x: tan(x)-x,x0)
    me2=radf(kk)**2/(nrmf(kk)*R0n)/3
    return 2*me2/((kk**2-mpi**2)/2)
s=0
for N in [10,100,1000,4000]:
    pass
parts=[term(n) for n in range(1,3001)]
import itertools
acc=list(itertools.accumulate(parts))
for N in [10,100,1000,3000]: print(N, acc[N-1])
print('exact', (3+4*mpi**2)/(12*mpi**4))
print('tail estimate ratio', parts[-1]*3000**5, parts[1499]*1500**5)

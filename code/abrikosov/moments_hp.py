from gamma6 import moments
from mpmath import mp, mpf, mpc, sqrt, pi, gamma, log, pslq
import sys
mp.dps=45
for name,tau,G,extra in [('square',mpc(0,1),gamma(mpf(1)/4),log(3)),('triangular',mpc(mpf(1)/2,sqrt(3)/2),gamma(mpf(1)/3),log(3))]:
    for M in (70,90):
        r=moments(tau,[2,3,4,5],M=M); print(name,M,[mp.nstr(x,40) for x in r])
    for p,val in zip([2,3,4,5],r):
        print('  p=',p, pslq([log(val),log(pi),log(G),log(2),log(3),log(5),log(7)],maxcoeff=200,maxsteps=10**6))

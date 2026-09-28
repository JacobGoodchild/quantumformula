from betan import *
from mpmath import findroot
mp.dps=20
import itertools
def b(tau,n):
    return betan(sums(tau,N=12,J=2*n),n)
best={}
for n in (0,1,2):
    res=[]
    for x in [i/20 for i in range(0,11)]:
        for yv in [0.5+j*0.05 for j in range(0,40)]:
            t=mpc(x,yv)
            if abs(t)<1: continue
            res.append((b(t,n),x,yv))
    res.sort(); print('n=',n,'best grid:',[ (mp.nstr(r[0],10),r[1],round(r[2],3)) for r in res[:3]])

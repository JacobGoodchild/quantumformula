from gamma6 import moments
from mpmath import mp, mpf, mpc, sqrt, pi, gamma, log, pslq, cbrt, root
mp.dps=70
rs=moments(mpc(0,1),[2,3,4,5],M=120)
rt=moments(mpc(mpf(1)/2,sqrt(3)/2),[2,3,4,5],M=120)
G4=gamma(mpf(1)/4); G3=gamma(mpf(1)/3)
Os=G4**2/(2*pi**mpf(1.5))
v=rs[2]
print('square p4',v)
print(pslq([v,Os**3,Os**3*sqrt(2),Os**3*root(2,4),Os**3*root(8,4),1/Os,pi*Os**5],maxcoeff=10**6,maxsteps=10**6))
print(pslq([log(rs[3]),log(pi),log(G4),log(2),log(3),log(5),log(7)],maxcoeff=500,maxsteps=10**6))
print('tri p5',rt[3], pslq([log(rt[3]),log(pi),log(G3),log(2),log(3),log(5),log(7)],maxcoeff=500,maxsteps=10**6))
print('tri p4', pslq([log(rt[2]),log(pi),log(G3),log(2),log(3),log(5),log(7)],maxcoeff=500,maxsteps=10**6))
print('tri p3', pslq([log(rt[1]),log(pi),log(G3),log(2),log(3),log(5),log(7)],maxcoeff=500,maxsteps=10**6))
print('sq p3', pslq([log(rs[1]),log(pi),log(G4),log(2),log(3),log(5),log(7)],maxcoeff=500,maxsteps=10**6))

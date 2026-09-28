from mpmath import mp, mpf, ellipk, quad, pi, sqrt, zeta, log, pslq, gamma
mp.dps=40
from mpmath import agm
Kp=lambda k: pi/2/agm(1,k)
K=lambda k: ellipk(k*k)
for name,f in [('Kp^2',lambda k:Kp(k)**2),('K^2',lambda k:K(k)**2),('KKp',lambda k:K(k)*Kp(k)),('K^2/kp',lambda k:K(k)**2/sqrt(1-k*k))]:
    v=quad(f,[0,mpf(1)/2,1]); print(v)
    print(name,v, pslq([v,pi**3,pi*mp.catalan,zeta(3),pi**2*log(2),gamma(mpf(1)/4)**8/pi**2, pi*log(2)**2],maxcoeff=10**4,maxsteps=10**5))

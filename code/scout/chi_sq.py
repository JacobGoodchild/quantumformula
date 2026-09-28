from mpmath import mp, mpf, ellipk, quad, pi, sqrt, zeta, catalan, log, pslq, identify
mp.dps=40
K=lambda k: ellipk(k*k)
I = quad(lambda k: K(k)**2/sqrt(1-k*k), [0,mpf(1)/2,mpf(9)/10,1])
chi = 2/pi**3*I
print('I=',I,'chi=',chi)
basis=[I, zeta(3), pi**3, pi*catalan, pi**2*log(2), pi, log(2)*pi, catalan, pi**3*log(2)]
print(pslq(basis,maxcoeff=10**5,maxsteps=10**6))
from mpmath import gamma
g=gamma(mpf(1)/4); Kl=ellipk(mpf(1)/2)
b=[I, Kl**2*pi, Kl**4/pi, pi**3, pi**5/Kl**2, pi*log(2)*Kl**2, zeta(3), pi*catalan]
print(pslq(b,maxcoeff=10**5,maxsteps=10**6))
print(identify(I,['pi','gamma(1/4)','log(2)']))

from mpmath import mp, mpf, sqrt, ellipe, ellipk, quad, pi, cos, gamma, identify, pslq
mp.dps=50
def inner(a2):  # (1/2pi) int sqrt(a2+cos^2 x) dx
    return 2/pi*sqrt(1+a2)*ellipe(1/(1+a2))
s = quad(lambda y: inner(cos(y)**2), [0, pi/2])*2/pi
print(s)
g=gamma(mpf(1)/4)
K=ellipk(mpf(1)/2); E=ellipe(mpf(1)/2)
print(pslq([s, K/pi**2, E/pi**2, 1/(pi**2*K), 1/pi, 1/K, K/pi, E/pi, pi/K],maxcoeff=10**6,maxsteps=10**6))
print(identify(s,['pi','gamma(1/4)','sqrt(2)']))
g4=g**4
cands={'g4/pi^3':g4/pi**3,'pi/g4':pi/g4*pi,'1/pi':1/pi,'1':mpf(1),'g4/pi^3s2':g4/pi**3/sqrt(2),'pi2/g4 s2':pi**2/g4*sqrt(2), 's2/pi':sqrt(2)/pi,'s2':sqrt(2)}
names=list(cands); vals=[cands[n] for n in names]
r=pslq([s]+vals,maxcoeff=10**5,maxsteps=10**6); print(names, r)

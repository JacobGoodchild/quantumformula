from hp import *
mp.dps=60
for d in [mpf(1)/3, mpf(1)/5, mpf(2)/3, mpf(3)/4, mpf(1)/7]:
    v=Ig(d); k=(1-d*d)/(1+d*d); K=ellipk(k*k); E=ellipe(k*k)
    print(d, pslq([v,1,K/pi,E/pi],maxcoeff=10**8,maxsteps=10**6), [ (v*4*pi - (2*K-E)) ])

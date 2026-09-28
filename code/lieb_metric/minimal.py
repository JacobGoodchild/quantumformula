from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, findroot
mp.dps=30
def minimal(d):
    k=(1-d*d)/(1+d*d); kp=sqrt(1-k*k); K=ellipk(k*k); E=ellipe(k*k)
    return ((2*K-E)*E-kp**2*K**2)/(4*pi*E), mpf(1)/2-kp*K/(2*E)
if __name__=='__main__':
    import numpy as np
    from check_numeric2 import integ
    from scipy.optimize import minimize_scalar
    for d in [0.5,0.25]:
        tm,sstar=minimal(mpf(d))
        r=minimize_scalar(lambda s: integ(d,d,s,N=120),bounds=(0.0,0.5),method='bounded',options={'xatol':1e-7})
        print(d,'formula min',tm,'s*',sstar,' numeric min',r.fun,'s*',r.x)

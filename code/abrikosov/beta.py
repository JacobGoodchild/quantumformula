# Abrikosov parameter beta_A = <|psi|^4>/<|psi|^2>^2 for lowest-Landau-level vortex lattices.
# Standard formula: beta = sum_{k in reciprocal lattice} exp(-|k|^2 l^2 /2) with cell area 2 pi l^2
# (Kleiner-Roth-Autler). Write lattice by tau = x + i y (Im tau > 0): beta = sqrt(y) * sum |theta-type| ...
from mpmath import mp, mpf, exp, pi, sqrt, gamma, pslq, jtheta, identify
mp.dps=40
def beta(tau):
    # cell with sides 1 and tau, area y -> scale so area = 2 pi l^2 with l=1: side a with a^2 y = 2pi
    x,y=tau.real,tau.imag
    s=0
    N=30
    for m in range(-N,N+1):
        for n in range(-N,N+1):
            # reciprocal lattice vector norm^2 for lattice {a(m + n tau)}: |k|^2 = (2pi/(a^2 y))^2 |m tau - n|^2 * a^2 ... use formula:
            # sum_k exp(-|k|^2/2) with |k|^2 = (2 pi / y) |n - m tau|^2 / (a^2/ (2pi) ...)
            s+=exp(-pi*abs(n-m*tau)**2/y)
    return s
from mpmath import mpc
bs=beta(mpc(0,1)); bt=beta(mpc(mpf(1)/2,sqrt(3)/2))
print('square',bs,' triangular',bt)
g4=gamma(mpf(1)/4); g3=gamma(mpf(1)/3)
print(identify(bs,['pi','gamma(1/4)','sqrt(2)']), pslq([bs**2, g4**4/pi**3, g4**8/pi**6, 1],maxcoeff=10**6,maxsteps=10**6))
print(pslq([bt**2, g3**6/pi**4, g3**12/pi**8, g3**18/pi**12, 1],maxcoeff=10**6,maxsteps=10**6))
print(pslq([bt**4, g3**12/pi**8, g3**24/pi**16, 1],maxcoeff=10**6,maxsteps=10**6))

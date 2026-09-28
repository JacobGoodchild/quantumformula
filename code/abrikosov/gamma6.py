# <|psi|^{2p}>/<|psi|^2>^p for the LLL Abrikosov lattice, computed in real space from the theta-function
# form psi(x,y) = sum_n exp(i pi tau n^2 + 2 pi i n z/a)... use standard: psi = exp(-pi y^2/(a^2 Im tau)) theta
from mpmath import mp, mpf, mpc, exp, pi, sqrt, jtheta, quad, gamma, log, pslq
mp.dps=25
def moments(tau,p_list,M=60):
    # lattice spanned by 1 and tau (area Im tau); LLL function with one zero per cell:
    # psi(z) = exp(pi (z^2 - |z|^2)/(2 Im tau)) * theta_1(pi z | tau) up to shift; |psi|^2 doubly periodic.
    q=exp(1j*pi*tau); y0=tau.imag
    def rho(u,v):   # point z = u + v tau, u,v in [0,1)
        z=u+v*tau
        return abs(jtheta(1,pi*z,q))**2*exp(-pi*(2*z.imag**2)/y0)
    # midpoint rule on the torus (periodic smooth -> exponential convergence)
    vals=[rho(mpf(i)/M,mpf(j)/M) for i in range(M) for j in range(M)]
    m1=sum(vals)/len(vals)
    return [ (sum(v**p for v in vals)/len(vals))/m1**p for p in p_list]
for name,tau in [('square',mpc(0,1)),('triangular',mpc(mpf(1)/2,sqrt(3)/2))]:
    r=moments(tau,[2,3],M=50); print(name,'beta',r[0],'gamma6',r[1])

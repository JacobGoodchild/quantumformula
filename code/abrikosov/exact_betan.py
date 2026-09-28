# Exact higher-Landau-level Abrikosov parameters on square and triangular lattices.
from betan import sums
from mpmath import mp, mpf, mpc, pi, sqrt, gamma, pslq, jtheta, exp, diff
import sympy as sp
mp.dps=80
X,P,OM=sp.symbols('x P Om')
out={}
for name,tau,Om in [('square',mpc(0,1),gamma(mpf(1)/4)**2/(2*pi**mpf(1.5))),
                    ('triangular',mpc(mpf(1)/2,sqrt(3)/2),3*gamma(mpf(1)/3)**3/(2**(mpf(7)/3)*pi**2))]:
    S=sums(tau,N=45,J=8); Pv=pi**2*Om**4
    polys=[]
    for j in range(0,9):
        r=S[j]/S[0]; B=[Pv**b for b in range(0,j//2+1)]
        rel=pslq([r]+B,maxcoeff=10**12,maxsteps=10**7)
        co=[sp.Rational(-c,rel[0]) for c in rel[1:]]
        polys.append(sum(c*P**b for b,c in enumerate(co)))
    print(name,'S_j/S_0 =',polys)
    for n in range(0,5):
        L2=sp.Poly(sp.expand(sp.assoc_laguerre(n,0,X)**2),X)
        expr=sp.expand(sum(c*polys[k] for (k,),c in zip(L2.monoms(),L2.coeffs())))
        val=sum(mpf(sp.Rational(c).p)/sp.Rational(c).q*Pv**k for (k,),c in zip(sp.Poly(expr,P).monoms(),sp.Poly(expr,P).coeffs()))*Om
        # direct
        from betan import betan
        print('  beta_%d = Om*('%n, expr, ')  =', mp.nstr(val,30), ' direct-diff', mp.nstr(val-betan(S,n),3))

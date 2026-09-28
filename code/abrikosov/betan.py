from mpmath import mp, mpf, mpc, exp, pi, sqrt, gamma, pslq, laguerre
mp.dps=50
def sums(tau,N=40,J=6):
    x0,y=tau.real,tau.imag
    S=[mpf(0)]*(J+1)
    for m in range(-N,N+1):
        for n in range(-N,N+1):
            x=pi*abs(n-m*tau)**2/y; e=exp(-x); p=mpf(1)
            for j in range(J+1):
                S[j]+=p*e; p*=x
    return S
def betan(S,n):
    # sum e^{-x} L_n(x)^2 ; expand L_n(x)^2 as polynomial in x
    from mpmath import taylor
    import sympy as sp
    X=sp.Symbol('x'); poly=sp.Poly(sp.expand(sp.assoc_laguerre(n,0,X)**2),X)
    return sum(mpf(sp.Rational(c).p)/mpf(sp.Rational(c).q)*S[k] for (k,),c in zip(poly.monoms(),poly.coeffs()))
if __name__=='__main__':
    for name,tau,Om in [('square',mpc(0,1),gamma(mpf(1)/4)**2/(2*pi**mpf(1.5))),('triangular',mpc(mpf(1)/2,sqrt(3)/2),3*gamma(mpf(1)/3)**3/(2**(mpf(7)/3)*pi**2))]:
        S=sums(tau,J=8)
        print(name,'S0',S[0],'check closed',S[0]-Om, 'S1/S0',S[1]/S[0])
        for n in range(0,4):
            b=betan(S,n)
            rel=pslq([b,Om,Om/pi,Om/pi**2,Om**3,Om**3/pi, Om**5, 1/Om, 1/(pi*Om)],maxcoeff=10**6,maxsteps=10**6)
            print('  n=',n,'beta_n=',mp.nstr(b,25),rel)

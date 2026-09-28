from mpmath import mp, mpf, mpc, exp, pi, sqrt, gamma, pslq, diff
mp.dps=50
def B(x,y,sign,N=30):
    tau=mpc(x,y); s=0
    for m in range(-N,N+1):
        for n in range(-N,N+1):
            s+=(sign**(m+n))*exp(-pi*abs(n-m*tau)**2/y)
    return s
Om=gamma(mpf(1)/4)**2/(2*pi**mpf(1.5))
b_xx=diff(lambda x: B(x,1,1),0,2); b_yy=diff(lambda y: B(0,y,1),1,2)
c0=B(0,1,-1); c_xx=diff(lambda x: B(x,1,-1),0,2); c_yy=diff(lambda y: B(0,y,-1),1,2)
c_x=diff(lambda x: B(x,1,-1),0,1); c_y=diff(lambda y: B(0,y,-1),1,1)
print('B12(i)=',c0,' Om/sqrt2=',Om/sqrt(2)); print('grad B12', mp.nstr(c_x,3), mp.nstr(c_y,3))
print('b_xx',b_xx,' b_yy',b_yy); print('c_xx',c_xx,' c_yy',c_yy)
bas=lambda v: pslq([v,Om,sqrt(2)*Om,pi*Om**3,sqrt(2)*pi*Om**3,pi**2*Om**5,sqrt(2)*pi**2*Om**5],maxcoeff=10**8,maxsteps=10**6)
print('c_xx rel',bas(c_xx)); print('c_yy rel',bas(c_yy))
a2=-b_xx/c_xx; a3=-b_yy/c_yy
print('alpha2 (rhombic->square) =',a2); print('alpha3 (square->rectangular) =',a3)

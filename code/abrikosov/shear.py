# Hessian of the LLL Abrikosov parameter beta(tau) at the triangular and square points.
from mpmath import mp, mpf, mpc, exp, pi, sqrt, gamma, pslq, diff
mp.dps=50
def beta(x,y,N=30):
    tau=mpc(x,y); s=0
    for m in range(-N,N+1):
        for n in range(-N,N+1):
            s+=exp(-pi*abs(n-m*tau)**2/y)
    return s
for name,x0,y0,Om in [('triangular',mpf(1)/2,sqrt(3)/2,3*gamma(mpf(1)/3)**3/(2**(mpf(7)/3)*pi**2)),
                      ('square',mpf(0),mpf(1),gamma(mpf(1)/4)**2/(2*pi**mpf(1.5)))]:
    bxx=diff(lambda x: beta(x,y0),x0,2); byy=diff(lambda y: beta(x0,y),y0,2); bxy=diff(lambda x,y: beta(x,y),(x0,y0),(1,1))
    bx=diff(lambda x: beta(x,y0),x0,1); by=diff(lambda y: beta(x0,y),y0,1)
    print(name,'grad',mp.nstr(bx,5),mp.nstr(by,5),' bxx',bxx,' byy',byy,' bxy',mp.nstr(bxy,5))
    for v,lab in ((bxx,'bxx'),(byy,'byy')):
        print('  ',lab,pslq([v,Om,pi*Om**3,pi**2*Om**5,pi**3*Om**7],maxcoeff=10**8,maxsteps=10**6))

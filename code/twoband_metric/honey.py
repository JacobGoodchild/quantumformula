# Integrated quantum metric of lower band of gapped graphene H = [[D, f],[f*, -D]],
# physical orbital positions: A at 0, B at (a1+a2)/3; Bravais constant a = 1.
from mpmath import mp, mpf, mpc, exp, cos, sin, pi, quad, sqrt
I=mpc(0,1)
mp.dps=25
def trg(k1,k2,D):
    ph=(k1+k2)/3
    f = exp(I*ph)*(1+exp(-I*k1)+exp(-I*k2))
    f1= I/3*f - I*exp(I*ph)*exp(-I*k1)        # d f/d k1
    f2= I/3*f - I*exp(I*ph)*exp(-I*k2)
    fx=f1+f2/2; fy=sqrt(3)/2*f2
    d=[f.real,f.imag,D]
    s=0
    for g in (fx,fy):
        dd=[g.real,g.imag,0]
        d2=sum(t*t for t in d); 
        s+=(sum(t*t for t in dd)*d2-sum(a*b for a,b in zip(d,dd))**2)/(4*d2*d2)
    return s
def M(D):
    D=mpf(D)
    return quad(lambda a: quad(lambda b: trg(a,b,D),[0,pi,2*pi]),[0,pi,2*pi])/(4*pi**2)
if __name__=='__main__':
    import sys
    print(M(sys.argv[1]))

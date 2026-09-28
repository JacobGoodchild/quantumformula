# Lieb lattice with four independent hoppings: x-bonds a1 (intra), a2 (inter); y-bonds b1, b2.
# Tr g = [Px b/2 + Py a/2 - a b/2]/S^2, a=|f|^2, b=|g|^2, S=a+b, Px=a1^2+a2^2, Py=b1^2+b2^2.
from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, pslq
def M4(a1,a2,b1,b2):
    Px=a1**2+a2**2; Py=b1**2+b2**2
    def fx(x):
        a=Px+2*a1*a2*cos(x); A=a+Py; B=2*b1*b2; D=A*A-B*B
        s1=1/sqrt(D); s2=A/D**mpf(1.5)
        # <b/S^2> = <(S-a)/S^2> = s1 - a s2 ; <a b/S^2> = a s1 - a^2 s2
        return Px/2*(s1-a*s2) + Py/2*a*s2 - (a*s1-a*a*s2)/2
    return quad(fx,[0,pi])/pi
if __name__=='__main__':
    mp.dps=50
    from fractions import Fraction as F
    for h in [(1,F(1,2),F(3,2),1),(2,1,1,F(1,3)),(F(3,2),F(1,2),F(2,3),F(1,3))]:
        a1,a2,b1,b2=[mpf(x.numerator)/x.denominator if isinstance(x,F) else mpf(x) for x in h]
        v=M4(a1,a2,b1,b2)
        Px=a1**2+a2**2; Py=b1**2+b2**2; kx=2*a1*a2/Px; ky=2*b1*b2/Py
        # Green-fn modulus for <1/(Px+Py+2a1a2 cx+2b1b2 cy)>
        A0=Px+Py; B1=2*a1*a2; B2=2*b1*b2; m=4*B1*B2/(A0**2-(B1-B2)**2); s=sqrt(A0**2-(B1-B2)**2)
        K=ellipk(m);E=ellipe(m)
        print([str(x) for x in h],'v=',mp.nstr(v,20),'m=',mp.nstr(m,15), pslq([v,K/pi,E/pi,K/pi*s,E/pi*s,K/pi/s,E/pi/s],maxcoeff=10**8,maxsteps=10**6))

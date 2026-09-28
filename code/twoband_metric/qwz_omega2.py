# <Omega^2> (BZ average of squared Berry curvature) of the QWZ lower band, exact inner integral
import sympy as sp
from mpmath import mp, mpf, quad, pi, cos
cx,sx,c,m,u,A,B,sy=sp.symbols('cx sx c m u A B sy')
d=sp.Matrix([sx,sy,m+cx+c]); dX=sp.Matrix([cx,0,-sx]); dY=sp.Matrix([0,c,-sy])
trip=d.dot(dX.cross(dY))
num=sp.expand(trip**2)
num=sp.expand(num.subs(sy**4,(1-c**2)**2).subs(sy**2,1-c**2))
num=sp.expand(num.subs(sx**4,(1-cx**2)**2).subs(sx**2,1-cx**2))
assert not num.has(sy) and not num.has(sx), num
numu=sp.expand(num.subs(c,(u-A)/B)); poly=sp.Poly(numu,u)
def avg_neg(n):   # <u^-n>, u = A + B cos t ; closed forms via derivatives of D^{-1/2}
    Ds=sp.Symbol('Ds')
    base=(A**2-B**2)**sp.Rational(-1,2)   # <1/u>
    # <u^-(n)> = (-1)^{n-1}/(n-1)! d^{n-1}/dA^{n-1} <1/u>
    return sp.simplify((-1)**(n-1)/sp.factorial(n-1)*sp.diff(base,A,n-1))
def avg_pos(n):
    return sum(sp.binomial(n,2*k)*A**(n-2*k)*B**(2*k)*sp.binomial(2*k,k)/4**k for k in range(n//2+1))
tot=0
for mon,co in zip(poly.monoms(),poly.coeffs()):
    j=mon[0]-3   # divide by u^3 (|d|^6 = u^3)
    tot+= co*(avg_pos(j) if j>=0 else avg_neg(-j))
inner=(tot/4).subs({A:2+m**2+2*m*cx,B:2*(m+cx)})
f=sp.lambdify((cx,m),inner,'mpmath')
def Om2(mm,dps=40):
    mp.dps=dps; mm=mpf(mm)
    return quad(lambda x: f(cos(x),mm),[0,pi])/pi
if __name__=='__main__':
    print(Om2(3,30),Om2('0.5',30))

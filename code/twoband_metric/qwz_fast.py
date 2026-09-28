import sympy as sp
from mpmath import mp, mpf, quad, pi, sqrt, cos, ellipk, ellipe, pslq
cx,sx,c,m,u,A,B=sp.symbols('cx sx c m u A B')
# ky-dependence: cy=c, sy^2=1-c^2. Build Tr g numerator/denominator symbolically.
sy=sp.Symbol('sy')
d=[sx,sy,m+cx+c]; dX=[cx,0,-sx]; dY=[0,c,-sy]
d2=sum(t*t for t in d)
num=0
for dd in (dX,dY):
    num+= (sum(t*t for t in dd)*d2 - sum(a*b for a,b in zip(d,dd))**2)
num=sp.expand(num)
# replace sy^2 -> 1-c^2 (only even powers appear), sx^2 -> 1-cx^2
num=sp.expand(num.subs(sy**4,(1-c**2)**2).subs(sy**2,1-c**2))
num=sp.expand(num.subs(sx**4,(1-cx**2)**2).subs(sx**2,1-cx**2))
assert not num.has(sy) and not num.has(sx)
Aexpr=2+m**2+2*m*cx; Bexpr=2*(m+cx)
# integrand = num/(4 u^2), u = A + B c
numu=sp.expand(num.subs(c,(u-A)/B))
poly=sp.Poly(numu,u)
coeffs={mon[0]:co for mon,co in zip(poly.monoms(),poly.coeffs())}
# <u^j> for j = -2,-1 ... we get terms u^{j-2}
def avg_upow(n):  # n integer
    if n==-1: return 1/sp.sqrt(A**2-B**2)
    if n==-2: return A/(A**2-B**2)**sp.Rational(3,2)
    # n>=0: <(A+Bc)^n> = sum_k binom(n,2k) A^{n-2k} B^{2k} binom(2k,k)/4^k
    return sum(sp.binomial(n,2*k)*A**(n-2*k)*B**(2*k)*sp.binomial(2*k,k)/4**k for k in range(n//2+1))
inner=sum(co*avg_upow(j-2) for j,co in coeffs.items())/4
inner=inner.subs({A:Aexpr,B:Bexpr})
f=sp.lambdify((cx,m),inner,'mpmath')
def M(mm,dps=40):
    mp.dps=dps
    mm=mpf(mm)
    return quad(lambda x: f(cos(x),mm),[0,pi])/pi
if __name__=='__main__':
    print(M(3,30), M(mpf(1)/2,30))

from aniso3 import parts
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp, itertools, math
mp.dps=60
def issq(x): return math.isqrt(x.numerator)**2==x.numerator and math.isqrt(x.denominator)**2==x.denominator
pts=[]
cands=[F(a,b) for b in range(2,10) for a in range(1,b)]
for k1,k2 in itertools.product(cands,cands):
    if k1==k2 or issq((1+k1)*(1+k2)): continue
    if len(pts)>=22: break
    if any(abs(k1-p[0])<F(1,50) and abs(k2-p[1])<F(1,50) for p in pts): continue
    K1=mpf(k1.numerator)/k1.denominator; K2=mpf(k2.numerator)/k2.denominator
    d1=sqrt((1-K1)/(1+K1)); d2=sqrt((1-K2)/(1+K2))
    gx,gy=parts(d1,d2); m=K1*K2; K=ellipk(m);E=ellipe(m); r=sqrt((1+K1)*(1+K2))
    rel=pslq([gx,K/pi,E/pi,K*r/pi,E*r/pi,1],maxcoeff=10**9,maxsteps=10**6)
    if rel and rel[0]:
        pts.append((k1,k2,[F(-rel[i],rel[0]) for i in range(1,6)]))
print(len(pts),'points')
x,y=sp.symbols('k1 k2')
mons=[x**i*y**j for i in range(5) for j in range(5) if i+j<=4]
for idx,name in enumerate(["K","E","K r","E r","const"]):
    cs=sp.symbols('c0:%d'%len(mons))
    eqs=[sum(c*mo.subs({x:sp.Rational(p[0].numerator,p[0].denominator),y:sp.Rational(p[1].numerator,p[1].denominator)}) for c,mo in zip(cs,mons))
         - sp.Rational(p[2][idx].numerator,p[2][idx].denominator)*(1-sp.Rational(p[0].numerator,p[0].denominator)*sp.Rational(p[1].numerator,p[1].denominator)) for p in pts]
    sol=sp.solve(eqs,cs,dict=True)
    if sol:
        poly=sum(sol[0].get(c,c)*mo for c,mo in zip(cs,mons))
        print(name, ':', sp.factor(poly), '/ (1 - k1 k2)')
    else: print(name,'no fit')
for p in pts[:10]:
    k1,k2,c=p
    print(k1,k2,'(1-m)=',1-k1*k2,' coeffs:',c, ' denoms:',[q.denominator for q in c])

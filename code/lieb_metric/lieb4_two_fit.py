from lieb4 import M4
from lieb4_family import basis
from mpmath import mp, mpf, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp, itertools
mp.dps=60
rows=[]
vals=[F(2),F(3),F(4),F(3,2),F(5),F(1,2),F(5,2),F(1,3),F(7,2),F(6)]
for a,b in itertools.product(vals,vals):
    if a==b or a==1 or b==1: continue
    A=mpf(a.numerator)/a.denominator; Bb=mpf(b.numerator)/b.denominator
    v=M4(A,mpf(1),Bb,mpf(1)); k2,pre=basis(A,1,Bb,1)
    rel=pslq([v,ellipk(k2)/pi/pre,ellipe(k2)/pi/pre],maxcoeff=10**12,maxsteps=10**6)
    if rel and rel[0]: rows.append((a,b,F(-rel[1],rel[0]),F(-rel[2],rel[0])))
print(len(rows),'points')
x,y=sp.symbols('a b')
okA=all(r[2]==(r[0]**2+r[1]**2+2)/(4*r[0]) for r in rows); print('alpha=(a^2+b^2+2)/(4a) on all points:',okA)
R=[(sp.Rational(r[0].numerator,r[0].denominator),sp.Rational(r[1].numerator,r[1].denominator),
    sp.Rational(r[3].numerator,r[3].denominator)*4*sp.Rational(r[0].numerator,r[0].denominator)/(sp.Rational(r[0].numerator,r[0].denominator)**2+sp.Rational(r[1].numerator,r[1].denominator)**2+2)) for r in rows]
for deg in (2,3,4):
    mons=[x**i*y**j for i in range(deg+1) for j in range(deg+1) if i+j<=deg]
    p=sp.symbols('p0:%d'%len(mons)); q=sp.symbols('q0:%d'%len(mons))
    eqs=[sum(pi_*m.subs({x:a,y:b}) for pi_,m in zip(p,mons)) - r*sum(qi*m.subs({x:a,y:b}) for qi,m in zip(q,mons)) for a,b,r in R]
    M=sp.Matrix([[e.coeff(v) for v in p+q] for e in eqs])
    ns=M.nullspace()
    print('deg',deg,'nullspace dim',len(ns))
    if len(ns)==1:
        v=ns[0]; P=sum(v[i]*mons[i] for i in range(len(mons))); Q=sum(v[len(mons)+i]*mons[i] for i in range(len(mons)))
        print('beta*4a/(a^2+b^2+2) =',sp.factor(P/Q)); break

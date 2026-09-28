from lieb4 import M4
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp, itertools, random
mp.dps=60
def basis(a1,a2,b1,b2):
    Px=a1**2+a2**2; Py=b1**2+b2**2; B1=2*a1*a2; B2=2*b1*b2
    c1=(B2-Px-Py)/B1; c2=(-B2-Px-Py)/B1
    return 2*(c1-c2)/((1-c1)*(-1-c2)), sqrt((1-c1)*(-1-c2))
random.seed(1)
vals=[F(2),F(3),F(1,2),F(3,2),F(5,2),F(1,3),F(4),F(2,3)]
rows=[]
for a,c,b in itertools.product(vals,vals,vals):
    if random.random()>0.5: continue
    A=mpf(a.numerator)/a.denominator; C=mpf(c.numerator)/c.denominator; Bb=mpf(b.numerator)/b.denominator
    v=M4(A,C,Bb,mpf(1)); k2,pre=basis(A,C,Bb,mpf(1))
    rel=pslq([v,ellipk(k2)/pi/pre,ellipe(k2)/pi/pre],maxcoeff=10**14,maxsteps=10**6)
    if rel and rel[0]: rows.append((a,c,b,F(-rel[1],rel[0]),F(-rel[2],rel[0])))
print(len(rows),'points')
X=sp.symbols('a c b')
def fit(idx,deg):
    mons=[X[0]**i*X[1]**j*X[2]**k for i in range(deg+1) for j in range(deg+1) for k in range(deg+1) if i+j+k<=deg]
    Mrows=[]
    for r in rows:
        sub={X[0]:sp.Rational(r[0].numerator,r[0].denominator),X[1]:sp.Rational(r[1].numerator,r[1].denominator),X[2]:sp.Rational(r[2].numerator,r[2].denominator)}
        val=sp.Rational(r[idx].numerator,r[idx].denominator)
        mv=[m.subs(sub) for m in mons]
        Mrows.append(mv+[-val*x for x in mv])
    ns=sp.Matrix(Mrows).nullspace()
    if len(ns)==1:
        v=ns[0]; n=len(mons)
        return sp.factor(sum(v[i]*mons[i] for i in range(n))/sum(v[n+i]*mons[i] for i in range(n))), len(ns)
    return None,len(ns)
for idx,name in ():
    for deg in (2,3,4):
        r,nd=fit(idx,deg); print(name,'deg',deg,'nullity',nd)
        if r is not None: print('  ',name,'=',r); break
bad=0
for a,c,b,al,be in rows:
    N=(a*a-c*c)**2+(b*b-1)**2; S1=(a-c)**2+(b-1)**2; S2=(a+c)**2+(b+1)**2
    if be != -al*N/(S1*S2): bad+=1
print('guess fails on',bad,'of',len(rows))

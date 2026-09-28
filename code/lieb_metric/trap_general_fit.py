from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellippi, pslq
from fractions import Fraction as F
import sympy as sp, itertools, random
mp.dps=50
def wB4(a,c,b,d):   # <|g|^2/(|f|^2+|g|^2)>, f: a,c (x);  g: b,d (y)
    Px=a*a+c*c; Py=b*b+d*d; B1=2*a*c; B2=2*b*d
    def fx(x):
        A=Px+B1*cos(x); return A/sqrt((A+Py)**2-B2**2)
    return 1-quad(fx,[0,pi])/pi
def basis(a,c,b,d):
    Px=a*a+c*c; Py=b*b+d*d; B1=2*a*c; B2=2*b*d
    c1=(B2-Px-Py)/B1; c2=(-B2-Px-Py)/B1
    return 2*(c1-c2)/((1-c1)*(-1-c2)), 2/(1-c1), sqrt((1-c1)*(-1-c2))
random.seed(2)
vals=[F(2),F(3),F(1,2),F(3,2),F(5,2),F(1,3),F(4),F(2,3)]
rows=[]
for a,c,b in itertools.product(vals,vals,vals):
    if random.random()>0.4: continue
    A,C,Bb=[mpf(x.numerator)/x.denominator for x in (a,c,b)]
    w=wB4(A,C,Bb,mpf(1)); m,n,pre=basis(A,C,Bb,mpf(1))
    rel=pslq([w-1,ellipk(m)/pi/pre,ellippi(n,m)/pi/pre],maxcoeff=10**14,maxsteps=10**6)
    if rel and rel[0]: rows.append((a,c,b,F(-rel[1],rel[0]),F(-rel[2],rel[0])))
print(len(rows),'points')
X=sp.symbols('a c b')
def fit(idx,deg):
    mons=[X[0]**i*X[1]**j*X[2]**k for i in range(deg+1) for j in range(deg+1) for k in range(deg+1) if i+j+k<=deg]
    Mr=[]
    for r in rows:
        sub={X[0]:sp.Rational(r[0].numerator,r[0].denominator),X[1]:sp.Rational(r[1].numerator,r[1].denominator),X[2]:sp.Rational(r[2].numerator,r[2].denominator)}
        val=sp.Rational(r[idx].numerator,r[idx].denominator); mv=[mm.subs(sub) for mm in mons]
        Mr.append(mv+[-val*x for x in mv])
    ns=sp.Matrix(Mr).nullspace()
    if len(ns)==1:
        v=ns[0]; k=len(mons); return sp.factor(sum(v[i]*mons[i] for i in range(k))/sum(v[k+i]*mons[i] for i in range(k))),1
    return None,len(ns)
for idx,name in ((3,'A (K/pre)'),(4,'B (Pi/pre)')):
    for deg in (1,2,3):
        r,nd=fit(idx,deg); print(name,'deg',deg,'nullity',nd)
        if r is not None: print('   ',name,'=',r); break

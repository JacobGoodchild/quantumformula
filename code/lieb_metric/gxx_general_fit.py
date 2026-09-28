from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, pslq
from fractions import Fraction as F
import sympy as sp, itertools, random
mp.dps=50
def gxx4(a,c,b,d):
    # g_xx = |f'|^2 |g|^2 / S^2, |f'|^2 = Px/2 - A/4 (A=|f|^2); ky-average exact
    Px=a*a+c*c; Py=b*b+d*d; B1=2*a*c; B2=2*b*d
    def fx(x):
        A=Px+B1*cos(x); U=A+Py; D=U*U-B2*B2; s1=1/sqrt(D); s2=U/D**mpf(1.5)
        return (Px/2-A/4)*(s1-A*s2)      # <b/S^2> = <1/S> - A<1/S^2>
    return quad(fx,[0,pi])/pi
def pre_m(a,c,b,d):
    R1=(a+c)**2+(b-d)**2; R2=(a-c)**2+(b+d)**2; return 16*a*b*c*d/(R1*R2), sqrt(R1*R2)
random.seed(3)
vals=[F(2),F(3),F(1,2),F(3,2),F(5,2),F(1,3),F(4),F(2,3)]
rows=[]
for a,c,b in itertools.product(vals,vals,vals):
    if random.random()>0.4: continue
    A,C,Bb=[mpf(x.numerator)/x.denominator for x in (a,c,b)]
    v=gxx4(A,C,Bb,mpf(1)); m,sq=pre_m(A,C,Bb,mpf(1))
    rel=pslq([v,ellipk(m)/pi/sq,ellipe(m)/pi/sq],maxcoeff=10**14,maxsteps=10**6)
    if rel and rel[0]: rows.append((a,c,b,F(-rel[1],rel[0]),F(-rel[2],rel[0])))
print(len(rows),'points')
X=sp.symbols('a c b')
def fit(idx,deg,den=1):
    mons=[X[0]**i*X[1]**j*X[2]**k for i in range(deg+1) for j in range(deg+1) for k in range(deg+1) if i+j+k<=deg]
    Mr=[]
    for r in rows:
        sub={X[0]:sp.Rational(r[0].numerator,r[0].denominator),X[1]:sp.Rational(r[1].numerator,r[1].denominator),X[2]:sp.Rational(r[2].numerator,r[2].denominator)}
        val=sp.Rational(r[idx].numerator,r[idx].denominator)*den.subs(sub) if den!=1 else sp.Rational(r[idx].numerator,r[idx].denominator)
        mv=[mm.subs(sub) for mm in mons]; Mr.append(mv+[-val*x for x in mv])
    ns=sp.Matrix(Mr).nullspace()
    if len(ns)==1:
        v=ns[0]; k=len(mons); return sp.factor(sum(v[i]*mons[i] for i in range(k))/sum(v[k+i]*mons[i] for i in range(k))/ (den if den!=1 else 1)),1
    return None,len(ns)
a,c,b=X; S1=(a-c)**2+(b-1)**2; S2=(a+c)**2+(b+1)**2
for idx,name in ():
    done=False
    for den in (1,S1*S2):
        for deg in (1,2,3):
            r,nd=fit(idx,deg,den); 
            if r is not None: print(name,'=',r,'(deg',deg,')'); done=True; break
        if done: break
    if not done: print(name,'no fit')
# polynomial fit of  (E coefficient) * S1*S2  with total degree <= 6
deg=6
mons=[X[0]**i*X[1]**j*X[2]**k for i in range(deg+1) for j in range(deg+1) for k in range(deg+1) if i+j+k<=deg]
Mr=[];rhs=[]
for r in rows:
    sub={X[0]:sp.Rational(r[0].numerator,r[0].denominator),X[1]:sp.Rational(r[1].numerator,r[1].denominator),X[2]:sp.Rational(r[2].numerator,r[2].denominator)}
    Mr.append([mm.subs(sub) for mm in mons]); rhs.append(sp.Rational(r[4].numerator,r[4].denominator)*(S1*S2).subs(sub))
Mm=sp.Matrix(Mr); bb=sp.Matrix(rhs)
sol=(Mm.T*Mm).LUsolve(Mm.T*bb)
res=Mm*sol-bb
print('max residual', max(abs(x) for x in res))
P=sum(sol[i]*mons[i] for i in range(len(mons)))
print('E coef * S1*S2 =', sp.factor(P))

from trap import wB
from mpmath import mp, mpf, sqrt, pi, ellipk, ellippi, pslq
from fractions import Fraction as F
import sympy as sp, math
mp.dps=50
rows=[]
for k1 in [F(1,3),F(1,2),F(2,3),F(1,4),F(3,4),F(2,5)]:
    for k2 in [F(1,5),F(2,7),F(3,5),F(5,6)]:
        if k1==k2: continue
        K1=mpf(k1.numerator)/k1.denominator;K2=mpf(k2.numerator)/k2.denominator
        d1=sqrt((1-K1)/(1+K1)); d2=sqrt((1-K2)/(1+K2))
        w=wB(d1,d2); m=K1*K2; n=K1*(1+K2)/(1+K1); s=sqrt((1+K2)/(1+K1))
        rel=pslq([w-1,s*ellipk(m)/pi,s*ellippi(n,m)/pi],maxcoeff=10**8,maxsteps=10**6)
        A=F(-rel[1],rel[0]);B=F(-rel[2],rel[0]); rows.append((k1,k2,A,B)); print(k1,k2,A,B)
x,y=sp.symbols('k1 k2')
for idx,name in ((2,'A'),(3,'B')):
    pts=[(r[0],r[1],sp.Rational(r[idx].numerator,r[idx].denominator)) for r in rows]
    # try A = poly(k1,k2)/((1+k1)*(1-k1*k2)) etc with small total degree
    for den in [1,(1+x),(1-x*y),(1+x)*(1-x*y),(1+y),(1+x)*(1+y),x,x*(1+y)]:
        mons=[x**i*y**j for i in range(4) for j in range(4) if i+j<=3]
        cs=sp.symbols('c0:%d'%len(mons))
        eqs=[sum(c*mo for c,mo in zip(cs,mons)).subs({x:sp.Rational(a.numerator,a.denominator),y:sp.Rational(b.numerator,b.denominator)})-v*den.subs({x:sp.Rational(a.numerator,a.denominator),y:sp.Rational(b.numerator,b.denominator)}) for a,b,v in pts]
        sol=sp.solve(eqs,cs,dict=True)
        if sol and all(s_.is_number for s_ in sol[0].values()) and len(sol[0])==len(cs):
            print(name,'=',sp.factor(sum(sol[0][c]*mo for c,mo in zip(cs,mons))/den)); break
    else: print(name,'no fit')

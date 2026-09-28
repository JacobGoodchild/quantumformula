from lieb4 import M4
from lieb4_family import basis
from mpmath import mp, mpf, pi, ellipk, ellipe, pslq
from fractions import Fraction as F
mp.dps=60
rows=[]
for a in [F(2),F(3),F(4),F(3,2),F(5)]:
    for b in [F(2),F(3),F(1,2),F(5,2)]:
        if a==b: continue
        A=mpf(a.numerator)/a.denominator; Bb=mpf(b.numerator)/b.denominator
        v=M4(A,mpf(1),Bb,mpf(1)); k2,pre=basis(A,1,Bb,1)
        rel=pslq([v,ellipk(k2)/pi/pre,ellipe(k2)/pi/pre],maxcoeff=10**10,maxsteps=10**6)
        if not rel or not rel[0]: print(a,b,'no rel',rel); continue
        al,be=F(-rel[1],rel[0]),F(-rel[2],rel[0]); rows.append((a,b,al,be)); print(a,b,al,be)
import pickle; pickle.dump(rows,open('lieb4_two.pkl','wb'))

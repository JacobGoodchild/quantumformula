import mpmath as mp, sympy as sp
from diamond import g_joyce
from fractions import Fraction as F
vals={4:'0.00061796260033398792810807558533012',5:'0.0004452715342203737731045364091334363',
      6:'0.0003377741377871855214956671944739493',7:'0.0002658128271949650938395312935021365',
      8:'0.0002150535512246489429953842654524992',10:'0.0001495985186874003856974715934279817',
      12:'0.0001103378647147750004610598626430863'}
mp.mp.dps=32
G=lambda s: -g_joyce(-(4+s)/4)/4
rows={}
for s,v in vals.items():
    v=mp.mpf(v); S=mp.mpf(s)
    rel=mp.pslq([v,1,G(S),mp.diff(G,S),mp.diff(G,S,2)],maxcoeff=10**6,maxsteps=10**6,tol=mp.mpf(10)**-22)
    rows[s]=[F(-rel[i],rel[0]) for i in (1,2,3,4)]
    print(s, rows[s])
x=sp.Symbol('s')
for idx,name in [(1,'a (G)'),(2,'b (G1)'),(3,'c (G2)')]:
    pts=[(k,sp.Rational(v[idx].numerator,v[idx].denominator)) for k,v in rows.items()]
    for deg in range(1,4):
        try: r=sp.polys.polyfuncs.rational_interpolate(pts[:2*deg+2],deg,X=x)
        except Exception: continue
        if all(r.subs(x,k)==y for k,y in pts): print(name, sp.factor(r)); break
    else:
        print(name,'poly fit', sp.factor(sp.interpolate(pts,x)))

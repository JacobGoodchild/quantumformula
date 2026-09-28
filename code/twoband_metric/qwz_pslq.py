from qwz_fast import M
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe, ellippi, pslq, identify
mp.dps=60
def roots(m):
    r1=-(m*m-2*m+2)/(2*(m-1)); r2=-(m*m+2*m+2)/(2*(m+1)); return r1,r2
for mm in [mpf(3),mpf(5)/2,mpf(4)]:
    v=M(mm,60)
    r=sorted([mpf(1),mpf(-1)]+list(roots(mm)),reverse=True)
    e1,e2,e3,e4=r
    k2=(e1-e2)*(e3-e4)/((e1-e3)*(e2-e4)); pre=sqrt((e1-e3)*(e2-e4))
    K=ellipk(k2);E=ellipe(k2)
    print(mm, v, 'k2',k2, 'roots',r)
    print(pslq([v,1,K/pi,E/pi,pre*K/pi,pre*E/pi],maxcoeff=10**6,maxsteps=10**6)); print(pslq([v,1,K/pi,E/pi,pre*K/pi,pre*E/pi,ellippi(k2*0+(e1-e2)/(e1-e3),k2)/pi,pre*ellippi((e1-e2)/(e1-e3),k2)/pi],maxcoeff=10**6,maxsteps=10**6))

# Integrated quantum metric of the lower band of a two-band model H = d(k).sigma
from mpmath import mp, mpf, sin, cos, pi, quad, ellipk, ellipe, pslq, sqrt
mp.dps=30
def trg_qwz(kx,ky,m):
    d=[sin(kx),sin(ky),m+cos(kx)+cos(ky)]
    dx=[cos(kx),0,-sin(kx)]; dy=[0,cos(ky),-sin(ky)]
    d2=sum(t*t for t in d)
    s=0
    for dd in (dx,dy):
        s+=(sum(t*t for t in dd)*d2-sum(a*b for a,b in zip(d,dd))**2)/(4*d2*d2)
    return s
def M(m):
    return quad(lambda x: quad(lambda y: trg_qwz(x,y,m),[0,pi]),[0,pi])/pi**2
if __name__=='__main__':
    for m in [mpf(3),mpf(1)/2,mpf(3)/2]:
        v=M(m); print(m,v)

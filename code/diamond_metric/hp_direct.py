# 30-digit direct BZ average (periodic trapezoid) of Tr g for gapped diamond, analytic derivatives
import mpmath as mp, sys
mp.mp.dps=32
I=mp.mpc(0,1)
nn=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
A=mp.matrix([[0,.5,.5],[.5,0,.5],[.5,.5,0]]); B=2*mp.pi*(A**-1).T
def trg(k,D):
    ph=[ (k[0]*a+k[1]*b+k[2]*c)/4 for a,b,c in nn]
    e=[mp.exp(I*p) for p in ph]
    f=sum(e); grads=[sum(I*v[x]/4*ei for v,ei in zip(nn,e)) for x in range(3)]
    d2=abs(f)**2+D*D; s=0
    for g in grads:
        s+=(abs(g)**2*d2-(mp.re(mp.conj(f)*g))**2)/(4*d2*d2)
    return s
def M(D,N):
    D=mp.mpf(D); tot=0
    for i in range(N):
        for j in range(N):
            for l in range(N):
                q=[mp.mpf(i)/N,mp.mpf(j)/N,mp.mpf(l)/N]
                k=[sum(q[r]*B[r,x] for r in range(3)) for x in range(3)]
                tot+=trg(k,D)
    return tot/N**3
if __name__=='__main__':
    from diamond import M_formula
    D=sys.argv[1]
    for N in map(int,sys.argv[2:]):
        v=M(D,N); print(N, v, v-M_formula(D)[0], flush=True)

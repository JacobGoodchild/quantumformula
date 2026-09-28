# <|Omega|^2> = sum_{i<j} <Omega_ij^2> for gapped diamond, periodic trapezoid in reduced coords
import mpmath as mp, sys
I=mp.mpc(0,1)
nn=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def run(s,N,dps):
    mp.mp.dps=dps; s=mp.mpf(s); D=mp.sqrt(s)
    A=mp.matrix([[0,.5,.5],[.5,0,.5],[.5,.5,0]]); B=2*mp.pi*(A**-1).T
    tot=0
    for i in range(N):
        for j in range(N):
            for l in range(N):
                q=[mp.mpf(i)/N,mp.mpf(j)/N,mp.mpf(l)/N]
                k=[sum(q[r]*B[r,x] for r in range(3)) for x in range(3)]
                e=[mp.exp(I*(k[0]*a+k[1]*b+k[2]*c)/4) for a,b,c in nn]
                f=sum(e); g=[sum(I*v[x]/4*ei for v,ei in zip(nn,e)) for x in range(3)]
                d2=abs(f)**2+s
                o=0
                for (x,y) in ((0,1),(1,2),(0,2)):
                    J=mp.im(mp.conj(g[x])*g[y]); o+=J*J
                tot+= s*o/(4*d2**3)
    return tot/N**3
if __name__=='__main__':
    s=sys.argv[1]
    for N in map(int,sys.argv[2:]): print(N, run(s,N,34), flush=True)

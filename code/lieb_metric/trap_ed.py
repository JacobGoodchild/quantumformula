# Independent physical check: finite L x L staggered Lieb torus, exact diagonalisation;
# long-time average of the CTQW return probability to an edge site B = sum over distinct
# eigenvalues E of (sum_{v in E} |v_B|^2)^2  -> w_B^2 as L grows (plus O(1/L^2) dispersive part).
import numpy as np, sys
def lat(L,d1,d2):
    idx=lambda x,y,s: 3*((x%L)*L+(y%L))+s
    N=3*L*L; H=np.zeros((N,N))
    for x in range(L):
        for y in range(L):
            A=idx(x,y,0); Bx=idx(x,y,1); Cy=idx(x,y,2)
            H[A,Bx]=H[Bx,A]=1+d1; H[Bx,idx(x+1,y,0)]=H[idx(x+1,y,0),Bx]=1-d1
            H[A,Cy]=H[Cy,A]=1+d2; H[Cy,idx(x,y+1,0)]=H[idx(x,y+1,0),Cy]=1-d2
    return H
def pbar(L,d1,d2):
    e,v=np.linalg.eigh(lat(L,d1,d2)); b=1
    w=v[b,:]**2; tot=0; i=0
    while i<len(e):
        j=i
        while j+1<len(e) and abs(e[j+1]-e[i])<1e-9: j+=1
        tot+=w[i:j+1].sum()**2; i=j+1
    flat=w[abs(e)<1e-9].sum()
    return tot, flat
if __name__=='__main__':
    d1,d2=float(sys.argv[1]),float(sys.argv[2])
    for L in map(int,sys.argv[3:]):
        t,f=pbar(L,d1,d2); print(L,'Pbar',t,'flat weight on B',f,flush=True)

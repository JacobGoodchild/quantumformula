import numpy as np
# pi-flux square lattice, gauge: hopping along x = -1, along y = -(-1)^x
def energy(L):
    N=L*L; H=np.zeros((N,N))
    idx=lambda x,y:(x%L)*L+(y%L)
    for x in range(L):
        for y in range(L):
            i=idx(x,y)
            H[i,idx(x+1,y)]+=-1; H[idx(x+1,y),i]+=-1
            s=(-1)**x
            H[i,idx(x,y+1)]+=-s; H[idx(x,y+1),i]+=-s
    e=np.linalg.eigvalsh(H)
    return e[:N//2].sum()/N
for L in [8,16,24,32]:
    print(L,energy(L))

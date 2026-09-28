# Numerical check that physical orbital positions minimise the integrated metric of gapped graphene:
# keep bonds fixed, place the B orbital at dA + r (fictitious), recompute M.
import numpy as np
a1=np.array([1.0,0]); a2=np.array([0.5,np.sqrt(3)/2]); dA=(a1+a2)/3
B=2*np.pi*np.linalg.inv(np.array([a1,a2])).T
def M(D,r,N=90,h=1e-5):
    nn=[dA,dA-a1,dA-a2]
    def P(k):
        f=sum(np.exp(1j*k@d) for d in nn)*np.exp(1j*k@r)   # extra phase = orbital shift r
        e,v=np.linalg.eigh(np.array([[D,f],[np.conj(f),-D]])); w=v[:,:1]; return w@w.conj().T
    t=0
    for i in range(N):
        for j in range(N):
            k=((i+0.5)/N)*B[0]+((j+0.5)/N)*B[1]
            for dk in (np.array([h,0]),np.array([0,h])):
                dP=(P(k+dk)-P(k-dk))/(2*h); t+=0.5*np.trace(dP@dP).real
    return t/N**2
m0=M(1.0,np.zeros(2)); print('r=0', m0)
for ang in [0,1,2,3,4,5]:
    r=0.05*np.array([np.cos(ang),np.sin(ang)]); print('r at angle',ang, M(1.0,r)-m0)

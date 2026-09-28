# Independent check: real-space-built Bloch Hamiltonian of gapped graphene (Cartesian k),
# numerical lower-band projector, metric from finite differences, BZ mesh in reduced coords.
import numpy as np
a1=np.array([1.0,0]); a2=np.array([0.5,np.sqrt(3)/2]); dA=(a1+a2)/3
nn=[dA, dA-a1, dA-a2]           # vectors from A to its three B neighbours
B=2*np.pi*np.linalg.inv(np.array([a1,a2])).T   # rows = reciprocal vectors b1,b2
def P(k,D):
    f=sum(np.exp(1j*k@d) for d in nn)
    H=np.array([[D,f],[np.conj(f),-D]])
    e,v=np.linalg.eigh(H); w=v[:,:1]; return w@w.conj().T
def integ(D,N=150,h=1e-5):
    t=0
    for i in range(N):
        for j in range(N):
            k=((i+0.5)/N)*B[0]+((j+0.5)/N)*B[1]
            for dk in (np.array([h,0]),np.array([0,h])):
                dP=(P(k+dk,D)-P(k-dk,D))/(2*h); t+=0.5*np.trace(dP@dP).real
    return t/N**2
if __name__=='__main__':
    for D in [1.0,3.0,0.5]: print(D,integ(D))

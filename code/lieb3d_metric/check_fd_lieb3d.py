# Independent: 4x4 Bloch Hamiltonian of 3D Lieb (corner A + three edge sites), numerical
# 2-fold flat-band projector, Tr g = (1/2) sum_i Tr(dP_i dP_i), midpoint mesh.
import numpy as np, sys
def P(k,d):
    f=[(1+d)*np.exp(1j*ki/2)+(1-d)*np.exp(-1j*ki/2) for ki in k]
    H=np.zeros((4,4),complex)
    for i in range(3): H[0,i+1]=f[i]; H[i+1,0]=np.conj(f[i])
    e,v=np.linalg.eigh(H); idx=np.argsort(abs(e))[:2]; w=v[:,idx]; return w@w.conj().T
def integ(d,N,h=1e-5):
    g=(np.arange(N)+0.5)*2*np.pi/N; t=0
    for a in g:
        for b in g:
            for c in g:
                k=np.array([a,b,c])
                for i in range(3):
                    e=np.zeros(3); e[i]=h
                    dP=(P(k+e,d)-P(k-e,d))/(2*h); t+=0.5*np.trace(dP@dP).real
    return t/N**3
if __name__=='__main__':
    d=float(sys.argv[1])
    for N in map(int,sys.argv[2:]): print(d,N,integ(d,N),flush=True)

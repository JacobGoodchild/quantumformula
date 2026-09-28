# Independent check: lower-band projector from numerical diagonalisation of H = d.sigma,
# metric = (1/2)Tr(dP dP) by central differences, midpoint k-mesh.
import numpy as np
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
def P(kx,ky,m):
    H=np.sin(kx)*sx+np.sin(ky)*sy+(m+np.cos(kx)+np.cos(ky))*sz
    e,v=np.linalg.eigh(H); w=v[:,:1]; return w@w.conj().T
def integ(m,N=200,h=1e-5):
    ks=(np.arange(N)+0.5)*2*np.pi/N; t=0
    for a in ks:
        for b in ks:
            for dx,dy in [(h,0),(0,h)]:
                dP=(P(a+dx,b+dy,m)-P(a-dx,b-dy,m))/(2*h); t+=0.5*np.trace(dP@dP).real
    return t/N**2
if __name__=='__main__':
    for m in [1.0,0.5,3.0,1.7]: print(m,integ(m))

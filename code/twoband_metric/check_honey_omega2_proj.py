# Independent: Berry curvature from the projector, Omega = -i Tr(P [dxP, dyP]) (gauge invariant),
# physical-position Bloch Hamiltonian, finite differences, midpoint mesh.
import numpy as np
from check_fd_honey import P, B
def om2(D,N=150,h=1e-5):
    t=0; tot=0
    ex=np.array([h,0]); ey=np.array([0,h])
    for i in range(N):
        for j in range(N):
            k=((i+0.5)/N)*B[0]+((j+0.5)/N)*B[1]
            p=P(k,D); dx=(P(k+ex,D)-P(k-ex,D))/(2*h); dy=(P(k+ey,D)-P(k-ey,D))/(2*h)
            om=(-1j*np.trace(p@(dx@dy-dy@dx))).real
            tot+=om*om; t+=om
    return tot/N**2, t/N**2
if __name__=='__main__':
    for s in [1.0,3.0]: print(s, om2(np.sqrt(s)))

# Independent: Fukui-Hatsugai plaquette Berry phases from numerical eigenvectors; <Omega^2> ~ mean (F/dk^2)^2
import numpy as np
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
def vecs(N,m):
    k=(np.arange(N))*2*np.pi/N
    V=np.empty((N,N,2),complex)
    for i,a in enumerate(k):
        for j,b in enumerate(k):
            H=np.sin(a)*sx+np.sin(b)*sy+(m+np.cos(a)+np.cos(b))*sz
            V[i,j]=np.linalg.eigh(H)[1][:,0]
    return V
def om2(N,m):
    V=vecs(N,m); dk=2*np.pi/N
    U1=np.einsum('ijk,ijk->ij',V.conj(),np.roll(V,-1,0)); U2=np.einsum('ijk,ijk->ij',V.conj(),np.roll(V,-1,1))
    F=np.angle(U1*np.roll(U2,-1,0)/np.roll(U1,-1,1)/U2)
    return np.mean((F/dk**2)**2), F.sum()/(2*np.pi)
if __name__=='__main__':
    for m in [1.0,0.5,3.0]:
        r=[om2(N,m)[0] for N in (200,400)]
        print(m, r, 'Richardson', (4*r[1]-r[0])/3)

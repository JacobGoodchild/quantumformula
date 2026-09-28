import numpy as np
a1=np.array([1.0,0]); a2=np.array([0.5,np.sqrt(3)/2]); dA=(a1+a2)/3; nn=[dA,dA-a1,dA-a2]
B=2*np.pi*np.linalg.inv(np.array([a1,a2])).T
def om2(N,D):
    V=np.empty((N,N,2),complex)
    for i in range(N):
        for j in range(N):
            k=(i/N)*B[0]+(j/N)*B[1]; f=sum(np.exp(1j*k@d) for d in nn)
            V[i,j]=np.linalg.eigh(np.array([[D,f],[np.conj(f),-D]]))[1][:,0]
    U1=np.einsum('ijk,ijk->ij',V.conj(),np.roll(V,-1,0)); U2=np.einsum('ijk,ijk->ij',V.conj(),np.roll(V,-1,1))
    F=np.angle(U1*np.roll(U2,-1,0)/np.roll(U1,-1,1)/U2)
    Acell=abs(np.linalg.det(B))/N**2          # k-space plaquette area
    return np.mean((F/Acell)**2)
if __name__=='__main__':
    for s in [1.0,3.0]:
        r=[om2(N,np.sqrt(s)) for N in (200,400)]; print(s,r,(4*r[1]-r[0])/3)

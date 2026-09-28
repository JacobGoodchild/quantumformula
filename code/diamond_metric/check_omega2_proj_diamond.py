# Independent: Berry curvature components from the numerical 2x2 projector, Omega_ij = -i Tr(P[d_iP, d_jP]),
# physical-position Hamiltonian of gapped diamond, finite differences, FCC primitive mesh.
import numpy as np
A=np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]]); B=2*np.pi*np.linalg.inv(A).T
nn=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/4.0
def Pk(K,D):
    f=np.exp(1j*K@nn.T).sum(1); E=np.sqrt(abs(f)**2+D*D)
    n=np.stack([f.real/E, -f.imag/E, D/E*np.ones_like(E)],1)   # H=[[D,f],[f*,-D]] -> d=(Re f,-Im f,D)
    return n
def om2(D,N=40,h=1e-5):
    g=(np.arange(N)+0.5)/N; I,J,L=np.meshgrid(g,g,g,indexing='ij')
    K=np.stack([I.ravel(),J.ravel(),L.ravel()],1)@B
    n=Pk(K,D); dn=[]
    for a in range(3):
        e=np.zeros(3); e[a]=h; dn.append((Pk(K+e,D)-Pk(K-e,D))/(2*h))
    tot=0
    for (x,y) in ((0,1),(1,2),(0,2)):
        om=0.5*np.einsum('ij,ij->i',n,np.cross(dn[x],dn[y]))   # lower band: Omega = -(1/2) n.(dn x dn); squared anyway
        tot+=om**2
    return tot.mean()
if __name__=='__main__':
    for s in [4.0,2.0]: print(s, om2(np.sqrt(s)), om2(np.sqrt(s),56))

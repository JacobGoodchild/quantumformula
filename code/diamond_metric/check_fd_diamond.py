import numpy as np
A=np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])          # FCC primitive vectors (cubic a=1)
B=2*np.pi*np.linalg.inv(A).T
nn=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/4.0
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
def Pm(K,D):   # vectorised projectors for array of k (n,3)
    f=np.exp(1j*K@nn.T).sum(1)
    E=np.sqrt(abs(f)**2+D*D)
    # lower band projector = (1 - dhat.sigma)/2
    dx,dy,dz=f.real/E,f.imag/E,D/E   # H = Re f sx - Im f sy? use H=[[D,f],[f*,-D]] -> d=(Re f, -Im f, D)
    return dx,-dy,dz
def integ(D,N=40,h=1e-5):
    g=(np.arange(N)+0.5)/N
    I,J,L=np.meshgrid(g,g,g,indexing='ij')
    K=np.stack([I.ravel(),J.ravel(),L.ravel()],1)@B
    tot=0
    for a in range(3):
        e=np.zeros(3);e[a]=h
        p=np.array(Pm(K+e,D)); m=np.array(Pm(K-e,D))
        dn=(p-m)/(2*h)                 # derivative of unit vector n
        tot+=(dn**2).sum(0)/4          # Tr g = |dn|^2/4 for P=(1-n.sigma)/2
    return tot.mean()
if __name__=='__main__':
    for D in [1.0,2.0,0.5,np.sqrt(8)]:
        print(D,integ(D),integ(D,56))

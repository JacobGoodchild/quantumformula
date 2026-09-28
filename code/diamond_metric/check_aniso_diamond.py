import numpy as np
A=np.array([[0,.5,.5],[.5,0,.5],[.5,.5,0]]); B=2*np.pi*np.linalg.inv(A).T
nn=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/4.0
def n_of(K,t,D):
    f=(np.exp(1j*K@nn.T)*t).sum(1); E=np.sqrt(abs(f)**2+D*D)
    return np.stack([f.real/E,-f.imag/E,D/E*np.ones_like(E)])
def M(t,D,N=48,h=1e-5):
    g=(np.arange(N)+0.5)/N; I,J,L=np.meshgrid(g,g,g,indexing='ij')
    K=np.stack([I.ravel(),J.ravel(),L.ravel()],1)@B; tot=0
    for a in range(3):
        e=np.zeros(3); e[a]=h; dn=(n_of(K+e,t,D)-n_of(K-e,t,D))/(2*h); tot+=(dn**2).sum(0)/4
    return tot.mean()
t=np.array([1.2,0.7,1.0,0.9]); T=(t**2).sum()
print('magic D^2=2T:', M(t,np.sqrt(2*T)), ' expect 1/64 =',1/64)
print('non-magic D=1:', M(t,1.0))

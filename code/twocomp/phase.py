# Two-component LLL vortex lattices (Mueller-Ho): minimise f = beta(tau) + alpha*beta12(tau, r0),
# beta12 = sum_k exp(-|k|^2/2) cos(k.r0), unit cell area 2*pi (l=1).
import numpy as np
from scipy.optimize import minimize
def recip(tau):
    x,y=tau.real,tau.imag
    a=np.sqrt(2*np.pi/y)          # side so that area a^2 y = 2 pi
    a1=np.array([a,0.]); a2=a*np.array([x,y])
    A=np.array([a1,a2]); B=2*np.pi*np.linalg.inv(A).T
    return A,B
N=8
mm,nn=np.meshgrid(np.arange(-N,N+1),np.arange(-N,N+1)); mm=mm.ravel(); nn=nn.ravel()
def f(p,alpha):
    x,y,u,v=p
    if y<=0.2: return 1e9
    A,B=recip(complex(x,y))
    K=np.outer(mm,B[0])+np.outer(nn,B[1]); k2=(K**2).sum(1)
    r0=u*A[0]+v*A[1]
    w=np.exp(-k2/2)
    return w.sum()+alpha*(w*np.cos(K@r0)).sum()
def best(alpha):
    res=[]
    for x0 in (0.0,0.25,0.5):
        for y0 in (0.87,1.0,1.4,2.0):
            for u0,v0 in ((0.5,0.5),(0.5,0.0),(1/3,1/3),(0.3,0.6)):
                r=minimize(f,[x0,y0,u0,v0],args=(alpha,),method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-13,'maxiter':20000})
                res.append((r.fun,r.x))
    res.sort(key=lambda t:t[0]); return res[0]
if __name__=='__main__':
    for alpha in [0.05,0.15,0.17,0.18,0.2,0.3,0.37,0.38,0.5,0.9,0.92,0.93,0.95]:
        fv,x=best(alpha)
        t=complex(x[0],x[1]); 
        print(alpha, round(fv,10), 'tau=',np.round(t,5),'|tau|=',round(abs(t),5),' r0=',np.round(x[2:]%1,4))

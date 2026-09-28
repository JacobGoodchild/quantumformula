# Gapped honeycomb with arbitrary bond hoppings t1,t2,t3 (ideal geometry, physical positions).
import numpy as np
from mpmath import mp, mpf, quad, pi, cos, sqrt
a1=np.array([1.0,0]); a2=np.array([0.5,np.sqrt(3)/2]); dA=(a1+a2)/3
nn=[dA, dA-a1, dA-a2]; B=2*np.pi*np.linalg.inv(np.array([a1,a2])).T
def M_fd(t,D,N=150,h=1e-5):
    def P(k):
        f=sum(tj*np.exp(1j*k@d) for tj,d in zip(t,nn))
        e,v=np.linalg.eigh(np.array([[D,f],[np.conj(f),-D]])); w=v[:,:1]; return w@w.conj().T
    s=0
    for i in range(N):
        for j in range(N):
            k=((i+0.5)/N)*B[0]+((j+0.5)/N)*B[1]
            for dk in (np.array([h,0]),np.array([0,h])):
                dP=(P(k+dk)-P(k-dk))/(2*h); s+=0.5*np.trace(dP@dP).real
    return s/N**2
def Gmp(t,D):
    # <1/(|f|^2 + D^2)>, f = t1 + t2 e^{i k1} + t3 e^{i k2}; k2-average exact
    t1,t2,t3=[mpf(x) for x in t]; D=mpf(D)
    def fk(k1):
        A=t1**2+t2**2+t3**2+D**2+2*t1*t2*cos(k1); R2=4*t3**2*(t1**2+t2**2+2*t1*t2*cos(k1))
        return 1/sqrt(A*A-R2)
    return quad(fk,[0,pi])/pi
def M_formula(t,D):
    S=sum(mpf(x)**2 for x in t); return mpf(1)/48+(3*S-mpf(D)**2)/48*Gmp(t,D)
if __name__=='__main__':
    mp.dps=30
    for t,D in [((1.3,0.8,1.0),1.0),((0.5,1.0,1.7),0.7),((1.3,0.8,1.0),float(np.sqrt(3*(1.69+0.64+1))))]:
        print(t,D,'FD',M_fd(t,D),'formula',M_formula(t,D))

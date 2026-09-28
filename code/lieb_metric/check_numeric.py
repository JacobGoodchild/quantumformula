# Independent check: integrated quantum metric of Lieb flat band from numerical diagonalisation
import numpy as np
def H(kx,ky,d):
    f=(1+d)*np.exp(1j*kx/2)+(1-d)*np.exp(-1j*kx/2)
    g=(1+d)*np.exp(1j*ky/2)+(1-d)*np.exp(-1j*ky/2)
    return np.array([[0,f,g],[np.conj(f),0,0],[np.conj(g),0,0]])
def P(kx,ky,d):
    e,v=np.linalg.eigh(H(kx,ky,d)); i=np.argmin(abs(e)); w=v[:,i:i+1]; return w@w.conj().T
def trg(kx,ky,d,h=1e-5):
    s=0
    for dx,dy in [(h,0),(0,h)]:
        dP=(P(kx+dx,ky+dy,d)-P(kx-dx,ky-dy,d))/(2*h)
        s+=0.5*np.trace(dP@dP).real
    return s
def integ(d,N=200):
    ks=(np.arange(N)+0.5)*2*np.pi/N
    return np.mean([trg(a,b,d) for a in ks for b in ks])
import sys
for d in [0.5,0.3]:
    print(d,integ(d,int(sys.argv[1]) if len(sys.argv)>1 else 120))

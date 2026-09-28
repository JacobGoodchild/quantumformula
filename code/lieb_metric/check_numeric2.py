# Independent check (no analytic reduction): flat-band projector from numerical
# diagonalisation of the 3x3 Bloch Hamiltonian, metric from finite differences of P,
# BZ average by midpoint rule. Optional orbital position s (B at (s,0), C at (0,s)).
import numpy as np, sys
def P(kx,ky,d1,d2,s):
    f=(1+d1)*np.exp(1j*kx*s)+(1-d1)*np.exp(-1j*kx*(1-s))
    g=(1+d2)*np.exp(1j*ky*s)+(1-d2)*np.exp(-1j*ky*(1-s))
    H=np.array([[0,f,g],[np.conj(f),0,0],[np.conj(g),0,0]])
    e,v=np.linalg.eigh(H); i=np.argmin(abs(e)); w=v[:,i:i+1]; return w@w.conj().T
def integ(d1,d2,s=0.5,N=160,h=1e-5):
    ks=(np.arange(N)+0.5)*2*np.pi/N; tot=0
    for a in ks:
        for b in ks:
            for dx,dy in [(h,0),(0,h)]:
                dP=(P(a+dx,b+dy,d1,d2,s)-P(a-dx,b-dy,d1,d2,s))/(2*h)
                tot+=0.5*np.trace(dP@dP).real
    return tot/N**2
if __name__=='__main__':
    print('aniso 3/7,5/11 :',integ(3/7,5/11))
    print('aniso 0.37,0.81:',integ(0.37,0.81))
    for s in [0.3,0.4,0.5,0.6]:
        print('position s',s, integ(0.5,0.5,s,N=100))

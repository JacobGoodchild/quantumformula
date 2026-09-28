from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, ellippi, pslq
mp.dps=50
def Ig(d1,d2):
    p1=1+d1*d1;q1=1-d1*d1;p2=1+d2*d2;q2=1-d2*d2
    def fx(x):
        a=2*p1+2*q1*cos(x)
        # b = 2p2+2q2 cos y ; S=a+b = A + B cos y
        A=a+2*p2; B=2*q2; D=A*A-B*B
        s1=1/sqrt(D); s2=A/D**mpf(1.5)          # <1/S>, <1/S^2>
        # term x: (p1-a/4)*(S-a)/S^2 = (p1-a/4)*(1/S - a/S^2)
        tx=(p1-a/4)*(s1-a*s2)
        # term y: (p2-b/4)*a/S^2 with b=S-a: (p2+a/4)*a/S^2 - a/(4S)
        ty=(p2+a/4)*a*s2 - a/4*s1
        return tx+ty
    return quad(fx,[0,pi])/pi
if __name__=='__main__':
    d1=mpf(1)/2; d2=mpf(1)/3
    v=Ig(d1,d2); print(v, Ig(d1,d1))
    # square lattice Green fn modulus: <1/(A0 + B1 cx + B2 cy)>, A0=2p1+2p2, B1=2q1,B2=2q2
    A0=2*(1+d1**2)+2*(1+d2**2); B1=2*(1-d1**2); B2=2*(1-d2**2)
    m=4*B1*B2/(A0**2-(B1-B2)**2)
    K=ellipk(m);E=ellipe(m)
    pre=1/sqrt(A0**2-(B1-B2)**2)
    print(pslq([v,1,K/pi,E/pi],maxcoeff=10**7,maxsteps=10**6))
    n1=(B1-B2)**2/(A0**2)  # guesses
    for n in [ -B1*B2*0+ (2*B2/(A0+B1-B2))**2 ]:
      pass
    s=sqrt(A0**2-(B1-B2)**2)
    print('s^2',s**2)
    print(pslq([v,1,K/pi,E/pi,s*K/pi,s*E/pi,K/pi/s,E/pi/s],maxcoeff=10**7,maxsteps=10**6))
    for n in [B1**2/A0**2, (B1-B2)/(A0), 2*B1/(A0+B1-B2),2*B2/(A0+B2-B1), (A0-B1+B2)/(A0+B1-B2)]:
      P=ellippi(n,m)
      print(n, pslq([v,1,K/pi,E/pi,P/pi],maxcoeff=10**7,maxsteps=10**6))

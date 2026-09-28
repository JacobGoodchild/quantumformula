from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, pslq
mp.dps=50
def parts(d1,d2):
    p1=1+d1*d1;q1=1-d1*d1;p2=1+d2*d2;q2=1-d2*d2
    def fx(x,which):
        a=2*p1+2*q1*cos(x); A=a+2*p2; B=2*q2; D=A*A-B*B
        s1=1/sqrt(D); s2=A/D**mpf(1.5)
        if which==0: return (p1-a/4)*(s1-a*s2)
        return (p2+a/4)*a*s2 - a/4*s1
    return [quad(lambda x: fx(x,w),[0,pi])/pi for w in (0,1)]
def formula(d1,d2):
    k1=(1-d1**2)/(1+d1**2); k2=(1-d2**2)/(1+d2**2); m=k1*k2
    return (2+k1+k2)/(8*sqrt((1+k1)*(1+k2)))*(2*ellipk(m)-ellipe(m))/pi, m, k1, k2
for d1,d2 in [(mpf(3)/7,mpf(5)/11),(mpf(1)/10,mpf(9)/10),(mpf('0.37'),mpf('0.81'))]:
    gx,gy=parts(d1,d2); f,m,k1,k2=formula(d1,d2)
    print(d1,d2,'sum',gx+gy,'formula',f,'diff',gx+gy-f)
    K=ellipk(m);E=ellipe(m); r=sqrt((1+k1)*(1+k2))
    print('  gx pslq',pslq([gx,1,K/pi,E/pi,K/pi/r,E/pi/r,K*r/pi,E*r/pi],maxcoeff=10**6,maxsteps=10**6))

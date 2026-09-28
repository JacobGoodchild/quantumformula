from mpmath import mp, mpf, sqrt, quad, pi, cos, ellipk, ellipe, pslq
mp.dps=40
def Ig(d):
    d=mpf(d); p=1+d*d; q=1-d*d; B=2*q
    def fx(x):
        a=2*p+2*q*cos(x); A=a+2*p; D=A*A-B*B
        return (p-a/2)/sqrt(D) + a*a/2*A/D**mpf(1.5)
    return quad(fx,[0,pi])/pi
for d in ['0.5','0.3']:
    print(d, Ig(d))
d=mpf(1)/2; v=Ig(d); k=(1-d*d)/(1+d*d)
K=ellipk(k*k); E=ellipe(k*k)
print(pslq([v,1,K/pi,E/pi],maxcoeff=10**6,maxsteps=10**6))

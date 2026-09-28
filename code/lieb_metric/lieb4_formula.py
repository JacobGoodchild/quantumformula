from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
def M4formula(a1,a2,b1,b2):
    A0=a1**2+a2**2+b1**2+b2**2; B1=2*a1*a2; B2=2*b1*b2
    s=sqrt(A0**2-(B1-B2)**2); m=4*B1*B2/s**2
    return A0/(4*s)*(2*ellipk(m)-ellipe(m))/pi
if __name__=='__main__':
    from lieb4 import M4
    mp.dps=50
    for h in [('1','0.5','1.5','1'),('2','1','1','0.333'),('1.5','0.5','0.667','0.333'),('0.37','1.91','2.3','0.2'),('1','0.001','1','0.001')]:
        a=[mpf(x) for x in h]; print(h, mp.nstr(M4formula(*a),25), mp.nstr(M4formula(*a)-M4(*a),3))

# Fully general staggered Lieb lattice: x-bonds a (intra), c (inter); y-bonds b (intra), d (inter).
from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
def Mgen(a,c,b,d):
    a,c,b,d=[mpf(x) for x in (a,c,b,d)]
    S1=(a-c)**2+(b-d)**2; S2=(a+c)**2+(b+d)**2; R1=(a+c)**2+(b-d)**2; R2=(a-c)**2+(b+d)**2
    k2=16*a*b*c*d/(R1*R2)
    return (a*a+b*b+c*c+d*d)/(2*pi*sqrt(R1*R2))*(ellipk(k2)-((a*a-c*c)**2+(b*b-d*d)**2)/(S1*S2)*ellipe(k2))
if __name__=='__main__':
    from lieb4 import M4
    from aniso import Ig
    mp.dps=50
    for h in [('1.37','0.61','2.2','0.93'),('0.3','1.7','1.1','2.9'),('3.3','0.2','0.5','0.45')]:
        v=Mgen(*h); print(h, mp.nstr(v,25), mp.nstr(v-M4(*[mpf(x) for x in h]),3))
    dx,dy=mpf('0.37'),mpf('0.81')
    print('reduces to F1(b):', mp.nstr(Mgen(1+dx,1-dx,1+dy,1-dy)-Ig(dx,dy),3))

from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
from lieb4_family import basis
def Mfam(a):
    a=mpf(a); k2,pre=basis(a,1,1,1)
    al=(a*a+3)/(4*a); be=-(a+1)**2*(a*a+3)/(4*a*(a*a+2*a+5))
    return (al*ellipk(k2)+be*ellipe(k2))/(pi*pre), k2, pre
if __name__=='__main__':
    from lieb4 import M4
    mp.dps=50
    for a in ['1.37','2.718281828','0.6','0.25']:
        v,k2,pre=Mfam(a); print(a, mp.nstr(v,25), mp.nstr(v-M4(mpf(a),1,1,1),3), 'k2',mp.nstr(k2,12),'pre^2',mp.nstr(pre**2,15))

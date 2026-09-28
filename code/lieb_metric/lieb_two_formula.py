from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
def Mab(a,b):
    a=mpf(a); b=mpf(b)
    R1=(a+1)**2+(b-1)**2; R2=(a-1)**2+(b+1)**2; S1=(a-1)**2+(b-1)**2; S2=(a+1)**2+(b+1)**2
    k2=16*a*b/(R1*R2)
    return (a*a+b*b+2)/(2*pi*sqrt(R1*R2))*(ellipk(k2)-((a*a-1)**2+(b*b-1)**2)/(S1*S2)*ellipe(k2))
if __name__=='__main__':
    from lieb4 import M4
    mp.dps=50
    for a,b in [('1.7','0.45'),('3.14159','2.71828'),('0.3','0.8'),('5','0.2')]:
        v=Mab(a,b); print(a,b, mp.nstr(v,25), mp.nstr(v-M4(mpf(a),1,mpf(b),1),3), mp.nstr(v-Mab(b,a),3))

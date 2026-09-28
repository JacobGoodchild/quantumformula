from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
def gxx_gen(a,c,b,d):
    a,c,b,d=[mpf(x) for x in (a,c,b,d)]
    S1=(a-c)**2+(b-d)**2; S2=(a+c)**2+(b+d)**2; R1=(a+c)**2+(b-d)**2; R2=(a-c)**2+(b+d)**2
    m=16*a*b*c*d/(R1*R2); T=a*a+b*b+c*c+d*d
    return ((a*a+c*c)/2*ellipk(m)-(a*a-c*c)**2*T/(2*S1*S2)*ellipe(m))/(pi*sqrt(R1*R2))
if __name__=='__main__':
    from gxx_general_fit import gxx4
    from lieb_general_formula import Mgen
    mp.dps=40
    for h in [('1.37','0.61','2.2','0.93'),('0.3','1.7','1.1','2.9')]:
        x=gxx_gen(*h); y=gxx_gen(h[2],h[3],h[0],h[1])
        print(h, mp.nstr(x,22), mp.nstr(x-gxx4(*[mpf(t) for t in h]),3), 'x+y-total', mp.nstr(x+y-Mgen(*h),3))

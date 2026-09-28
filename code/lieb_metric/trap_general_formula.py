from mpmath import mp, mpf, sqrt, pi, ellipk, ellippi
def wBgen(a,c,b,d):
    a,c,b,d=[mpf(x) for x in (a,c,b,d)]
    S1=(a-c)**2+(b-d)**2; R1=(a+c)**2+(b-d)**2; R2=(a-c)**2+(b+d)**2
    m=16*a*b*c*d/(R1*R2); n=4*a*c/R1
    return 1+2*((b-d)**2*ellipk(m)-S1*ellippi(n,m))/(pi*sqrt(R1*R2))
if __name__=='__main__':
    from trap_general_fit import wB4
    from trap_formula import wB_formula
    mp.dps=40
    for h in [('1.37','0.61','2.2','0.93'),('0.3','1.7','1.1','2.9'),('3.3','0.2','0.5','0.45')]:
        v=wBgen(*h); print(h, mp.nstr(v,25), mp.nstr(v-wB4(*[mpf(x) for x in h]),3))
    dx,dy=mpf('0.37'),mpf('0.81'); print('reduces to F9:', mp.nstr(wBgen(1+dx,1-dx,1+dy,1-dy)-wB_formula(dx,dy),3))

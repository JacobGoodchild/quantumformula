from mpmath import mp, mpf, sqrt, pi, ellipk
def Mhoney(D):
    s=mpf(D)**2; P=sqrt((1+s)*(9+s))
    k2=mpf(1)/2+(3-6*s-s*s)/(2*(1+s)*P)
    return mpf(1)/48+(9-s)*ellipk(k2)/(24*pi*(1+s)**mpf(0.75)*(9+s)**mpf(0.25))
if __name__=='__main__':
    mp.dps=50
    from honey_trap import M
    for D,N in [('1',120),('0.5',200),('2',90),('3',90),('3.7',90),('0.3',300)]:
        a=Mhoney(D); b=M(D,N,50); mp.dps=50
        print(D, mp.nstr(a,40), 'diff', mp.nstr(a-b,3))

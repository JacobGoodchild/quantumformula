from mpmath import mp, mpf, sqrt, pi, ellipk, ellippi
def wB_formula(d1,d2):
    k1=(1-d1**2)/(1+d1**2); k2=(1-d2**2)/(1+d2**2); m=k1*k2; n=k1*(1+k2)/(1+k1)
    return 1+((1-k2)*sqrt((1+k1)/(1+k2))*ellipk(m)-2*(1-m)/sqrt((1+k1)*(1+k2))*ellippi(n,m))/pi
if __name__=='__main__':
    from trap import wB
    mp.dps=40
    for d1,d2 in [(mpf('0.37'),mpf('0.81')),(mpf('0.1'),mpf('0.9')),(mpf('0.6'),mpf('0.6')),(mpf('0.05'),mpf('0.3'))]:
        print(d1,d2, mp.nstr(wB_formula(d1,d2),25), mp.nstr(wB_formula(d1,d2)-wB(d1,d2),3))

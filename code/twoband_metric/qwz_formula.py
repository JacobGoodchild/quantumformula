from mpmath import mp, mpf, pi, ellipk, ellippi
def Mformula(m):
    m=mpf(m); k2=8*(m*m-2)/m**4; n=4*(m-1)/m**2
    return mpf(1)/8 - (m*m-4*m+2)/(4*m*m*(m-1))*ellipk(k2)/pi - (m-2)**2/(4*m*(m-1))*ellippi(n,k2)/pi
if __name__=='__main__':
    from qwz_fast import M
    mp.dps=50
    for m in ['0.5','0.333','0.9','1.1','1.7','2.5','7','0.05','-0.5','-3']:
        a=Mformula(m); b=M(m,50)
        print(m, mp.nstr(a,25), mp.nstr(a-b,3))

from mpmath import mp, mpf, pi, sqrt, ellipk, diff
def Gs(s):
    s=mpf(s); P=sqrt((1+s)*(9+s)); k2=mpf(1)/2+(3-6*s-s*s)/(2*(1+s)*P)
    return 2*ellipk(k2)/(pi*(1+s)**mpf(0.75)*(9+s)**mpf(0.25))
def Om2formula(s):
    s=mpf(s); return -(s+9)/384*(Gs(s)+(s+3)*diff(Gs,s))
if __name__=='__main__':
    from honey_omega2_pslq import Om2s
    for s,N in [('0.3',260),('2.7',130),('0.05',500)]:
        v=Om2s(s,N,50); mp.dps=50; print(s, mp.nstr(Om2formula(s),30), mp.nstr(Om2formula(s)-v,3))

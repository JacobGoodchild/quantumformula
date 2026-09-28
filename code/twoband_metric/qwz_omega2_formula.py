from mpmath import mp, mpf, pi, ellipk, ellipe
def Om2formula(m):
    m=mpf(m); k2=8*(m*m-2)/m**4
    a=(m**4-2*m**2-4)/(32*m**2*(m**2-2)); b=-(m**6-6*m**4-4*m**2+16)/(32*(m*m-4)**2*(m*m-2))
    return (a*ellipk(k2)+b*ellipe(k2))/pi
if __name__=='__main__':
    from qwz_omega2 import Om2
    for mm in ['2.5','1.5','0.5','0.9','1.3','0.2','-0.7','-3.3','1.41']:
        v=Om2(mm,50); mp.dps=50; print(mm, mp.nstr(Om2formula(mm),25), mp.nstr(Om2formula(mm)-v,3))

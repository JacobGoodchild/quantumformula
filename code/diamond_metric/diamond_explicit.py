import mpmath as mp
def M_explicit(D):
    D=mp.mpf(D); s=D*D
    k2=mp.mpf(1)/2-((8+s)*mp.sqrt(16+s)-16*mp.sqrt(4+s))/(2*D**3)
    return mp.mpf(1)/64+(8-s)*(2*mp.sqrt(4+s)-mp.sqrt(16+s))*mp.ellipk(k2)**2/(16*mp.pi**2*D**3), k2
if __name__=='__main__':
    from diamond import M_formula
    mp.mp.dps=40
    for D in ['1','2','0.5','0.1','5']:
        a,k2=M_explicit(D); print(D, a, a-M_formula(D)[0], 'k2=',k2)
    for D in ['1e-2','1e-3']:
        print(D, M_explicit(D)[0]*mp.mpf(D))
    print('large', M_explicit(100)[0]*10**4)

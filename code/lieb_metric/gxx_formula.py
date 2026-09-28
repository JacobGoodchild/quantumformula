from mpmath import mp, mpf, sqrt, pi, ellipk, ellipe
def gxx(d1,d2):
    k1=(1-d1**2)/(1+d1**2); k2=(1-d2**2)/(1+d2**2); m=k1*k2
    return sqrt((1+k2)/(1+k1))/(16*pi)*(4*ellipk(m)-(1-k1)*(2+k1+k2)/(1-m)*ellipe(m))
if __name__=='__main__':
    from aniso3 import parts, formula
    mp.dps=50
    for d1,d2 in [(mpf('0.37'),mpf('0.81')),(mpf('0.1'),mpf('0.9')),(mpf(3)/7,mpf(5)/11),(mpf('0.05'),mpf('0.6'))]:
        gx,gy=parts(d1,d2)
        print(mp.nstr(gxx(d1,d2),25), mp.nstr(gxx(d1,d2)-gx,3), mp.nstr(gxx(d2,d1)-gy,3), mp.nstr(gxx(d1,d2)+gxx(d2,d1)-formula(d1,d2)[0],3))

# Gapped diamond lattice (staggered mass Delta), physical orbital positions, cubic a = 1.
import mpmath as mp
def g_joyce(t):   # Joyce: <1/(t - (c1c2+c2c3+c3c1))>, analytic continuation to t < -1 tested below
    x = 3 / t
    k2 = mp.mpf(1)/2 - (mp.sqrt(3)/2) * (3 + x)**(-1.5) * (4*x + (3 - x)*mp.sqrt(1 - x))
    return 3**mp.mpf(1.5) * (3 + x)**(-1.5) * (2 - mp.sqrt(1 - x)) * (2*mp.ellipk(k2)/mp.pi)**2 / t
def g_num(t):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2); A = t - c1*c2; B = c1 + c2
        return mp.sign(A)/mp.sqrt(A*A - B*B)
    return mp.quad(f, [0, mp.pi/2, mp.pi], [0, mp.pi/2, mp.pi]) / mp.pi**2
def M_formula(D):
    D=mp.mpf(D); w=4+D*D
    G=-g_joyce(-w/4)/4          # <1/(4u + w)>
    return mp.mpf(1)/64 + (8-D*D)*G/64, G
if __name__=='__main__':
    mp.mp.dps=25
    for t in ['-1.5','-3']:
        t=mp.mpf(t); print('joyce vs num at t=',t, g_joyce(t)-g_num(t), g_joyce(t))
    for D in ['1','2','0.5']:
        print(D, M_formula(D))

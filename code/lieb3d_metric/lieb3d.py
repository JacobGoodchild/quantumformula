import mpmath as mp
def gsq(z):  # <1/(z - c1 - c2)>, z > 2
    return 2/(mp.pi*z)*mp.ellipk(4/z**2)
def G(w):    # simple cubic <1/(w - c1 - c2 - c3)>, w >= 3
    return mp.quad(lambda t: gsq(w-mp.cos(t)), [0, mp.pi])/mp.pi
def M3(d):
    d=mp.mpf(d); p=1+d*d; q=1-d*d; w=3*p/q
    g0=G(w); dg=mp.diff(G,w)
    return w/4*g0 + (w*w-9)/12*dg
def watson():
    g=mp.gamma
    return mp.sqrt(6)/(32*mp.pi**3)*g(mp.mpf(1)/24)*g(mp.mpf(5)/24)*g(mp.mpf(7)/24)*g(mp.mpf(11)/24)
if __name__=='__main__':
    mp.mp.dps=25
    for d in ['0.5','0.3','0.1','0.01']:
        print(d, M3(d))
    print('3/4 Watson =', 3*watson()/4, ' G(3)=',G(3),watson())

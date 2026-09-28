import mpmath as mp
from diamond import g_joyce
def G(s): return -g_joyce(-(4+s)/4)/4
def Om2(s):
    s=mp.mpf(s)
    return ((s+8)*G(s)+(5*s*s+76*s+128)*mp.diff(G,s)+2*s*(s+4)*(s+16)*mp.diff(G,s,2))/1024
if __name__=='__main__':
    mp.mp.dps=32
    for s,v in [(2,'0.001594810651686037898371758075373164'),(3,'0.0009273628522294665432807546306368192')]:
        print(s, Om2(s), Om2(s)-mp.mpf(v))

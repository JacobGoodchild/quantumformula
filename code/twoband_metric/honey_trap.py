from mpmath import mp, mpf, pi
from honey import trg
import sys
def M(D,N,dps):
    mp.dps=dps; D=mpf(D)
    h=2*pi/N; s=0
    for i in range(N):
        a=i*h
        for j in range(N):
            s+=trg(a,j*h,D)
    return s/N**2
if __name__=='__main__':
    D=sys.argv[1]
    for N in [int(x) for x in sys.argv[2:]]:
        print(N, M(D,N,40))

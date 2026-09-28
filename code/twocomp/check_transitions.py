# Independent check: direct numerical minimisation (numpy, double precision) along the relevant
# one-parameter families on either side of the predicted transition points.
import numpy as np
from scipy.optimize import minimize_scalar
from phase import f
a2=0.37321550790897782596; a3=0.92564969740393552904
def rhombic_min(alpha):   # tau on unit circle, tau = e^{i theta}, r0 = cell centre
    g=lambda th: f([np.cos(th),np.sin(th),0.5,0.5],alpha)
    r=minimize_scalar(g,bounds=(np.pi/2-0.6,np.pi/2+0.6),method='bounded',options={'xatol':1e-12})
    return np.cos(r.x)
def rect_min(alpha):      # tau = i y, r0 = centre
    g=lambda y: f([0.0,y,0.5,0.5],alpha)
    r=minimize_scalar(g,bounds=(0.6,1.6),method='bounded',options={'xatol':1e-12})
    return r.x
for d in (-2e-3,-5e-4,5e-4,2e-3):
    print('alpha2%+g: optimal Re(tau) on |tau|=1 = %.6f' % (d, abs(rhombic_min(a2+d))))
for d in (-2e-3,-5e-4,5e-4,2e-3):
    y=rect_min(a3+d); print('alpha3%+g: optimal y = %.6f (1/y=%.6f)' % (d, y, 1/y))

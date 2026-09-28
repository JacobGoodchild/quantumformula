import sympy as sp
r=sp.symbols('r',positive=True); pi=sp.pi
j1=sp.sin(pi*r)/(pi*r)**2-sp.cos(pi*r)/(pi*r)
phi=-r*sp.cos(pi*r)/(2*pi)-j1/2
# check ODE and BC
L=sp.diff(phi,r,2)+2*sp.diff(phi,r)/r-2*phi/r**2+pi**2*phi-2*sp.sin(pi*r)
print('ODE residual:',sp.simplify(L),' phi(1)=',sp.simplify(phi.subs(r,1)),' phi(0+)=',sp.limit(phi,r,0))
alpha=sp.simplify(-sp.Rational(4,3)*sp.integrate(r**2*sp.sin(pi*r)*phi,(r,0,1)))
print('alpha (units m e^2 R^4/hbar^2) =',alpha, '=', sp.N(alpha,40))

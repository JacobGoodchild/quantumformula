import sympy as sp
r=sp.symbols('r',positive=True); pi=sp.pi
cs=sp.symbols('c0:12')
basis=[r**k*sp.sin(pi*r) for k in range(-2,4)]+[r**k*sp.cos(pi*r) for k in range(-2,4)]
phi=sum(c*b for c,b in zip(cs,basis))
L=sp.diff(phi,r,2)+2*sp.diff(phi,r)/r-2*phi/r**2+pi**2*phi - 2*sp.sin(pi*r)
L=sp.expand(sp.simplify(L*r**4))
eqs=[]
for trig in (sp.sin(pi*r),sp.cos(pi*r)):
    part=sp.expand(L).coeff(trig)
    eqs+=sp.Poly(part,r).coeffs()
sol=sp.solve(eqs,cs,dict=True)
print(sol)
if sol:
    print(sp.simplify(phi.subs(sol[0])))

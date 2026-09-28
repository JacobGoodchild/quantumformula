# NOTES — working log for quantumformula

## Status summary
Phase 1 (2026-09-28, extended to ~1h15 at user's request): 7 verified results, all in QUANTUM
GEOMETRY of lattice bands (quantum metric = minimal Wannier spread / "quantum weight"; Berry-curvature
fluctuations). Formulas 1-7 in formulas.txt. Key technique: (i) inner BZ integral exact when the
denominator is linear in one cosine; (ii) DIVERGENCE-THEOREM trick: if |grad f|^2 is affine in
|f|^2 and Laplacian(eps) = -L^2 eps, averages of |grad eps|^2 h(eps) reduce to lattice Green
functions -> closed forms (honeycomb, diamond, Lieb in any d).
Plan for Phase 2: (1) push the quantum-geometry line to other flat-band lattices (dimerised
kagome / dice / checkerboard / 3D perovskite-Lieb, Lieb with next-neighbour terms), and to
related band-geometric quantities with rational-in-Bloch-function integrands (Wannier spreads
of filled bands, flat-band Berry-curvature-free metrics, flat-band exciton sizes);
(2) CTQW long-time averages on flat-band and decorated lattices; (3) impurity/vacancy
quantities from lattice Green functions; (4) quick scans of spin-chain and quantum-info items.

## Log

### 2026-09-28 Phase 1
- Scout: pi-flux square lattice half-filling energy = 3F2(-1/4,1/4,1/2;1,1;1) = 0.958091398682850128...
  (series from binomial expansion; no Gamma(1/4)/K/E closed form found by PSLQ). Probably known
  numerically in the flux-phase literature; parked.
- Scout: local static susceptibility chi_00 of half-filled square lattice = (2/pi^3) Int_0^1 K(k)^2/k' dk
  = 0.559318193582093534...; Int K^2/k' = 2 Int K K'. 1D chain gives exactly 1/2. No closed form
  found (these are L-value type constants; cf. Rogers-Wan-Zucker moments of elliptic integrals).
  Time-integrated CTQW return probability on Z^2 = pi Int rho^2 = (4/pi^3) Int_0^1 K^2 dk = same class.
- HIT: Lieb flat-band quantum metric -> FORMULA 1 (isotropic, anisotropic, minimal over positions).
  Leads: individual components g_xx for anisotropic case are K,E with algebraic coefficients
  (PSLQ found one); 3D Lieb (perovskite) flat bands -> SC lattice Green function; other gapped
  flat-band lattices (decorated/dimerised kagome, dice with staggering, checkerboard variants).

- HITS (all verified, see formulas.txt):
  F2 QWZ integrated metric: 1/8 + K, Pi (third kind) terms, k^2 = 8(m^2-2)/m^4, n = 4(m-1)/m^2, all m.
  F3 gapped graphene metric = 1/48 + (9-D^2) G_tri(3+D^2)/48, exactly 1/48 at D = 3.
  F4 general simplex bipartite lattices; magic mass D*^2 = z^2/(z-2); diamond via Joyce FCC (K^2).
  F5 3D Lieb: M = (w/4)G + ((w^2-9)/12)G'; = Watson/4 at delta = 0 (finite in 3D). General d formula.
  F6 QWZ <Omega^2> in K, E only (rational coefficient functions by exact interpolation).
  F7 gapped graphene <Omega^2> = -(s+9)/384 [G + (s+3)G'].
- Pitfalls met: python literal 1j/3 silently in double precision inside mpmath code (capped
  agreement at 1e-34) -> use mpc(0,1). Plaquette Berry-phase method fails with physical-position
  (non-periodic) gauges -> use projector formula. Joyce FCC formula continues correctly below band.
- LEADS / TODO:
  * Derive F7 by hand (Jacobian of k -> f); F6 likewise.
  * <Omega^2>, metric for gapped diamond (3D curvature vector) and hyperdiamond.
  * Anisotropic components g_xx separately for Lieb (superfluid-weight tensor) - PSLQ found K,E with
    algebraic coefficients at one point.
  * Haldane / Kane-Mele metric (d_z k-dependent; harder).
  * Breathing kagome / other flat bands if gapped; checkerboard; dice with staggering.
  * Minimal (over orbital positions) metric for gapped graphene and QWZ (embedding dependence).
  * Other directions not yet touched: CTQW long-time averages, impurity phase shifts, spin chains,
    quantum-info items (see Directions list in the mission).

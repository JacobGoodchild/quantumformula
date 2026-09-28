# NOTES — working log for quantumformula

## Status summary
END OF PHASE 1 (2026-09-28, ~17:00-18:00 UTC; the user extended it to about 1h15).
9 verified results, all in formulas.txt, none found in the literature searched:
  F1 Lieb flat band quantum metric: (2K-E)/(4pi); anisotropic total and x/y components; minimal metric.
  F2 Qi-Wu-Zhang integrated quantum metric, one formula (K, Pi) valid for all m.
  F3 gapped graphene metric = 1/48 + (9-D^2) G_tri/48, exactly 1/48 at D = 3 (also the minimal metric).
  F4 simplex bipartite lattices: general reduction, magic mass D*^2 = z^2/(z-2), gapped diamond (K^2).
  F5 3D Lieb flat bands via the simple-cubic Green function; Watson/4 at zero staggering; general d.
  F6 QWZ Berry-curvature fluctuations <Omega^2> (K, E only).
  F7 gapped graphene <Omega^2> = -(s+9)/384 [G + (s+3)G'].
  F8 gapped diamond <|Omega|^2> = [(s+8)G + (5s^2+76s+128)G' + 2s(s+4)(s+16)G'']/1024.
  F9 quantum-walk trapping fraction on the anisotropic staggered Lieb lattice (K, Pi).
Main technique: the inner BZ integral is exact when the denominator is linear in one cosine;
divergence-theorem reduction to lattice Green functions; PSLQ over {1,K,E,Pi} or {1,G,G',G''};
exact rational interpolation of coefficient functions; checks at held-out parameters and with
independent finite-difference / plaquette / exact-diagonalisation codes.
Negative results: pi-flux energy (a 3F2, no nicer closed form); square-lattice local susceptibility and
time-integrated CTQW return probability (L-value-type constants, no closed form);
breathing kagome flat band is NOT isolated (touches the next band), so skipped.
Plan for Phase 2 (overnight): hand derivations of F6-F8; Haldane/Kane-Mele metric; anisotropic
honeycomb (strained graphene) metric; minimal metrics with embedding freedom; then branch to
CTQW long-time averages on other flat-band lattices, impurity phase shifts from Green functions,
XX/XY chain items and quantum-information items.

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
- Breathing kagome (t_up=1, t_down=0.5 or 0.8): flat band at -(t_up+t_down) still touches the middle
  band -> integrated metric diverges; no gapped formula. (negative)
- Late checks: diamond Delta=1 direct N=50 diff 2.4e-22; 3D Lieb delta=0 N=40 extrapolates to 0.37912.

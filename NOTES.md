# NOTES — working log for quantumformula

## Status summary
(Phase 1 started 2026-09-28.)

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

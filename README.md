# quantumformula

Computational search for new, exactly verified formulas in quantum mechanics and
quantum physics (quantum walks, free fermions on lattices, lattice Green functions,
spin chains, quantum information, one-body spectral problems, solvable disorder).

**This work is AI-assisted:** the computations, derivations and write-ups were carried
out by Claude Code (Anthropic's AI coding agent), directed by Jacob Goodchild.

- **[formulas.txt](formulas.txt)** — every result, with its index, set-up, exact closed
  form, derivation sketch, independent verification numbers and an honest novelty status.
- **[NOTES.md](NOTES.md)** — working log: what was tried, negative results, open leads.
- `code/verify_all.py` — re-checks every formula against an independent computation (~10 s).
- `code/` — the Python scripts (numpy / scipy / mpmath / sympy) that produced and checked
  each result.

Method: compute a quantity to 30–80 digits, identify a closed form with integer-relation
detection (PSLQ) over a sensible basis of constants, then confirm it by at least one
independent computation (exact diagonalisation on growing systems, simulation, or a
separately written numerical integral). Novelty labels come from web searches only and
should be confirmed by experts.

## Results so far (details and honest novelty labels in formulas.txt)
1. Lieb-lattice flat band: exact quantum metric for arbitrary bond strengths; minimal version.
2. Qi-Wu-Zhang Chern insulator: exact integrated quantum metric for every mass m.
3. Gapped graphene: exact quantum metric; exactly 1/48 at gap parameter Delta = 3t.
4. General "simplex" lattices and the gapped diamond lattice; a universal "magic mass".
5. 3D Lieb (perovskite) flat bands: metric = Watson's integral / 4 at zero staggering.
6-8. Exact Berry-curvature fluctuations for Qi-Wu-Zhang, gapped graphene and gapped diamond.
9. Exact fraction of a quantum walker trapped forever on the anisotropic Lieb lattice.
10. Strained graphene / distorted diamond: the magic-mass value is universal (any bond strengths).
11. Abrikosov vortex lattices: exact beta in all Landau levels, higher moments, shear stiffness.
12. Rotating two-component condensates: exact vortex-lattice phase boundaries.

Predecessor project: https://github.com/JacobGoodchild/findformula

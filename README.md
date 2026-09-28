# quantumformula

Computational search for new, exactly verified formulas in quantum mechanics and
quantum physics (quantum walks, free fermions on lattices, lattice Green functions,
spin chains, quantum information, one-body spectral problems, solvable disorder).

**This work is AI-assisted:** the computations, derivations and write-ups were carried
out by Claude Code (Anthropic's AI coding agent), directed by Jacob Goodchild.

- **[formulas.txt](formulas.txt)** — every result, with its index, set-up, exact closed
  form, derivation sketch, independent verification numbers and an honest novelty status.
- **[NOTES.md](NOTES.md)** — working log: what was tried, negative results, open leads.
- `code/` — the Python scripts (numpy / scipy / mpmath / sympy) that produced and checked
  each result.

Method: compute a quantity to 30–80 digits, identify a closed form with integer-relation
detection (PSLQ) over a sensible basis of constants, then confirm it by at least one
independent computation (exact diagonalisation on growing systems, simulation, or a
separately written numerical integral). Novelty labels come from web searches only and
should be confirmed by experts.

Predecessor project: https://github.com/JacobGoodchild/findformula

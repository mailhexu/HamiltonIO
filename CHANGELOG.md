# Changelog

## Unreleased

### Wannier90 Wigner-Seitz weights (ws-weights epic, stories 026-030)

- **Fixed**: `WannierHam.gen_ham` now divides by the Wannier-Seitz degeneracy
  `ndegen(R)` (previously the division was commented out — a latent bug).
  `_hr.dat` and `_tb.dat` both store raw `ham_r`; the reader must divide, per
  the Wannier90 convention (`hamiltonian.F90:478`, `plot.F90:346`). **This
  changes every Wannier90-derived band structure and TB2J exchange value; the
  new values are correct.**
- **Added**: per-orbital-pair Wigner-Seitz correction from `_wsvec.dat`
  (Wannier90 `use_ws_distance`). `read_from_wannier_dir` auto-detects
  `{prefix}_wsvec.dat` and applies scheme 2 when the header is `.true.` — no
  user flag or API change. Matters when Wannier centres are off high-symmetry
  positions.
- **Added**: `parse_wsvec(filename)` and `validate_ws_weights(hr_path, mp_grid,
  wsvec_path=None)` in `HamiltonIO.wannier`.
- **Internal**: `WannierHam` gains `use_ws`/`ws_shifts` attributes (set by
  auto-detection, not a public knob). `shift_position` propagates ws state by
  re-keying to shifted R-vectors.

### Bug fixes

- `gen_ham`: the commented-out division used the integer loop index, but
  `R_degens` is dict-keyed by R-tuple. Now uses `R_degens.get(R, 1)` (correct
  key, and no KeyError on partial dicts / after `shift_position`).
- `wannier/__init__.py`: fixed pre-existing F401 on re-exports via redundant
  aliases.

### Known limitation

`gen_ham` `convention=1` has a separate pre-existing shape bug
(`np.dot(k, R + rjminusri)` raises for nbasis > 0), unrelated to ws weights.
`convention=2` (the default, used by TB2J) is correct.

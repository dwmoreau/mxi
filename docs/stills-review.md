# Stills indexing: where it departs from dials

`mxi_stills` reimplements, in new C++, the algorithms `dials.stills_process` uses to index a
still. Each entry here is a place where it does something different from dials, or from the
paper the algorithm comes from. An entry is written when the change is made, and gives its
reason: a measured cost, with the run that measured it, or a stated weakness of the model.

| # | date | what differs | from (dials code or paper) | reason | measured effect |
|---|---|---|---|---|---|
| 1 | 2026-10-11 | Refinement stops, crystal unchanged, when the equations number fewer than the free parameters (3N < P) | dials' LevMar `dof = 3N - P` is a `size_t` (scitbx `normal_equations.h`), so it wraps and dials iterates on an underdetermined system | The wrap is an unsigned-arithmetic accident, and the result is a model fitted to fewer equations than parameters | 1 LevMar run on each of lysozyme and luciferase in the R6 replay (spec: 16 and 37 on the full sets). Each is the first refinement of an outlier pass with N = 1 or 2 spots, whose Sauter-Poon step then fails for N < 10, so no frame's outcome changes |
| 2 | 2026-10-11 | No detector parameterisation is built, so the panel frames are not recomposed on every refinement | dials builds a hierarchical detector parameterisation on every refiner build, whose `compose()` rewrites the copied panels at about 1e-12 mm and is written back to the experiment (spec 9.2) | The detector is fixed on every configuration (STATUS 3, decided 2026-10-10); the parameterisation exists only to be discarded | None resolvable here: on 5,186 captured runs the starting objective and RMSDs agree with dials to 4.4e-12 relative |

## What follows dials exactly

**Spot centroids to reciprocal space** (`src/stills/geometry.cc`). This is dials'
`centroid_px_to_mm_panel` and `map_centroids_to_reciprocal_space` for a still, with
dxtbx's `parallax_correction_inv` and `get_lab_coord`. The arithmetic is the same and
in the same order, and it is compiled without fused multiply-add, as dials is. On all
26,639 strong spots of the 460-frame corpus (detectors of 256, 32, 8 and 1 panels),
every millimetre position, variance, s1 and reciprocal lattice point is bit-identical to
dials. Evidence: `stills_process_profiling/claude/rewrite/runs/r5-plumbing/`.

**Stills refinement** (`src/stills/crystal`, `predict`, `refine`). dials' `e_refine` for one
crystal with the beam and detector fixed:
- the prediction: q = UBh rotated about e1 = q x s0 onto the Ewald sphere by Delta-Psi
  (Sauter et al. 2014, eqs. 4-12), in dials' operation order;
- the residuals X, Y in mm and Delta-Psi, weights 1/var(X), 1/var(Y) and 1e6 (Brewster et al.
  2018, eq. 5), and the reflection filter (an index, not overloaded, not unsuitable);
- U moved by three angles in mrad about the laboratory axes, and B as the lower-triangular
  factor of G* = B^T B. The free elements of G*, scaled by 1e5, are the non-pivot columns of
  the reduced row-echelon form of the group's invariance equations, which is the choice dials
  makes in all six groups the corpus holds (P1, P2_1, P2_12_12_1, P4_32_12, P4/mmm, P6/mmm);
- parameters no spot moves by more than 1e-10 are fixed, as dials' AutoReduce does;
- dials' Levenberg-Marquardt: mu from 1e-3 times the largest diagonal, Nielsen's update, a zero
  gradient element refused, the stop on every RMSD changing by less than 1e-4 relative, 100
  trial steps, a rejected step's journal row deleted, and on the 8192 exit the model left at
  the last rejected trial.

On the 5,186 single-crystal LevMar runs dials makes on the R3 corpus, replayed from their
captured inputs: every outcome is the same, including the five failures (four crystals that
cannot be composed, one empty table); the starting objective and RMSDs agree to 4.4e-12
relative; at dials' final model the predictions agree to 1.25e-12 mm and 6e-16 rad.

The end points do not agree to the last digits, for two measured reasons. dials' rotation-axis
term (`drde`) is a central difference with step 1e-8, whose rounding perturbs the Jacobian by up
to 7e-6 relative on the normal matrix's diagonal; the same arithmetic in another order does not
reproduce that rounding. And dials stops on a 1e-4 relative change in RMSD, in a valley of the
objective flat enough that points whose objectives agree to 1e-8 lie 1e-4 degrees apart. So
end points differ by up to a few 1e-4 degrees (one psii run 1e-2), and the objective at the two
end points, evaluated by one function, agrees to 5.4e-8 relative on every run except the two
in row 1. Evidence: `stills_process_profiling/claude/rewrite/runs/r6-refine/`.

Not yet departures, but named as candidates:
- **The rotation-axis derivative.** dials' `drde` moves the unit axis e1 along the part of
  dq0 x s0 perpendicular to it; the derivative of e1 = normalise(q0 x s0) divides that by
  |q0 x s0|. Against a central difference of the whole prediction, dials' form is off by up to
  4e-4 of the orientation columns, the analytic form by under 1e-8 (tests/stills). With
  `axis_derivative = analytic` the replay meets the 1e-4 degree tolerance on 97.6 % of runs
  (5,059 of 5,186), and every miss but the two in row 1 is in the flat valley; the C++ ends at a
  lower objective than dials on 71-85 % of runs per dataset. The default stays `dials`.
 replacing Nave's mosaicity with a stills
partiality model, and replacing Niggli reduction with Selling reduction if a profile shows that
reduction costs time.

# Stills indexing: where it departs from dials

`mxi_stills` reimplements, in new C++, the algorithms `dials.stills_process` uses to index a
still. Each entry here is a place where it does something different from dials, or from the
paper the algorithm comes from. An entry is written when the change is made, and gives its
reason: a measured cost, with the run that measured it, or a stated weakness of the model.

| # | date | what differs | from (dials code or paper) | reason | measured effect |
|---|---|---|---|---|---|
| 1 | 2026-10-11 | Refinement stops, crystal unchanged, when the equations number fewer than the free parameters (3N < P) | dials' LevMar `dof = 3N - P` is a `size_t` (scitbx `normal_equations.h`), so it wraps and dials iterates on an underdetermined system | The wrap is an unsigned-arithmetic accident, and the result is a model fitted to fewer equations than parameters | 1 LevMar run on each of lysozyme and luciferase in the R6 replay (spec: 16 and 37 on the full sets). Each is the first refinement of an outlier pass with N = 1 or 2 spots, whose Sauter-Poon step then fails for N < 10, so no frame's outcome changes |
| 2 | 2026-10-11 | No detector parameterisation is built, so the panel frames are not recomposed on every refinement | dials builds a hierarchical detector parameterisation on every refiner build, whose `compose()` rewrites the copied panels at about 1e-12 mm and is written back to the experiment (spec 9.2) | The detector is fixed on every configuration (STATUS 3, decided 2026-10-10); the parameterisation exists only to be discarded | None resolvable here: on 5,186 captured runs the starting objective and RMSDs agree with dials to 4.4e-12 relative |
| 3 | 2026-10-11 | d-spacings, 2-theta and the cell volume for the Nave fit come from B: d = 1/\|Bh\|, V = 1/det B | dials takes them from cctbx's unit cell, rebuilt from the crystal's parameters (`crystal.get_unit_cell().d()`, `.two_theta()`, `.volume()`) | No cctbx in the C++; the unit-cell object exists only to answer these three questions | On the 2,403 single-crystal Nave fits of the R7 replay: d within 7.3e-16 relative and the volume within 7.0e-16; every outcome and every acceptance flag the same as dials'. The 2-theta order (ties between reflections of equal spacing) differs in 59 fits, none of them among the 13 below. alpha and eta agree within 1e-6 relative on 2,204 of the 2,217 fits that complete; the other 13 take the same L-BFGS path, with the same iteration and call counts, to points up to 4.5e-4 apart (eta, one psii fit) whose objectives agree to 3.8e-9 relative: a flat direction of the likelihood, the same class as the refinement replay's flat valley |
| 4 | 2026-10-11 | Where serialtbx's Nave likelihood would overflow (`exp(B (psi +- psi_model)) > 1e100`), the same likelihood and gradient are evaluated in log-sum-exp form (`NaveOptions::safe_likelihood`, default on) | serialtbx `max_like.py:117-122` raises "Model has diverged", which fails the candidate, or the frame after indexing | The model does not diverge, the evaluation does: on the corpus all 25 such fits fail at a line-search trial point 3-6 evaluations in, after a first step from a mis-sized seed to a band 23-64 times too narrow; the exception stops the search from stepping back. dials' gradient is the true derivative (central differences, 5.5e-7), so the form is the same function | Replay: every fit that completes in dials completes bit-identically; all 25 converge (5-10 iterations) and pass the Ewald limit. End to end (`runs/r7-nave-fixes/`, against pr23_B): +6 cytochrome and +22 lysozyme frames integrated, none lost, 3 frames on another lattice; merged CC1/2 +0.0013 +- 0.0013 (cytochrome) and +0.0013 +- 0.0011 (lysozyme), paired over six CC1/2 splits; merged <I/sigma> +0.008 and +0.003 |

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

**Sauter-Poon outlier rejection** (`src/stills/outliers`, with `lbfgs`). rstbx's
`find_outliers_from_matches` on a first pass, as dials' stills indexer calls it (Sauter & Poon
2010): radial distances between observed and calculated positions in millimetres (each pixel
position times panel 0's pixel size, as dials forms them), two equal distances an error, the
empirical CDF k/N, the Rayleigh fit on the smallest int(0.4 N) from the last of them, its L1
objective with the L2 gradient, scitbx's L-BFGS with scitbx.lbfgs.run's termination and
exception handling, the 0.655 sigma bar, and every spot from the first crossing outward
rejected. The calculated pixel position is dxtbx's `millimeter_to_pixel` with the forward
parallax correction (`StillsPanel::mm_to_px`).

**The Nave fit** (`src/stills/nave`). dials' `NaveParameters` for one crystal (Sauter et al.
2014, eq. 16 and appendix B): the binned seed in 2-theta order, the 2x2 solve as scitbx.matrix
forms it, the 3 V^(1/3) floor on the domain size, serialtbx's maximum-likelihood objective and
its coded gradient, L-BFGS with gtol 0.1 and the drop test, the acceptance band, the 1.4
widening and the Ewald proximal volume reject at 0.0025. The L-BFGS is scitbx's (Liu & Nocedal
1989; More & Thuente 1994) followed line by line, because the fits stop on tests that a
different minimiser would meet at different points.

On every Sauter-Poon pass and Nave fit dials makes on the R3 corpus (2,728 passes and 2,403
fits on single crystals), replayed from their captured inputs: every outcome is the same,
including dials' 100 too-few rejections, 25 diverged fits and 161 Ewald-volume rejections;
every kept and accepted flag is the same; the fitted sigma, the seed, alpha, eta, the crystal's
half-mosaicity and domain size and the Ewald volume are bit-identical, and the minimiser takes
the same number of iterations and function calls. The calculated pixel positions are
bit-identical to dials' on every reflection. The minimiser alone matches scitbx's final point,
iteration and call counts on the Rosenbrock function to the last bit (tests/stills). Evidence:
`stills_process_profiling/claude/rewrite/runs/r7-outliers-nave/`.

**Candidate crystal models** (`src/stills/candidates`, `reduction`, `cell`). dials'
`candidate_orientation_matrices` (Gildea et al. 2014): triplets of the first 100 basis vectors
in combination order, stable-sorted by i^2 + j^2 + k^2; the 20 degree tests on angle(a, b) and
on a x b against c; b and c flipped to make the angles obtuse, all three negated for a
left-handed triplet; dxtbx's `Crystal(a, b, c)`; Niggli reduction; dxtbx's `change_basis`; and
the V > abc/100 cut. The Niggli reduction is cctbx's Krivy-Gruber with the epsilon comparisons
of Grosse-Kunstleve, Sauter & Adams (2004): epsilon = V^(1/3) 1e-5 against the squared
quantities, the same step order, the last epsilon-zero element taking the sign in A4, exact
sign tests in A5-A7, and the 1000-step limit. The unit cell and the crystal are computed as
cctbx's `unit_cell` and dxtbx's `Crystal` compute them, through the same cell parameters,
metric and fractionalisation matrix.

**Global index assignment** (`src/stills/candidates`). dials' `AssignIndicesGlobal` for one
image: the strict d > d_min cut on unindexed rows, h = A^-1 rlp rounded half away from zero,
the first crystal with the smallest squared distance to its integer, the 0.3 tolerance, (0, 0, 0)
refused, and of rows given the same h on the same crystal the nearer kept, ties to the earlier
row (Gildea et al. 2014, section 3.1). The rotation-angle test between duplicates never applies
to a still and is not implemented.

On every candidate search and index assignment dials makes on the R3 corpus (491 searches, 1,423
assignments), replayed from their captured inputs: the 547,597 candidate models are the same,
in the same order, with A and the cell parameters bit-identical, and every miller index and
crystal id is the same. Evidence: `stills_process_profiling/claude/rewrite/runs/r8a-candidates/`.

Not yet departures, but named as candidates:
- **The rotation-axis derivative.** dials' `drde` moves the unit axis e1 along the part of
  dq0 x s0 perpendicular to it; the derivative of e1 = normalise(q0 x s0) divides that by
  |q0 x s0|. Against a central difference of the whole prediction, dials' form is off by up to
  4e-4 of the orientation columns, the analytic form by under 1e-8 (tests/stills). With
  `axis_derivative = analytic` the replay meets the 1e-4 degree tolerance on 97.6 % of runs
  (5,060 of 5,186), and every miss but the two in row 1 is in the flat valley; the C++ ends at a
  lower objective than dials on 71-85 % of runs per dataset. The default stays `dials`.
- **The Sauter-Poon bar.** The paper rejects a spot more than sigma beyond the fitted curve
  (eq. 6, Fig. 3c); rstbx uses the Rayleigh standard deviation, 0.655 sigma, and rejects every
  spot from the first crossing outward, where eq. 6 classifies each spot on its own.
- **Measured and not adopted** (`runs/r7-nave-fixes/`, R7-FIXES-PLAN.md): no Ewald proximal volume
  limit (more frames index, but merged CC1/2 falls, -0.0047 +- 0.0042 on cytochrome, and the
  added crystals carry ~1 degree half-mosaicity); and the integration predictor at eta/2 below
  (completeness +6 to +8 points and multiplicity x1.3-1.6, but lysozyme CC1/2 -0.0091 +- 0.0057
  and 22-65 % more reflections to integrate). Both wait for the replacement model.
- **The integration predictor's mosaic term.** dials' `NaveStillsReflectionPredictor` adds a
  quarter of the fitted full mosaicity where eq. 16 adds a half, so integration accepts a
  narrower band than the fit describes. Not in this library yet; it matters when integration
  is.
- **Replacing Nave's mosaicity** with a stills partiality model. Its "mosaicity" and domain size
  are fitted parameters with a loose physical meaning, the 1.4 widening and the Ewald-volume
  limit are empirical, and a replacement is to be scored by merged statistics, not by agreement
  with Nave.
- **Selling reduction** in place of Niggli reduction, if a profile shows that reduction costs
  time. Measured so far (R8a replay): the C++ makes a search's whole candidate list, reductions
  included, in 0.14-0.16 of the time dials takes.

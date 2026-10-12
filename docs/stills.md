# Stills indexing for dials.stills_process

`mxi_stills` is a Python module. `dials.stills_process` calls it in place of
its own indexing step, which covers:

- the basis-vector search;
- candidate refinement and outlier rejection;
- the post-index refinement.

The code is C++, in a library of its own (`src/stills/`) built over the mxi
library. The module wraps it with pybind11 (`bindings/`).

## What it does so far

- **The frame model** (`src/stills/geometry.hh`). A multi-panel detector, with
  each panel's laboratory-frame vectors, pixel size and parallax coefficients.
  It is built once per process and kept. `map_spots` takes a still's spot
  centroids to millimetres, s1 and reciprocal lattice points, with dials'
  arithmetic in dials' order. The corpus check finds it bit-identical to dials
  on every spot.
- **Stills refinement** (`src/stills/refine.hh`, with `crystal.hh` and
  `predict.hh`). `refine_crystal` refines one still's crystal, its orientation
  and its cell constrained by the space group, against its indexed spots, with
  the beam and detector fixed. Prediction, target, weights and minimiser are
  dials' stills refinement (`docs/stills-review.md` lists where they differ).
  `predict_spots` and `linearise_crystal` expose the prediction and the
  least-squares problem at the starting model. `AxisDerivative` chooses the
  form of the rotation-axis term of the Jacobian: `dials` (the default),
  `dials_limit`, or `analytic`.
- **Outlier rejection and mosaicity** (`src/stills/outliers.hh`,
  `src/stills/nave.hh`, with `lbfgs.hh`). `sauter_poon` rejects spots by the
  Rayleigh fit of their distances from prediction, and `nave_parameters` fits a
  crystal's mosaicity and domain size to its Delta-Psi and returns the
  reflections inside the fitted band, as dials' stills indexer does both.
  By default its likelihood does not overflow where dials' raises "Model has
  diverged" (`NaveOptions::safe_likelihood`, ledger row 4); off, the fit is
  dials' to the last bit.
  `lbfgs::run` is scitbx's L-BFGS minimiser, which both fits use.
  `StillsPanel::mm_to_px` is dxtbx's millimetre-to-pixel conversion with the
  parallax correction.
- **Candidates and index assignment** (`src/stills/candidates.hh`, with
  `reduction.hh` and `cell.hh`). `candidate_orientation_matrices` makes dials'
  candidate crystal models from a search's basis vectors, Niggli-reduced as
  cctbx reduces them (`niggli_reduction`), and `assign_indices_global` indexes
  a still's spots against one or more crystals as dials'
  `AssignIndicesGlobal` does. `UnitCell` and `CrystalModel` compute a cell and
  a crystal as cctbx's `unit_cell` and dxtbx's `Crystal` do.
- **Known symmetry** (`src/stills/known_symmetry.hh`, with `sgtbx.hh`).
  `KnownSymmetry` and `filter_known_symmetry` keep the candidates dials keeps
  when the space group and cell are known, through `find_matching_symmetry`:
  the minimum cell, the lattice group from the 81 two-folds, its subgroups of
  the target Bravais type, their reference settings and best cells, and the
  obliquity of each. `sgtbx.hh` is the part of cctbx's sgtbx this needs:
  rational symmetry operations, space groups kept in cctbx's operation order,
  the Hall symbol parser, and the space-group type with its change of basis to
  the reference setting. Its reference tables (`sgtbx_tables.cc`) are copied
  from cctbx by a script and are not edited by hand.
- **Applying the target symmetry** (`src/stills/apply_symmetry.hh`).
  `SymmetryHandler` holds the user's target as dials' handler does, and its
  `apply_symmetry` puts a candidate into the target's lattice symmetry: the
  candidate's metric subgroups (`iotbx_converter`), the best of the target's
  Bravais type, the orientation constrained to that crystal system
  (`constrain_orientation`, rstbx's L-BFGS fit), and the change of basis to
  the user's setting.
- **The fft1d search** (`src/stills/fft1d.hh`, with `fftpack.hh`).
  `fft1d_basis_vectors` is dials' `FFT1D.find_basis_vectors`: rstbx's DPS
  search over the hemisphere, with its adaptive two-pass grid, refinement and
  ranking, and scitbx's FFTPACK transform.
- **The method fallback** (`src/stills/method_list.hh`).
  `index_with_method_list` tries `indexing.stills.method_list` in order, as
  stills_process does.
- **The entry point** (`src/stills/index.hh`). `index_still` maps the spots and
  returns "not indexed", because there is no search yet.
- **The Python side** (`bindings/mxi_stills/dxtbx_adapter.py`). It builds the
  arrays from live dxtbx models and keeps the C++ detector between stills. Panel
  vectors come from the live panel, never from an `.expt`, which stores them in
  the panel's local frame. Parallax coefficients come from the panel's px/mm
  strategy.
- **In dials.stills_process**, `indexing.stills.backend=mxi` selects it.

The tests in `tests/stills/` (`mxi_stills_tests`) check the geometry against
dials' own values for real spots on detectors of 256, 32, 8 and 1 panels, and
the refinement against dials' own runs on real stills in five space groups,
the outlier rejection and Nave fit against dials' own calls, and the cell
reduction, candidates and index assignment, the space-group types, the
known-symmetry filter, apply_symmetry and the fft1d search against cctbx's and
dials' own results, and the method fallback against stub searches. Those values are
in `tests/stills/real_stills.hh`, `tests/stills/real_refine.hh`,
`tests/stills/real_outliers.hh`, `tests/stills/real_candidates.hh` and
`tests/stills/real_symmetry.hh` and `tests/stills/real_fft1d.hh`, which are
generated.

## Building

The library is built by every plain CMake build. Only the Python module
needs a Python environment. `pip` builds it from the root `pyproject.toml`
with scikit-build-core, which turns on `MXI_PYTHON` and builds only that
target:

```sh
git submodule update --init --recursive
CXX=<the C++ compiler the target Python was built with> \
    <target python> -m pip install .
```

The build tree goes in `build-py/` and is kept between installs, so a rebuild
after an edit is incremental. Use the same compiler as the Python the module
is loaded into, so the C++ runtime matches. In a DIALS conda environment that
is `conda_base/bin/x86_64-conda-linux-gnu-c++`.

To check what is installed:

```sh
<target python> -c "import mxi_stills; print(mxi_stills.git_revision)"
<target python> -m pytest bindings/tests
```

`git_revision` is `git describe` at configure time, and ends in `-dirty` if
the tree had uncommitted changes. The module is installed with pip, so no
source checkout records which build is loaded. This attribute is the only
place to read it.

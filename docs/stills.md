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
- **The entry point** (`src/stills/index.hh`). `index_still` maps the spots and
  returns "not indexed", because there is no search yet.
- **The Python side** (`bindings/mxi_stills/dxtbx_adapter.py`). It builds the
  arrays from live dxtbx models and keeps the C++ detector between stills. Panel
  vectors come from the live panel, never from an `.expt`, which stores them in
  the panel's local frame. Parallax coefficients come from the panel's px/mm
  strategy.
- **In dials.stills_process**, `indexing.stills.backend=mxi` selects it.

The tests in `tests/stills/` (`mxi_stills_tests`) check the geometry against
dials' own values for real spots on detectors of 256, 32, 8 and 1 panels.
Those values are in `tests/stills/real_stills.hh`, which is generated.

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

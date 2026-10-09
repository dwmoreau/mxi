# Stills indexing for dials.stills_process

`mxi_stills` is a Python module. `dials.stills_process` calls it in place of
its own indexing step, which covers:

- the basis-vector search;
- candidate refinement and outlier rejection;
- the post-index refinement.

The code is C++, in a library of its own (`src/stills/`) built over the mxi
library. The module wraps it with pybind11 (`bindings/`).

So far it only reports which build it is. The indexer itself comes next.

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

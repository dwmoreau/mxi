// The Python module dials.stills_process calls the stills indexer through.

#include <pybind11/pybind11.h>

#include "stills/build_info.hh"

namespace py = pybind11;

PYBIND11_MODULE(_mxi_stills, m) {
  m.doc() = "mxi stills indexing for dials.stills_process";

  const auto info = mxi::stills::build_info();
  m.attr("__version__") = info.version;
  m.attr("git_revision") = info.git_revision;
  m.attr("compiler") = info.compiler;
}

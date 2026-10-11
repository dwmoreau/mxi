// The Python module dials.stills_process calls the stills indexer through.
// Plain numpy arrays cross the boundary in both directions.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>

#include "stills/build_info.hh"
#include "stills/geometry.hh"
#include "stills/index.hh"

namespace py = pybind11;
using namespace mxi::stills;

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IndexArray =
    py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

void require_shape(const py::array &a, py::ssize_t rows, py::ssize_t cols,
                   const char *name) {
  const bool ok =
      cols == 0 ? (a.ndim() == 1 && a.shape(0) == rows)
                : (a.ndim() == 2 && a.shape(0) == rows && a.shape(1) == cols);
  if (!ok)
    throw py::value_error(std::string(name) + ": expected shape (" +
                          std::to_string(rows) +
                          (cols ? ", " + std::to_string(cols) : "") + ")");
}

mxi::Vec3 vec3(const Array &a, const char *name) {
  require_shape(a, 3, 0, name);
  const double *p = a.data();
  return {p[0], p[1], p[2]};
}

StillsDetector make_detector(const Array &origin, const Array &fast,
                             const Array &slow, const Array &pixel_size,
                             const Array &mu, const Array &t0) {
  const py::ssize_t n = origin.ndim() == 2 ? origin.shape(0) : -1;
  require_shape(origin, n, 3, "origin");
  require_shape(fast, n, 3, "fast");
  require_shape(slow, n, 3, "slow");
  require_shape(pixel_size, n, 2, "pixel_size");
  require_shape(mu, n, 0, "mu");
  require_shape(t0, n, 0, "t0");
  StillsDetector d;
  d.panels.resize(static_cast<std::size_t>(n));
  for (py::ssize_t i = 0; i < n; ++i) {
    StillsPanel &p = d.panels[static_cast<std::size_t>(i)];
    p.origin = {origin.at(i, 0), origin.at(i, 1), origin.at(i, 2)};
    p.fast = {fast.at(i, 0), fast.at(i, 1), fast.at(i, 2)};
    p.slow = {slow.at(i, 0), slow.at(i, 1), slow.at(i, 2)};
    p.pixel_size[0] = pixel_size.at(i, 0);
    p.pixel_size[1] = pixel_size.at(i, 1);
    p.mu = mu.at(i);
    p.t0 = t0.at(i);
  }
  return d;
}

// Spots from the reflection table's panel, xyzobs.px.value and
// xyzobs.px.variance columns; z is ignored.
Spots make_spots(const IndexArray &panel, const Array &xyz_px,
                 const Array &var_px) {
  const py::ssize_t n = panel.ndim() == 1 ? panel.shape(0) : -1;
  require_shape(panel, n, 0, "panel");
  require_shape(xyz_px, n, 3, "xyzobs.px.value");
  require_shape(var_px, n, 3, "xyzobs.px.variance");
  Spots s;
  s.panel.assign(panel.data(), panel.data() + n);
  s.x_px.resize(static_cast<std::size_t>(n));
  s.y_px.resize(static_cast<std::size_t>(n));
  s.var_x_px.resize(static_cast<std::size_t>(n));
  s.var_y_px.resize(static_cast<std::size_t>(n));
  for (py::ssize_t i = 0; i < n; ++i) {
    const auto k = static_cast<std::size_t>(i);
    s.x_px[k] = xyz_px.at(i, 0);
    s.y_px[k] = xyz_px.at(i, 1);
    s.var_x_px[k] = var_px.at(i, 0);
    s.var_y_px[k] = var_px.at(i, 1);
  }
  return s;
}

Array to_array(const std::vector<mxi::Vec3> &v) {
  Array a({static_cast<py::ssize_t>(v.size()), py::ssize_t{3}});
  auto r = a.mutable_unchecked<2>();
  for (std::size_t i = 0; i < v.size(); ++i) {
    const auto k = static_cast<py::ssize_t>(i);
    r(k, 0) = v[i].x;
    r(k, 1) = v[i].y;
    r(k, 2) = v[i].z;
  }
  return a;
}

// Two columns and a zero third, as dials stores a still's mm position.
Array to_array(const std::vector<double> &x, const std::vector<double> &y) {
  Array a({static_cast<py::ssize_t>(x.size()), py::ssize_t{3}});
  auto r = a.mutable_unchecked<2>();
  for (std::size_t i = 0; i < x.size(); ++i) {
    const auto k = static_cast<py::ssize_t>(i);
    r(k, 0) = x[i];
    r(k, 1) = y[i];
    r(k, 2) = 0.0;
  }
  return a;
}

} // namespace

PYBIND11_MODULE(_mxi_stills, m) {
  m.doc() = "mxi stills indexing for dials.stills_process";

  const auto info = build_info();
  m.attr("__version__") = info.version;
  m.attr("git_revision") = info.git_revision;
  m.attr("compiler") = info.compiler;

  py::register_exception<GeometryError>(m, "GeometryError", PyExc_ValueError);

  py::class_<StillsDetector>(m, "Detector")
      .def(py::init(&make_detector), py::arg("origin"), py::arg("fast"),
           py::arg("slow"), py::arg("pixel_size"), py::arg("mu"), py::arg("t0"),
           "Panels in the laboratory frame: origin, fast and slow (n, 3), "
           "pixel_size (n, 2), and the parallax mu and t0 (n), zero for none.")
      .def("__len__", &StillsDetector::size);

  m.def(
      "map_spots",
      [](const StillsDetector &detector, const Array &s0, double wavelength,
         const IndexArray &panel, const Array &xyz_px, const Array &var_px) {
        const Spots spots = make_spots(panel, xyz_px, var_px);
        const MappedSpots mapped =
            map_spots(detector, vec3(s0, "s0"), wavelength, spots);
        py::dict out;
        out["xyzobs.mm.value"] = to_array(mapped.x_mm, mapped.y_mm);
        out["xyzobs.mm.variance"] = to_array(mapped.var_x_mm, mapped.var_y_mm);
        out["s1"] = to_array(mapped.s1);
        out["rlp"] = to_array(mapped.rlp);
        return out;
      },
      py::arg("detector"), py::arg("s0"), py::arg("wavelength"),
      py::arg("panel"), py::arg("xyz_px"), py::arg("var_px"),
      "Map spot centroids to millimetres, s1 and the reciprocal lattice "
      "point, as dials' centroid_px_to_mm and "
      "map_centroids_to_reciprocal_space do for a still.");

  py::class_<IndexResult>(m, "IndexResult")
      .def_property_readonly(
          "indexed",
          [](const IndexResult &r) { return r.status == IndexStatus::indexed; })
      .def_readonly("reason", &IndexResult::reason)
      .def_readonly("n_spots", &IndexResult::n_spots);

  m.def(
      "index_still",
      [](const StillsDetector &detector, const Array &s0, double wavelength,
         const IndexArray &panel, const Array &xyz_px, const Array &var_px) {
        const Spots spots = make_spots(panel, xyz_px, var_px);
        return index_still(detector, vec3(s0, "s0"), wavelength, spots);
      },
      py::arg("detector"), py::arg("s0"), py::arg("wavelength"),
      py::arg("panel"), py::arg("xyz_px"), py::arg("var_px"));
}

// Python bindings for the basis-vector searches.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <optional>

#include "stills/fft1d.hh"
#include "stills/rsgs.hh"

namespace py = pybind11;
using namespace mxi::stills;

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

std::vector<mxi::Vec3> vectors_from(const Array &a, const char *name) {
  if (a.ndim() != 2 || a.shape(1) != 3)
    throw py::value_error(std::string(name) + ": expected shape (n, 3)");
  std::vector<mxi::Vec3> v(static_cast<std::size_t>(a.shape(0)));
  for (py::ssize_t i = 0; i < a.shape(0); ++i)
    v[static_cast<std::size_t>(i)] = {a.at(i, 0), a.at(i, 1), a.at(i, 2)};
  return v;
}

Array vectors_array(const std::vector<mxi::Vec3> &v) {
  Array out({static_cast<py::ssize_t>(v.size()), py::ssize_t{3}});
  auto o = out.mutable_unchecked<2>();
  for (std::size_t i = 0; i < v.size(); ++i)
    for (py::ssize_t c = 0; c < 3; ++c)
      o(static_cast<py::ssize_t>(i), c) = v[i][static_cast<std::size_t>(c)];
  return out;
}

Array rows_array(const std::vector<std::array<double, 5>> &rows) {
  Array out({static_cast<py::ssize_t>(rows.size()), py::ssize_t{5}});
  auto o = out.mutable_unchecked<2>();
  for (std::size_t i = 0; i < rows.size(); ++i)
    for (std::size_t c = 0; c < 5; ++c)
      o(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(c)) = rows[i][c];
  return out;
}

} // namespace

void bind_basis_vector_search(py::module_ &m) {
  py::register_exception<TooFewSpots>(m, "TooFewSpots", PyExc_ValueError);

  m.def(
      "fft1d_basis_vectors",
      [](const Array &rlp, double max_cell,
         std::optional<double> characteristic_grid, bool trace, double max_grid,
         double quick_grid, double refine_divisor) -> py::object {
        Fft1dTrace t;
        const auto v = fft1d_basis_vectors(
            vectors_from(rlp, "rlp"), max_cell, characteristic_grid,
            trace ? &t : nullptr, {max_grid, quick_grid, refine_divisor});
        if (!trace)
          return vectors_array(v);
        py::dict d;
        d["characteristic_grid"] = t.characteristic_grid;
        d["incr"] = t.incr;
        d["adaptive"] = t.adaptive;
        d["n_quick"] = t.n_quick;
        d["n_quick_winners"] = t.n_quick_winners;
        d["n_directions"] = t.n_directions;
        d["seconds_quick"] = t.seconds_quick;
        d["seconds_scan"] = t.seconds_scan;
        d["seconds_refine"] = t.seconds_refine;
        d["unrefined"] = rows_array(t.unrefined);
        d["refined"] = rows_array(t.refined);
        return py::make_tuple(vectors_array(v), d);
      },
      py::arg("rlp"), py::arg("max_cell"),
      py::arg("characteristic_grid") = py::none(), py::arg("trace") = false,
      py::arg("max_grid") = Fft1dGrid{}.max_grid,
      py::arg("quick_grid") = Fft1dGrid{}.quick_grid,
      py::arg("refine_divisor") = Fft1dGrid{}.refine_divisor,
      "dials' FFT1D.find_basis_vectors: the basis vectors (n, 3), strongest "
      "first. With trace, also the intermediate steps: the grid, the "
      "directions searched, (kval, length, direction) before and after "
      "refinement, and the seconds of each stage. max_grid, quick_grid and "
      "refine_divisor set the grid; their defaults are rstbx's.");

  py::class_<RealSpaceGridSearch>(m, "RealSpaceGridSearch")
      .def(py::init([](const std::array<double, 6> &cell,
                       double characteristic_grid) {
             return RealSpaceGridSearch(cell, characteristic_grid);
           }),
           py::arg("target_unit_cell"), py::arg("characteristic_grid") = 0.02,
           "dials' RealSpaceGridSearch for a target cell (six parameters).")
      .def(
          "find_basis_vectors",
          [](const RealSpaceGridSearch &s, const Array &rlp,
             std::size_t max_vectors) {
            return vectors_array(
                s.find_basis_vectors(vectors_from(rlp, "rlp"), max_vectors));
          },
          py::arg("rlp"), py::arg("max_vectors") = 30,
          "The basis vectors (n, 3), strongest first.")
      .def(
          "scores",
          [](const RealSpaceGridSearch &s, const Array &rlp) {
            const auto w = s.scores(vectors_from(rlp, "rlp"));
            return py::array_t<double>(static_cast<py::ssize_t>(w.size()),
                                       w.data());
          },
          py::arg("rlp"), "One score per search vector, in their order.")
      .def_property_readonly("search_vectors",
                             [](const RealSpaceGridSearch &s) {
                               return vectors_array(s.search_vectors());
                             });

  m.def(
      "hemisphere_grid",
      [](double sampling) { return vectors_array(hemisphere_grid(sampling)); },
      py::arg("sampling"),
      "rstbx's construct_hemisphere_grid: the directions, in its order.");
}

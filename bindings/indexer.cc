// Python bindings for the stills indexer.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "stills/indexer.hh"

namespace py = pybind11;
using namespace mxi::stills;

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IndexArray =
    py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;
using FlagArray =
    py::array_t<std::uint64_t, py::array::c_style | py::array::forcecast>;

void require_rows(const py::array &a, py::ssize_t rows, py::ssize_t cols,
                  const char *name) {
  const bool ok =
      cols == 0 ? (a.ndim() == 1 && a.shape(0) == rows)
                : (a.ndim() == 2 && a.shape(0) == rows && a.shape(1) == cols);
  if (!ok)
    throw py::value_error(std::string(name) + ": expected shape (" +
                          std::to_string(rows) +
                          (cols ? ", " + std::to_string(cols) : "") + ")");
}

Array xyz_array(const std::vector<double> &x, const std::vector<double> &y) {
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

Array vec_array(const std::vector<mxi::Vec3> &v) {
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

Array mat_array(const mxi::Mat3 &m) {
  Array a({py::ssize_t{3}, py::ssize_t{3}});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < 3; ++i)
    for (py::ssize_t j = 0; j < 3; ++j)
      r(i, j) = m(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
  return a;
}

AxisDerivative axis_derivative_of(const std::string &name) {
  if (name == "dials")
    return AxisDerivative::dials;
  if (name == "dials_limit")
    return AxisDerivative::dials_limit;
  if (name == "analytic")
    return AxisDerivative::analytic;
  throw py::value_error(
      "axis_derivative: expected 'dials', 'dials_limit' or 'analytic'");
}

py::list stage_list(const std::vector<StageTrace> &stages) {
  py::list out;
  for (const StageTrace &t : stages)
    out.append(py::make_tuple(t.step, t.n_in, t.n_kept, t.a, t.b));
  return out;
}

py::dict to_dict(const IndexedStill &r) {
  py::dict out;
  out["indexed"] = r.indexed;
  out["reason"] = r.reason;
  out["method"] = r.method;

  py::list methods;
  for (const MethodTrace &m : r.methods) {
    py::dict d;
    d["method"] = m.method;
    d["error"] = m.error;
    d["n_basis_vectors"] = m.n_basis_vectors;
    d["d_min"] = m.d_min;
    d["best"] = m.best;
    d["seconds_search"] = m.seconds_search;
    d["seconds_candidates"] = m.seconds_candidates;
    d["seconds_refine"] = m.seconds_refine;
    py::list candidates;
    for (const CandidateTrace &c : m.candidates) {
      py::dict cd;
      cd["status"] = c.status;
      cd["n_indexed"] = c.n_indexed;
      cd["rmsd"] = c.rmsd;
      cd["cell"] = py::make_tuple(c.cell[0], c.cell[1], c.cell[2], c.cell[3],
                                  c.cell[4], c.cell[5]);
      cd["half_mosaicity_deg"] = c.half_mosaicity_deg;
      cd["domain_size_ang"] = c.domain_size_ang;
      cd["stages"] = stage_list(c.stages);
      candidates.append(cd);
    }
    d["candidates"] = candidates;
    d["post_index"] = stage_list(m.post_index);
    d["chosen_A"] = mat_array(m.chosen_A);
    methods.append(d);
  }
  out["methods"] = methods;

  out["xyzobs.mm.value"] = xyz_array(r.mapped.x_mm, r.mapped.y_mm);
  out["xyzobs.mm.variance"] = xyz_array(r.mapped.var_x_mm, r.mapped.var_y_mm);
  out["rlp"] = vec_array(r.mapped.rlp);
  if (!r.indexed)
    return out;

  out["U"] = mat_array(r.crystal.U);
  out["B"] = mat_array(r.crystal.B);
  out["B_raw"] = mat_array(r.raw_B);
  const auto &p = r.crystal.cell.params;
  out["unit_cell"] = py::make_tuple(p[0], p[1], p[2], p[3], p[4], p[5]);
  out["space_group_ops"] = r.space_group_ops;
  out["half_mosaicity_deg"] = r.half_mosaicity_deg;
  out["domain_size_ang"] = r.domain_size_ang;
  out["ml_half_mosaicity_deg"] = r.ml_half_mosaicity_deg;
  out["ml_domain_size_ang"] = r.ml_domain_size_ang;

  const auto n = static_cast<py::ssize_t>(r.row.size());
  IndexArray row(n);
  IndexArray miller({n, py::ssize_t{3}});
  Array delpsi(n);
  FlagArray flags(n);
  for (py::ssize_t i = 0; i < n; ++i) {
    const auto k = static_cast<std::size_t>(i);
    row.mutable_at(i) = static_cast<std::int64_t>(r.row[k]);
    for (py::ssize_t c = 0; c < 3; ++c)
      miller.mutable_at(i, c) = r.miller_index[k][static_cast<std::size_t>(c)];
    delpsi.mutable_at(i) = r.delpsi[k];
    flags.mutable_at(i) = r.flags[k];
  }
  out["row"] = row;
  out["miller_index"] = miller;
  out["xyzcal.mm"] = xyz_array(r.x_cal_mm, r.y_cal_mm);
  out["xyzcal.px"] = xyz_array(r.x_cal_px, r.y_cal_px);
  out["delpsical.rad"] = delpsi;
  out["s1"] = vec_array(r.s1);
  out["flags"] = flags;
  return out;
}

} // namespace

void bind_indexer(py::module_ &m) {
  py::class_<StillsIndexer>(m, "StillsIndexer")
      .def(py::init([](const StillsDetector &detector,
                       const std::array<double, 6> &unit_cell,
                       const std::string &space_group_hall,
                       const std::vector<std::string> &methods,
                       bool refine_candidates_with_known_symmetry,
                       std::optional<double> d_min_start,
                       std::size_t max_refine, double xy_rmsd_threshold,
                       std::optional<double> fft1d_characteristic_grid,
                       double rsgs_characteristic_grid,
                       const std::string &axis_derivative, bool safe_likelihood,
                       double ewald_proximal_volume_max,
                       double ewald_proximity_resolution_cutoff) {
             IndexerSettings s;
             s.unit_cell = unit_cell;
             s.space_group_hall = space_group_hall;
             s.methods = methods;
             s.refine_candidates_with_known_symmetry =
                 refine_candidates_with_known_symmetry;
             s.d_min_start = d_min_start.value_or(0.0);
             s.max_refine = max_refine;
             s.xy_rmsd_threshold = xy_rmsd_threshold;
             s.fft1d_characteristic_grid = fft1d_characteristic_grid;
             s.rsgs_characteristic_grid = rsgs_characteristic_grid;
             s.refine.axis_derivative = axis_derivative_of(axis_derivative);
             s.nave.safe_likelihood = safe_likelihood;
             s.nave.ewald_proximal_volume_max = ewald_proximal_volume_max;
             s.nave.ewald_proximity_resolution_cutoff =
                 ewald_proximity_resolution_cutoff;
             return StillsIndexer(detector, s);
           }),
           py::arg("detector"), py::arg("unit_cell"),
           py::arg("space_group_hall"),
           py::arg("methods") = std::vector<std::string>{"fft1d"},
           py::arg("refine_candidates_with_known_symmetry") = false,
           py::arg("d_min_start") = py::none(), py::arg("max_refine") = 5,
           py::arg("xy_rmsd_threshold") = 2.0,
           py::arg("fft1d_characteristic_grid") = py::none(),
           py::arg("rsgs_characteristic_grid") = 0.02,
           py::arg("axis_derivative") = "dials",
           py::arg("safe_likelihood") = true,
           py::arg("ewald_proximal_volume_max") = 0.0025,
           py::arg("ewald_proximity_resolution_cutoff") = 2.0,
           "dials' stills indexer for one detector and one target symmetry "
           "(the cell as (a, b, c, alpha, beta, gamma), the space group as a "
           "Hall symbol), kept for every still a process handles.")
      .def_property_readonly("max_cell", &StillsIndexer::max_cell)
      .def(
          "index",
          [](const StillsIndexer &self, const Array &s0, double wavelength,
             const IndexArray &panel, const Array &xyz_px, const Array &var_px,
             const FlagArray &flags) {
            const py::ssize_t n = panel.ndim() == 1 ? panel.shape(0) : -1;
            require_rows(panel, n, 0, "panel");
            require_rows(xyz_px, n, 3, "xyzobs.px.value");
            require_rows(var_px, n, 3, "xyzobs.px.variance");
            require_rows(flags, n, 0, "flags");
            require_rows(s0, 3, 0, "s0");
            StrongSpots strong;
            const auto k = static_cast<std::size_t>(n);
            strong.spots.panel.assign(panel.data(), panel.data() + n);
            strong.flags.assign(flags.data(), flags.data() + n);
            strong.spots.x_px.resize(k);
            strong.spots.y_px.resize(k);
            strong.spots.var_x_px.resize(k);
            strong.spots.var_y_px.resize(k);
            for (py::ssize_t i = 0; i < n; ++i) {
              const auto j = static_cast<std::size_t>(i);
              strong.spots.x_px[j] = xyz_px.at(i, 0);
              strong.spots.y_px[j] = xyz_px.at(i, 1);
              strong.spots.var_x_px[j] = var_px.at(i, 0);
              strong.spots.var_y_px[j] = var_px.at(i, 1);
            }
            const mxi::Vec3 s{s0.at(0), s0.at(1), s0.at(2)};
            IndexedStill r;
            {
              py::gil_scoped_release release;
              r = self.index(s, wavelength, strong);
            }
            return to_dict(r);
          },
          py::arg("s0"), py::arg("wavelength"), py::arg("panel"),
          py::arg("xyz_px"), py::arg("var_px"), py::arg("flags"),
          "Index one still from its strong spots: a dict with indexed, "
          "reason, the trace of each method tried, the spots mapped to "
          "millimetres and reciprocal space (every row), and when indexed the "
          "crystal and the refined reflections (row: each one's row in the "
          "input).");
}

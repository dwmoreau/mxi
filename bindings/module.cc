// The Python module dials.stills_process calls the stills indexer through.
// Plain numpy arrays cross the boundary in both directions.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <string>

#include "stills/build_info.hh"
#include "stills/geometry.hh"
#include "stills/index.hh"
#include "stills/nave.hh"
#include "stills/outliers.hh"
#include "stills/refine.hh"

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

mxi::Mat3 mat3(const Array &a, const char *name) {
  require_shape(a, 3, 3, name);
  mxi::Mat3 m;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      m(i, j) = a.at(i, j);
  return m;
}

Array to_array(const mxi::Mat3 &m) {
  Array a({py::ssize_t{3}, py::ssize_t{3}});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < 3; ++i)
    for (py::ssize_t j = 0; j < 3; ++j)
      r(i, j) = m(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
  return a;
}

using FlagArray =
    py::array_t<std::uint64_t, py::array::c_style | py::array::forcecast>;

// Indexed spots from the reflection table's miller_index, panel,
// xyzobs.mm.value, xyzobs.mm.variance and flags columns.
IndexedSpots make_indexed(const IndexArray &miller, const IndexArray &panel,
                          const Array &xyz_mm, const Array &var_mm,
                          const FlagArray &flags) {
  const py::ssize_t n = panel.ndim() == 1 ? panel.shape(0) : -1;
  require_shape(miller, n, 3, "miller_index");
  require_shape(panel, n, 0, "panel");
  require_shape(xyz_mm, n, 3, "xyzobs.mm.value");
  require_shape(var_mm, n, 3, "xyzobs.mm.variance");
  require_shape(flags, n, 0, "flags");
  IndexedSpots s;
  const auto k = static_cast<std::size_t>(n);
  s.miller.resize(k);
  s.panel.assign(panel.data(), panel.data() + n);
  s.flags.assign(flags.data(), flags.data() + n);
  s.x_mm.resize(k);
  s.y_mm.resize(k);
  s.var_x_mm.resize(k);
  s.var_y_mm.resize(k);
  for (py::ssize_t i = 0; i < n; ++i) {
    const auto j = static_cast<std::size_t>(i);
    for (py::ssize_t c = 0; c < 3; ++c)
      s.miller[j][static_cast<std::size_t>(c)] =
          static_cast<int>(miller.at(i, c));
    s.x_mm[j] = xyz_mm.at(i, 0);
    s.y_mm[j] = xyz_mm.at(i, 1);
    s.var_x_mm[j] = var_mm.at(i, 0);
    s.var_y_mm[j] = var_mm.at(i, 1);
  }
  return s;
}

AxisDerivative axis_derivative(const std::string &name) {
  if (name == "dials")
    return AxisDerivative::dials;
  if (name == "dials_limit")
    return AxisDerivative::dials_limit;
  if (name == "analytic")
    return AxisDerivative::analytic;
  throw py::value_error(
      "axis_derivative: expected 'dials', 'dials_limit' or 'analytic'");
}

Array to_array(const std::vector<double> &v, py::ssize_t cols) {
  const py::ssize_t rows = cols ? static_cast<py::ssize_t>(v.size()) / cols : 0;
  Array a({rows, cols});
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

py::dict to_dict(const RefineResult &r) {
  py::dict out;
  out["status"] = std::string(to_string(r.status));
  out["failed"] = failed(r.status);
  out["message"] = r.message;
  out["U"] = to_array(r.U);
  out["B"] = to_array(r.B);
  out["parameter_names"] = r.parameter_names;
  out["n_iterations"] = r.n_iterations;
  out["n_used"] = r.n_used;
  out["rmsd"] = py::make_tuple(r.rmsd[0], r.rmsd[1], r.rmsd[2]);
  py::list journal;
  for (const JournalRow &row : r.journal) {
    py::dict d;
    d["rmsd"] = py::make_tuple(row.rmsd[0], row.rmsd[1], row.rmsd[2]);
    d["objective"] = row.objective;
    d["mu"] = row.mu;
    d["nu"] = row.nu;
    d["parameters"] = row.parameters;
    journal.append(d);
  }
  out["journal"] = journal;
  IndexArray used(static_cast<py::ssize_t>(r.used.size()));
  for (std::size_t i = 0; i < r.used.size(); ++i)
    used.mutable_at(static_cast<py::ssize_t>(i)) =
        static_cast<std::int64_t>(r.used[i]);
  out["used"] = used;
  out["xyzcal.mm"] = to_array(r.x_calc, r.y_calc);
  Array delpsi(static_cast<py::ssize_t>(r.delpsi.size()));
  for (std::size_t i = 0; i < r.delpsi.size(); ++i)
    delpsi.mutable_at(static_cast<py::ssize_t>(i)) = r.delpsi[i];
  out["delpsical.rad"] = delpsi;
  out["s1"] = to_array(r.s1);
  return out;
}

std::vector<double> column(const Array &a, py::ssize_t rows, py::ssize_t col,
                           const char *name) {
  require_shape(a, rows, a.ndim() == 2 ? a.shape(1) : 0, name);
  std::vector<double> v(static_cast<std::size_t>(rows));
  for (py::ssize_t i = 0; i < rows; ++i)
    v[static_cast<std::size_t>(i)] = a.ndim() == 2 ? a.at(i, col) : a.at(i);
  return v;
}

py::array_t<bool> to_bool_array(const std::vector<bool> &v) {
  py::array_t<bool> a(static_cast<py::ssize_t>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i)
    a.mutable_at(static_cast<py::ssize_t>(i)) = v[i];
  return a;
}

py::dict to_dict(const lbfgs::Result &r) {
  py::dict out;
  out["failure"] = std::string(lbfgs::to_string(r.failure));
  out["raised"] = r.raised;
  out["iter"] = r.iter;
  out["nfun"] = r.nfun;
  return out;
}

std::vector<std::array<int, 3>> miller_of(const IndexArray &miller) {
  const py::ssize_t n = miller.ndim() == 2 ? miller.shape(0) : -1;
  require_shape(miller, n, 3, "miller_index");
  std::vector<std::array<int, 3>> h(static_cast<std::size_t>(n));
  for (py::ssize_t i = 0; i < n; ++i)
    for (py::ssize_t c = 0; c < 3; ++c)
      h[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
          static_cast<int>(miller.at(i, c));
  return h;
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

  m.def(
      "refine_still",
      [](const StillsDetector &detector, const Array &s0, const Array &U,
         const Array &B, const std::string &hall, const IndexArray &miller,
         const IndexArray &panel, const Array &xyz_mm, const Array &var_mm,
         const FlagArray &flags, const std::string &axis, double rmsd_tolerance,
         int max_iterations) {
        RefineOptions options;
        options.axis_derivative = axis_derivative(axis);
        options.rmsd_tolerance = rmsd_tolerance;
        options.max_iterations = max_iterations;
        const IndexedSpots spots =
            make_indexed(miller, panel, xyz_mm, var_mm, flags);
        RefineResult r;
        {
          py::gil_scoped_release release;
          r = refine_crystal(detector, vec3(s0, "s0"), mat3(U, "U"),
                             mat3(B, "B"), hall, spots, options);
        }
        return to_dict(r);
      },
      py::arg("detector"), py::arg("s0"), py::arg("U"), py::arg("B"),
      py::arg("hall"), py::arg("miller_index"), py::arg("panel"),
      py::arg("xyz_mm"), py::arg("var_mm"), py::arg("flags"),
      py::arg("axis_derivative") = "dials", py::arg("rmsd_tolerance") = 1e-4,
      py::arg("max_iterations") = 100,
      "Refine a still's crystal (U, B; the cell constrained by the space "
      "group named by its Hall symbol) against its indexed spots, as dials' "
      "stills refinement does, with the beam and detector fixed.");

  m.def(
      "linearise_still",
      [](const StillsDetector &detector, const Array &s0, const Array &U,
         const Array &B, const std::string &hall, const IndexArray &miller,
         const IndexArray &panel, const Array &xyz_mm, const Array &var_mm,
         const FlagArray &flags, const std::string &axis) {
        RefineOptions options;
        options.axis_derivative = axis_derivative(axis);
        const Linearisation l = linearise_crystal(
            detector, vec3(s0, "s0"), mat3(U, "U"), mat3(B, "B"), hall,
            make_indexed(miller, panel, xyz_mm, var_mm, flags), options);
        py::dict out;
        out["status"] = std::string(to_string(l.status));
        out["message"] = l.message;
        out["parameter_names"] = l.parameter_names;
        out["parameters"] = l.parameters;
        out["used"] = l.used;
        const auto n_param = static_cast<py::ssize_t>(l.parameters.size());
        out["residual"] = to_array(l.residual, 3);
        out["weight"] = to_array(l.weight, 3);
        out["jacobian"] = to_array(l.jacobian, n_param);
        return out;
      },
      py::arg("detector"), py::arg("s0"), py::arg("U"), py::arg("B"),
      py::arg("hall"), py::arg("miller_index"), py::arg("panel"),
      py::arg("xyz_mm"), py::arg("var_mm"), py::arg("flags"),
      py::arg("axis_derivative") = "dials",
      "The residuals, weights and Jacobian (rows X, Y, Delta-Psi per used "
      "spot; one column per parameter) at the starting model.");

  m.def(
      "predict_still",
      [](const StillsDetector &detector, const Array &s0, const Array &U,
         const Array &B, const IndexArray &miller, const IndexArray &panel) {
        const py::ssize_t n = panel.ndim() == 1 ? panel.shape(0) : -1;
        require_shape(miller, n, 3, "miller_index");
        require_shape(panel, n, 0, "panel");
        std::vector<std::array<int, 3>> h(static_cast<std::size_t>(n));
        for (py::ssize_t i = 0; i < n; ++i)
          for (py::ssize_t c = 0; c < 3; ++c)
            h[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
                static_cast<int>(miller.at(i, c));
        const std::vector<std::int64_t> p(panel.data(), panel.data() + n);
        const Predictions r = predict_spots(detector, vec3(s0, "s0"),
                                            mat3(U, "U"), mat3(B, "B"), h, p);
        py::dict out;
        IndexArray ok(n);
        Array delpsi(n);
        for (py::ssize_t i = 0; i < n; ++i) {
          const auto j = static_cast<std::size_t>(i);
          ok.mutable_at(i) = r.status[j] == PredictStatus::ok ? 1 : 0;
          delpsi.mutable_at(i) = r.delpsi[j];
        }
        out["predicted"] = ok;
        out["xyzcal.mm"] = to_array(r.x_mm, r.y_mm);
        out["delpsical.rad"] = delpsi;
        out["s1"] = to_array(r.s1);
        return out;
      },
      py::arg("detector"), py::arg("s0"), py::arg("U"), py::arg("B"),
      py::arg("miller_index"), py::arg("panel"),
      "Stills predictions (X, Y in mm, Delta-Psi, s1) of each reflection for "
      "the crystal A = U B; predicted is 0 where a reflection cannot be "
      "predicted or misses its panel's plane.");

  m.def(
      "mm_to_px",
      [](const StillsDetector &detector, const IndexArray &panel,
         const Array &xyz_mm) {
        const py::ssize_t n = panel.ndim() == 1 ? panel.shape(0) : -1;
        require_shape(panel, n, 0, "panel");
        require_shape(xyz_mm, n, xyz_mm.ndim() == 2 ? xyz_mm.shape(1) : 0,
                      "xyz_mm");
        Array out({n, py::ssize_t{2}});
        auto r = out.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < n; ++i) {
          const auto [x, y] = detector.panel(panel.at(i))
                                  .mm_to_px(xyz_mm.at(i, 0), xyz_mm.at(i, 1));
          r(i, 0) = x;
          r(i, 1) = y;
        }
        return out;
      },
      py::arg("detector"), py::arg("panel"), py::arg("xyz_mm"),
      "Millimetre positions on their panels to pixels, as dxtbx's "
      "millimeter_to_pixel does with the parallax correction.");

  m.def(
      "sauter_poon",
      [](const Array &obs_mm, const Array &calc_mm) {
        const py::ssize_t n = obs_mm.ndim() == 2 ? obs_mm.shape(0) : -1;
        const auto x_obs = column(obs_mm, n, 0, "obs_mm");
        const auto y_obs = column(obs_mm, n, 1, "obs_mm");
        const auto x_calc = column(calc_mm, n, 0, "calc_mm");
        const auto y_calc = column(calc_mm, n, 1, "calc_mm");
        SauterPoonResult r;
        {
          py::gil_scoped_release release;
          r = sauter_poon(x_obs, y_obs, x_calc, y_calc);
        }
        py::dict out;
        out["status"] = std::string(to_string(r.status));
        out["kept"] = to_bool_array(r.kept);
        out["sigma"] = r.sigma;
        out["sigma_start"] = r.sigma_start;
        out["sd"] = r.sd;
        out["n_fit"] = r.n_fit;
        out["first_rejected"] = r.first_rejected;
        out["lbfgs"] = to_dict(r.minimiser);
        return out;
      },
      py::arg("obs_mm"), py::arg("calc_mm"),
      "Sauter-Poon outlier rejection on observed and calculated spot "
      "positions (n, 2) in mm, as dials' stills indexer runs it.");

  m.def(
      "nave_parameters",
      [](const Array &delpsi, const Array &d, const Array &two_theta,
         double volume, double wavelength, double ewald_proximal_volume_max,
         double resolution_cutoff, double expansion, bool safe_likelihood) {
        const py::ssize_t n = delpsi.ndim() == 1 ? delpsi.shape(0) : -1;
        const auto psi = column(delpsi, n, 0, "delpsi");
        const auto dv = column(d, n, 0, "d");
        const auto tt = column(two_theta, n, 0, "two_theta");
        NaveOptions options;
        options.ewald_proximal_volume_max = ewald_proximal_volume_max;
        options.ewald_proximity_resolution_cutoff = resolution_cutoff;
        options.model_expansion_factor = expansion;
        options.safe_likelihood = safe_likelihood;
        NaveResult r;
        {
          py::gil_scoped_release release;
          r = nave_parameters(psi, dv, tt, volume, wavelength, options);
        }
        py::dict out;
        out["status"] = std::string(to_string(r.status));
        out["failed"] = failed(r.status);
        out["eta_seed"] = r.eta_seed;
        out["domain_size_seed"] = r.domain_size_seed;
        out["alpha"] = r.alpha;
        out["eta"] = r.eta;
        out["accepted"] = to_bool_array(r.accepted);
        out["half_mosaicity_deg"] = r.half_mosaicity_deg;
        out["domain_size_ang"] = r.domain_size_ang;
        out["ewald_volume"] = r.ewald_volume;
        out["lbfgs"] = to_dict(r.minimiser);
        return out;
      },
      py::arg("delpsi"), py::arg("d"), py::arg("two_theta"), py::arg("volume"),
      py::arg("wavelength"), py::arg("ewald_proximal_volume_max") = 0.0025,
      py::arg("resolution_cutoff") = 2.0, py::arg("expansion") = 1.4,
      py::arg("safe_likelihood") = false,
      "dials' Nave mosaicity and domain-size fit of a crystal's Delta-Psi "
      "(radians) against d (A) and 2-theta (degrees), with its acceptance "
      "flags and Ewald proximal volume test.");

  m.def(
      "cell_d_spacings",
      [](const Array &B, const IndexArray &miller) {
        const auto d = cell_d_spacings(mat3(B, "B"), miller_of(miller));
        return Array(static_cast<py::ssize_t>(d.size()), d.data());
      },
      py::arg("B"), py::arg("miller_index"), "d = 1 / |B h| (A).");

  m.def(
      "two_theta_deg",
      [](const Array &d, double wavelength) {
        const py::ssize_t n = d.ndim() == 1 ? d.shape(0) : -1;
        const auto t = two_theta_deg(column(d, n, 0, "d"), wavelength);
        return Array(static_cast<py::ssize_t>(t.size()), t.data());
      },
      py::arg("d"), py::arg("wavelength"), "2 asin(lambda / 2d) in degrees.");

  m.def(
      "cell_volume", [](const Array &B) { return cell_volume(mat3(B, "B")); },
      py::arg("B"), "The direct cell volume 1 / det B (A^3).");
}

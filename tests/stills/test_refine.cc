// Stills prediction, the crystal parameterisation and the refinement, against
// finite differences and against dials' own runs on real stills.

#include <cmath>
#include <string>
#include <vector>

#include "check.hh"
#include "real_refine.hh"
#include "stills/crystal.hh"
#include "stills/predict.hh"
#include "stills/refine.hh"
#include "symmetry.hh"

using namespace mxi::stills;
using mxi::Mat3;
using mxi::Vec3;

namespace {

StillsDetector detector_of(const real_refine::Run &r) {
  StillsDetector d;
  for (const auto &p : r.panels) {
    StillsPanel q;
    q.fast = p.fast;
    q.slow = p.slow;
    q.origin = p.origin;
    q.pixel_size[0] = p.pixel_size[0];
    q.pixel_size[1] = p.pixel_size[1];
    q.mu = p.mu;
    q.t0 = p.t0;
    d.panels.push_back(q);
  }
  return d;
}

IndexedSpots spots_of(const real_refine::Run &r) {
  IndexedSpots s;
  for (const auto &x : r.spots) {
    s.miller.push_back(x.h);
    s.panel.push_back(x.panel);
    s.x_mm.push_back(x.x_mm);
    s.y_mm.push_back(x.y_mm);
    s.var_x_mm.push_back(x.var_x_mm);
    s.var_y_mm.push_back(x.var_y_mm);
    s.flags.push_back(x.flags);
  }
  return s;
}

Mat3 add(const Mat3 &a, const Mat3 &b) {
  Mat3 r;
  for (std::size_t i = 0; i < 9; ++i)
    r.m[i] = a.m[i] + b.m[i];
  return r;
}

Mat3 d_inverse(const StillsPanel &p) {
  return Mat3::from_columns(p.fast, p.slow, p.origin).inverse();
}

double misorientation_deg(const Mat3 &U1, const Mat3 &U2) {
  const Mat3 R = U1 * U2.transpose();
  const double s = 0.5 * std::sqrt(std::pow(R(2, 1) - R(1, 2), 2) +
                                   std::pow(R(0, 2) - R(2, 0), 2) +
                                   std::pow(R(1, 0) - R(0, 1), 2));
  const double c = 0.5 * (R(0, 0) + R(1, 1) + R(2, 2) - 1.0);
  return std::atan2(s, c) * 180.0 / M_PI;
}

double objective(const Linearisation &l) {
  double f = 0.0;
  for (std::size_t e = 0; e < l.residual.size(); ++e)
    f += l.weight[e] * l.residual[e] * l.residual[e];
  return f / 2;
}

Mat3 symmetric(const Metric &g) {
  return {g[0], g[3], g[4], g[3], g[1], g[5], g[4], g[5], g[2]};
}

// dA for each parameter at the starting model, per mrad and per unit of the
// scaled metric element, as the refinement parameterises them.
std::vector<Mat3> parameter_directions(const real_refine::Run &r) {
  const CellConstraint c = CellConstraint::from_hall(r.hall);
  const Orientation o = compose_orientation(r.U0, {0.0, 0.0, 0.0});
  const Mat3 B = b_from_metric(metric_of(r.B0));
  std::vector<Mat3> dA;
  for (int k = 0; k < 3; ++k)
    dA.push_back(o.dU[static_cast<std::size_t>(k)] * B);
  for (const Metric &basis : c.basis)
    dA.push_back(r.U0 * db_from_metric(B, B.inverse(), basis) * 1.0e-5);
  return dA;
}

// The largest difference, relative to the largest entry, between ray
// derivatives in the given mode and central differences of the prediction,
// over every spot, parameter and quantity of a run.
double worst_against_finite_differences(const real_refine::Run &r,
                                        AxisDerivative mode) {
  const StillsDetector det = detector_of(r);
  const Beam beam(r.s0);
  const Mat3 A = r.U0 * b_from_metric(metric_of(r.B0));
  double worst = 0.0;
  for (const Mat3 &dA : parameter_directions(r)) {
    const double h = 1e-4;
    std::vector<std::array<double, 3>> fd, an;
    double scale = 0.0;
    for (const auto &s : r.spots) {
      const Mat3 D = d_inverse(det.panel(s.panel));
      Ray ray, plus, minus;
      check::is_true(predict(beam, A, s.h, D, &ray) == PredictStatus::ok &&
                         predict(beam, add(A, dA * h), s.h, D, &plus) ==
                             PredictStatus::ok &&
                         predict(beam, A - dA * h, s.h, D, &minus) ==
                             PredictStatus::ok,
                     r.name + ": predicted");
      const Vec3 hv{static_cast<double>(s.h[0]), static_cast<double>(s.h[1]),
                    static_cast<double>(s.h[2])};
      const RayDerivative d = ray_derivative(beam, ray, D, dA * hv, mode);
      an.push_back({d.dx, d.dy, d.ddelpsi});
      fd.push_back({(plus.x_mm - minus.x_mm) / (2 * h),
                    (plus.y_mm - minus.y_mm) / (2 * h),
                    (plus.delpsi - minus.delpsi) / (2 * h)});
      for (double v : fd.back())
        scale = std::fmax(scale, std::abs(v));
    }
    for (std::size_t i = 0; i < fd.size(); ++i)
      for (int k = 0; k < 3; ++k)
        worst = std::fmax(worst, std::abs(an[i][k] - fd[i][k]) / scale);
  }
  return worst;
}

} // namespace

// The free elements of G* are the non-pivot columns of the reduced row-echelon
// form, in column order (a*a*, b*b*, c*c*, a*b*, a*c*, b*c*): what dials
// refines, measured on every captured run (STATUS, R6).
TEST(stills_cell_constraint_free_elements) {
  const std::vector<std::pair<std::string, std::vector<int>>> cases = {
      {" P 1", {0, 1, 2, 3, 4, 5}},
      {" P 2yb", {0, 1, 2, 4}},
      {" P 2ac 2ab", {0, 1, 2}},
      {" P 4nw 2abw", {1, 2}},
      {"-P 4 2", {1, 2}},
      {"-P 6 2", {2, 3}},
      {" P 2 2 3", {2}},
  };
  for (const auto &[hall, columns] : cases) {
    const CellConstraint c = CellConstraint::from_hall(hall);
    check::equal(static_cast<long long>(c.column.size()),
                 static_cast<long long>(columns.size()), hall + ": count");
    for (std::size_t k = 0; k < columns.size(); ++k)
      check::equal(c.column[k], columns[k], hall + ": column");
    for (const auto &rotation : mxi::SpaceGroup::from_hall(hall).rotations()) {
      Mat3 R;
      for (int i = 0; i < 9; ++i)
        R.m[static_cast<std::size_t>(i)] =
            rotation[static_cast<std::size_t>(i)];
      for (const Metric &basis : c.basis) {
        const Mat3 G = symmetric(basis);
        const Mat3 diff = R * G * R.transpose() - G;
        for (double v : diff.m)
          check::close(v, 0.0, 1e-15, hall + ": basis invariant");
      }
    }
  }
  // Hexagonal: a*.b* is free and sets a*a* = b*b* = 2 a*.b*.
  const CellConstraint hex = CellConstraint::from_hall("-P 6 2");
  const Metric expect{2, 2, 0, 1, 0, 0};
  for (int e = 0; e < 6; ++e)
    check::close(hex.basis[1][static_cast<std::size_t>(e)],
                 expect[static_cast<std::size_t>(e)], 0.0, "hexagonal basis");
}

TEST(stills_b_from_metric_is_dxtbx_b) {
  for (const auto &r : real_refine::all()) {
    const Mat3 B = b_from_metric(metric_of(r.B0));
    for (std::size_t i = 0; i < 9; ++i)
      check::close(B.m[i], r.B0.m[i], 1e-17, r.name + ": B");
  }
}

TEST(stills_b_from_metric_refuses_impossible_cells) {
  const auto refuses = [](const Metric &g) {
    try {
      b_from_metric(g);
    } catch (const ComposeError &) {
      return true;
    }
    return false;
  };
  check::is_true(refuses({1e-4, 1e-4, -1e-4, 0, 0, 0}), "negative c*c*");
  check::is_true(refuses({1e-4, 1e-4, 1e-4, 0, 0, 2e-4}), "not definite");
  check::is_true(refuses({1, 1, 1, 0, 0, 0}), "a 1 A^3 cell");
  check::is_true(!refuses({1e-4, 1e-4, 1e-4, 0, 0, 0}), "a 1e6 A^3 cell");
}

TEST(stills_db_from_metric_matches_finite_differences) {
  for (const auto &r : real_refine::all()) {
    const Metric g = metric_of(r.B0);
    const Mat3 B = b_from_metric(g);
    for (int e = 0; e < 6; ++e) {
      Metric dg{};
      dg[static_cast<std::size_t>(e)] = 1.0;
      const double h = 1e-9;
      Metric gp = g, gm = g;
      gp[static_cast<std::size_t>(e)] += h;
      gm[static_cast<std::size_t>(e)] -= h;
      const Mat3 fd = (b_from_metric(gp) - b_from_metric(gm)) * (1.0 / (2 * h));
      const Mat3 an = db_from_metric(B, B.inverse(), dg);
      double scale = 0.0;
      for (double v : fd.m)
        scale = std::fmax(scale, std::abs(v));
      for (std::size_t i = 0; i < 9; ++i)
        check::close(an.m[i], fd.m[i], 1e-6 * scale, r.name + ": dB/dg");
    }
  }
}

TEST(stills_orientation_derivatives_match_finite_differences) {
  const Mat3 U0 = real_refine::cytochrome().U0;
  const std::array<double, 3> at{0.3, -0.2, 0.1};
  const Orientation o = compose_orientation(U0, at);
  for (std::size_t k = 0; k < 3; ++k) {
    auto p = at, m = at;
    const double h = 1e-4;
    p[k] += h;
    m[k] -= h;
    const Mat3 fd =
        (compose_orientation(U0, p).U - compose_orientation(U0, m).U) *
        (1.0 / (2 * h));
    for (std::size_t i = 0; i < 9; ++i)
      check::close(o.dU[k].m[i], fd.m[i], 1e-10, "dU/dphi");
  }
}

// The analytic derivative of X, Y and Delta-Psi is the derivative of the
// prediction. dials' form of the rotation-axis term is not: it leaves out a
// factor 1 / |q0 x s0| in the axis's change.
TEST(stills_ray_derivatives_against_finite_differences) {
  double worst_dials = 0.0;
  for (const auto &r : real_refine::all()) {
    check::close(worst_against_finite_differences(r, AxisDerivative::analytic),
                 0.0, 1e-7, r.name + ": analytic derivative");
    worst_dials = std::fmax(worst_dials, worst_against_finite_differences(
                                             r, AxisDerivative::dials_limit));
  }
  check::is_true(worst_dials > 1e-5,
                 "dials' axis term differs from the derivative by " +
                     std::to_string(worst_dials));
}

// The largest difference over every run of the R6 replay is 1.25e-12 mm and
// 6e-16 rad.
TEST(stills_predictions_match_dials) {
  for (const auto &r : real_refine::all()) {
    const IndexedSpots s = spots_of(r);
    const Predictions p =
        predict_spots(detector_of(r), r.s0, r.U1, r.B1, s.miller, s.panel);
    for (std::size_t i = 0; i < r.spots.size(); ++i) {
      check::is_true(p.status[i] == PredictStatus::ok, r.name + ": predicted");
      check::close(p.x_mm[i], r.spots[i].x_cal, 2e-12, r.name + ": X");
      check::close(p.y_mm[i], r.spots[i].y_cal, 2e-12, r.name + ": Y");
      check::close(p.delpsi[i], r.spots[i].delpsi, 1e-15,
                   r.name + ": Delta-Psi");
    }
  }
}

TEST(stills_starting_objective_matches_dials) {
  for (const auto &r : real_refine::all()) {
    const Linearisation l = linearise_crystal(detector_of(r), r.s0, r.U0, r.B0,
                                              r.hall, spots_of(r));
    check::close(objective(l), r.objective0, 1e-11 * r.objective0,
                 r.name + ": objective");
    check::equal(static_cast<long long>(l.parameters.size()), r.n_parameters,
                 r.name + ": parameters");
  }
}

// dials stops on a relative change in RMSD below 1e-4, short of the minimum,
// so the two end points are compared by the objective they reach.
TEST(stills_refinement_reaches_dials_objective) {
  for (const auto &r : real_refine::all()) {
    for (AxisDerivative mode :
         {AxisDerivative::dials, AxisDerivative::dials_limit,
          AxisDerivative::analytic}) {
      RefineOptions options;
      options.axis_derivative = mode;
      const RefineResult out = refine_crystal(detector_of(r), r.s0, r.U0, r.B0,
                                              r.hall, spots_of(r), options);
      check::is_true(out.status == RefineStatus::rmsd_converged,
                     r.name + ": " + to_string(out.status));
      const double f_dials = objective(linearise_crystal(
          detector_of(r), r.s0, r.U1, r.B1, r.hall, spots_of(r)));
      const double f_mxi = objective(linearise_crystal(
          detector_of(r), r.s0, out.U, out.B, r.hall, spots_of(r)));
      check::close(f_mxi, f_dials, 2e-4 * f_dials, r.name + ": objective");
      check::is_true(misorientation_deg(out.U, r.U1) < 1e-3,
                     r.name + ": orientation");
    }
  }
}

// Observations made from a known crystal, without Delta-Psi's restraint, have
// that crystal as their exact minimum.
TEST(stills_refinement_recovers_a_known_crystal) {
  for (const auto &r : real_refine::all()) {
    IndexedSpots s = spots_of(r);
    const StillsDetector det = detector_of(r);
    const Predictions truth =
        predict_spots(det, r.s0, r.U1, r.B1, s.miller, s.panel);
    s.x_mm = truth.x_mm;
    s.y_mm = truth.y_mm;

    const Mat3 U = mxi::rotation({0.3, 0.5, 0.8}, 0.02 * M_PI / 180) * r.U1;
    const CellConstraint c = CellConstraint::from_hall(r.hall);
    std::vector<double> g = c.independent(metric_of(r.B1));
    for (double &v : g)
      v *= 1.002;
    const Mat3 B = b_from_metric(c.metric(g));

    RefineOptions options;
    options.axis_derivative = AxisDerivative::analytic;
    options.delpsi_weight = 0.0;
    options.rmsd_tolerance = 0.0;
    const RefineResult out =
        refine_crystal(det, r.s0, U, B, r.hall, s, options);
    check::is_true(!failed(out.status), r.name + ": " + out.message);
    check::close(misorientation_deg(out.U, r.U1), 0.0, 1e-7,
                 r.name + ": orientation recovered");
    const Metric got = metric_of(out.B), want = metric_of(r.B1);
    for (std::size_t e = 0; e < 6; ++e)
      check::close(got[e], want[e], 1e-9 * want[0], r.name + ": cell");
  }
}

TEST(stills_refinement_failures) {
  const real_refine::Run r = real_refine::lysozyme();
  const StillsDetector det = detector_of(r);

  IndexedSpots none = spots_of(r);
  for (auto &h : none.miller)
    h = {0, 0, 0};
  check::is_true(refine_crystal(det, r.s0, r.U0, r.B0, r.hall, none).status ==
                     RefineStatus::no_reflections,
                 "unindexed spots are not used");

  IndexedSpots zero = spots_of(r);
  zero.var_x_mm[3] = 0.0;
  zero.var_y_mm[3] = 0.0;
  check::is_true(refine_crystal(det, r.s0, r.U0, r.B0, r.hall, zero).status ==
                     RefineStatus::zero_weights,
                 "zero variance in X and Y");

  IndexedSpots three = spots_of(r);
  for (auto *v : {&three.x_mm, &three.y_mm, &three.var_x_mm, &three.var_y_mm})
    v->resize(3);
  three.miller.resize(3);
  three.panel.resize(3);
  three.flags.resize(3);
  const RefineResult few = refine_crystal(det, r.s0, r.U0, r.B0, r.hall, three);
  check::is_true(few.status == RefineStatus::too_few_degrees_of_freedom,
                 std::string("nine equations, nine parameters: ") +
                     to_string(few.status));

  IndexedSpots far = spots_of(r);
  far.miller[0] = {400, 400, 400};
  check::is_true(refine_crystal(det, r.s0, r.U0, r.B0, r.hall, far).status ==
                     RefineStatus::prediction_failed,
                 "a reflection beyond the Ewald sphere");

  IndexedSpots overloaded = spots_of(r);
  for (auto &f : overloaded.flags)
    f |= std::uint64_t{1} << 10;
  check::is_true(
      refine_crystal(det, r.s0, r.U0, r.B0, r.hall, overloaded).status ==
          RefineStatus::no_reflections,
      "overloaded spots are not used");
}

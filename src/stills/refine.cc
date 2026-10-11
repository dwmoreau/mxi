#include "refine.hh"

#include <algorithm>
#include <cmath>

namespace mxi::stills {

namespace {

constexpr std::uint64_t kOverloaded = std::uint64_t{1} << 10;
constexpr std::uint64_t kNotSuitableForRefinement = std::uint64_t{1} << 27;
constexpr double kCellScale = 1.0e5;
// A parameter that moves no residual by more than this is fixed.
constexpr double kAutoReduceThreshold = 1e-10;

// The inverse as a cofactor matrix divided by the determinant.
Mat3 inverse_of(const Mat3 &d) {
  const auto &m = d.m;
  const double det = d.determinant();
  return {(m[4] * m[8] - m[5] * m[7]) / det, (m[2] * m[7] - m[1] * m[8]) / det,
          (m[1] * m[5] - m[2] * m[4]) / det, (m[5] * m[6] - m[3] * m[8]) / det,
          (m[0] * m[8] - m[2] * m[6]) / det, (m[2] * m[3] - m[0] * m[5]) / det,
          (m[3] * m[7] - m[4] * m[6]) / det, (m[1] * m[6] - m[0] * m[7]) / det,
          (m[0] * m[4] - m[1] * m[3]) / det};
}

Mat3 d_matrix_inverse(const StillsPanel &p) {
  return inverse_of(Mat3::from_columns(p.fast, p.slow, p.origin));
}

// Thrown inside a refinement to end it with a failure status.
struct Failure {
  RefineStatus status;
  std::string message;
};

// The model, its prediction of every used spot and the linearised equations,
// for one parameter vector.
class Problem {
public:
  Problem(const StillsDetector &detector, const Vec3 &s0, const Mat3 &U0,
          const CellConstraint &constraint, const IndexedSpots &spots,
          std::vector<std::size_t> used, const RefineOptions &options)
      : beam_(s0), U0_(U0), constraint_(constraint), spots_(spots),
        used_(std::move(used)), options_(options) {
    const std::size_t n = used_.size();
    D_.resize(n);
    std::vector<std::int64_t> seen_panel;
    std::vector<Mat3> seen_D;
    for (std::size_t k = 0; k < n; ++k) {
      const std::int64_t p = spots_.panel[used_[k]];
      auto it = std::find(seen_panel.begin(), seen_panel.end(), p);
      if (it == seen_panel.end()) {
        seen_panel.push_back(p);
        seen_D.push_back(d_matrix_inverse(detector.panel(p)));
        it = seen_panel.end() - 1;
      }
      D_[k] = seen_D[static_cast<std::size_t>(it - seen_panel.begin())];
    }
    rays_.resize(n);
    w_.resize(3 * n);
    for (std::size_t k = 0; k < n; ++k) {
      const std::size_t i = used_[k];
      const double vx = spots_.var_x_mm[i], vy = spots_.var_y_mm[i];
      w_[3 * k] = vx > 0.0 ? 1.0 / vx : 0.0;
      w_[3 * k + 1] = vy > 0.0 ? 1.0 / vy : 0.0;
      w_[3 * k + 2] = options_.delpsi_weight;
      if (w_[3 * k] == 0.0 && w_[3 * k + 1] == 0.0)
        throw Failure{RefineStatus::zero_weights,
                      "an indexed reflection has zero variance in X and Y"};
    }
    residual_.resize(3 * n);
  }

  std::size_t n_equations() const { return residual_.size(); }

  // Composes the model at x and predicts every used spot. With derivatives,
  // also fills the Jacobian columns of the free parameters.
  void evaluate(const std::vector<double> &x, const std::vector<int> &free,
                bool derivatives) {
    const Orientation o = compose_orientation(U0_, {x[0], x[1], x[2]});
    std::vector<double> g(constraint_.size());
    for (std::size_t k = 0; k < g.size(); ++k)
      g[k] = x[3 + k] * 1.0e-5;
    Mat3 B;
    try {
      B = b_from_metric(constraint_.metric(g));
    } catch (const ComposeError &e) {
      throw Failure{RefineStatus::cannot_compose, e.what()};
    }
    U_ = o.U;
    B_ = B;
    const Mat3 A = o.U * B;

    const std::size_t n = used_.size();
    for (std::size_t k = 0; k < n; ++k) {
      const std::size_t i = used_[k];
      const PredictStatus s =
          predict(beam_, A, spots_.miller[i], D_[k], &rays_[k]);
      if (s != PredictStatus::ok)
        throw Failure{RefineStatus::prediction_failed,
                      s == PredictStatus::misses_panel
                          ? "a predicted ray misses its panel's plane"
                          : "a reflection cannot be predicted"};
      residual_[3 * k] = rays_[k].x_mm - spots_.x_mm[i];
      residual_[3 * k + 1] = rays_[k].y_mm - spots_.y_mm[i];
      residual_[3 * k + 2] = rays_[k].delpsi;
    }
    if (!derivatives)
      return;

    // dq/dp for each free parameter: (dU/dphi) B h, or U (dB/dg) h.
    const Mat3 B_inverse = inverse_of(B);
    std::vector<Mat3> dA(free.size());
    for (std::size_t c = 0; c < free.size(); ++c) {
      const int p = free[c];
      if (p < 3) {
        dA[c] = o.dU[static_cast<std::size_t>(p)] * B;
      } else {
        const Metric dg = constraint_.basis[static_cast<std::size_t>(p - 3)];
        Mat3 dB = db_from_metric(B, B_inverse, dg);
        for (double &v : dB.m)
          v *= 1.0e-5;
        dA[c] = o.U * dB;
      }
    }
    jacobian_.assign(3 * n * free.size(), 0.0);
    for (std::size_t k = 0; k < n; ++k) {
      const auto &h = spots_.miller[used_[k]];
      const Vec3 hv{static_cast<double>(h[0]), static_cast<double>(h[1]),
                    static_cast<double>(h[2])};
      for (std::size_t c = 0; c < free.size(); ++c) {
        const RayDerivative d = ray_derivative(
            beam_, rays_[k], D_[k], dA[c] * hv, options_.axis_derivative);
        jacobian_[(3 * k) * free.size() + c] = d.dx;
        jacobian_[(3 * k + 1) * free.size() + c] = d.dy;
        jacobian_[(3 * k + 2) * free.size() + c] = d.ddelpsi;
      }
    }
  }

  double objective() const {
    double r_sq = 0.0;
    for (std::size_t e = 0; e < residual_.size(); ++e)
      r_sq += w_[e] * residual_[e] * residual_[e];
    return r_sq / 2;
  }

  std::array<double, 3> rmsd() const {
    std::array<double, 3> s{};
    for (std::size_t k = 0; k < used_.size(); ++k)
      for (int j = 0; j < 3; ++j)
        s[j] += residual_[3 * k + j] * residual_[3 * k + j];
    const double n = static_cast<double>(used_.size());
    return {std::sqrt(s[0] / n), std::sqrt(s[1] / n), std::sqrt(s[2] / n)};
  }

  // The normal matrix J^T W J (full, symmetric) and the right-hand side
  // -J^T W r, over the free columns of the last evaluate.
  void normal_equations(std::size_t n_free, std::vector<double> *normal,
                        std::vector<double> *rhs) const {
    normal->assign(n_free * n_free, 0.0);
    rhs->assign(n_free, 0.0);
    for (std::size_t e = 0; e < residual_.size(); ++e) {
      const double *row = &jacobian_[e * n_free];
      const double w = w_[e];
      for (std::size_t a = 0; a < n_free; ++a) {
        (*rhs)[a] += -residual_[e] * row[a] * w;
        for (std::size_t b = a; b < n_free; ++b)
          (*normal)[a * n_free + b] += row[a] * w * row[b];
      }
    }
    for (std::size_t a = 0; a < n_free; ++a)
      for (std::size_t b = 0; b < a; ++b)
        (*normal)[a * n_free + b] = (*normal)[b * n_free + a];
  }

  // Whether any spot's X, Y or Delta-Psi moves by more than the threshold for
  // free column c of the last evaluate.
  bool column_moves_anything(std::size_t c, std::size_t n_free) const {
    for (std::size_t e = 0; e < residual_.size(); ++e)
      if (std::abs(jacobian_[e * n_free + c]) > kAutoReduceThreshold)
        return true;
    return false;
  }

  const std::vector<double> &residual() const { return residual_; }
  const std::vector<double> &weight() const { return w_; }
  const std::vector<double> &jacobian() const { return jacobian_; }
  const Mat3 &U() const { return U_; }
  const Mat3 &B() const { return B_; }
  const std::vector<Ray> &rays() const { return rays_; }
  const std::vector<std::size_t> &used() const { return used_; }

private:
  Beam beam_;
  Mat3 U0_;
  const CellConstraint &constraint_;
  const IndexedSpots &spots_;
  std::vector<std::size_t> used_;
  RefineOptions options_;
  std::vector<Mat3> D_;
  std::vector<Ray> rays_;
  std::vector<double> w_, residual_, jacobian_;
  Mat3 U_, B_;
};

// The RMSDs and predictions of the model the refinement ends on: the last one
// evaluated.
void record_predictions(const Problem &problem, RefineResult *result) {
  result->rmsd = problem.rmsd();
  result->used = problem.used();
  for (const Ray &ray : problem.rays()) {
    result->x_calc.push_back(ray.x_mm);
    result->y_calc.push_back(ray.y_mm);
    result->delpsi.push_back(ray.delpsi);
    result->s1.push_back(ray.s1);
  }
}

std::vector<double> free_values(const std::vector<double> &x,
                                const std::vector<int> &free) {
  std::vector<double> out(free.size());
  for (std::size_t c = 0; c < free.size(); ++c)
    out[c] = x[static_cast<std::size_t>(free[c])];
  return out;
}

bool rmsds_converged(const std::array<double, 3> &last,
                     const std::array<double, 3> &previous, double tolerance) {
  for (int j = 0; j < 3; ++j) {
    if (previous[j] > 0.0 &&
        !(std::abs((previous[j] - last[j]) / previous[j]) < tolerance))
      return false;
  }
  return true;
}

} // namespace

const char *to_string(RefineStatus status) {
  switch (status) {
  case RefineStatus::rmsd_converged:
    return "RMSD no longer decreasing";
  case RefineStatus::max_iterations:
    return "Reached maximum number of iterations";
  case RefineStatus::max_trial_iterations:
    return "Reached maximum number of consecutive unsuccessful trial steps";
  case RefineStatus::objective_flat:
    return "Objective function is flat";
  case RefineStatus::too_few_degrees_of_freedom:
    return "Not enough degrees of freedom to refine";
  case RefineStatus::no_reflections:
    return "no reflections to refine against";
  case RefineStatus::zero_weights:
    return "zero weights";
  case RefineStatus::indeterminate:
    return "the normal equations have an indeterminate solution";
  case RefineStatus::not_positive_definite:
    return "the normal matrix is not positive definite";
  case RefineStatus::cannot_compose:
    return "unable to compose the crystal model";
  case RefineStatus::prediction_failed:
    return "prediction failed";
  }
  return "unknown";
}

bool failed(RefineStatus status) {
  return status >= RefineStatus::no_reflections;
}

// The spots refinement uses: an index, and neither overloaded nor marked
// unsuitable.
std::vector<std::size_t> usable(const IndexedSpots &spots) {
  std::vector<std::size_t> used;
  for (std::size_t i = 0; i < spots.size(); ++i) {
    const auto &h = spots.miller[i];
    const std::uint64_t f = spots.flags.empty() ? 0 : spots.flags[i];
    if ((h[0] == 0 && h[1] == 0 && h[2] == 0) || (f & kOverloaded) ||
        (f & kNotSuitableForRefinement))
      continue;
    used.push_back(i);
  }
  return used;
}

std::vector<std::string> parameter_names(const CellConstraint &constraint) {
  std::vector<std::string> names{"Crystal1Phi1", "Crystal1Phi2",
                                 "Crystal1Phi3"};
  for (std::size_t k = 0; k < constraint.size(); ++k)
    names.push_back("Crystal1g_param_" + std::to_string(k));
  return names;
}

// The starting parameters: zero angles, and the free metric elements of B,
// scaled.
std::vector<double> start_parameters(const CellConstraint &constraint,
                                     const Mat3 &B) {
  std::vector<double> x(3 + constraint.size(), 0.0);
  const std::vector<double> g = constraint.independent(metric_of(B));
  for (std::size_t k = 0; k < g.size(); ++k)
    x[3 + k] = g[k] * kCellScale;
  return x;
}

Predictions predict_spots(const StillsDetector &detector, const Vec3 &s0,
                          const Mat3 &U, const Mat3 &B,
                          const std::vector<std::array<int, 3>> &miller,
                          const std::vector<std::int64_t> &panel) {
  const Beam beam(s0);
  const Mat3 A = U * B;
  const std::size_t n = miller.size();
  Predictions out;
  out.status.resize(n);
  out.x_mm.assign(n, 0.0);
  out.y_mm.assign(n, 0.0);
  out.delpsi.assign(n, 0.0);
  out.s1.resize(n);
  std::vector<std::int64_t> seen_panel;
  std::vector<Mat3> seen_D;
  for (std::size_t i = 0; i < n; ++i) {
    auto it = std::find(seen_panel.begin(), seen_panel.end(), panel[i]);
    if (it == seen_panel.end()) {
      seen_panel.push_back(panel[i]);
      seen_D.push_back(d_matrix_inverse(detector.panel(panel[i])));
      it = seen_panel.end() - 1;
    }
    Ray ray;
    out.status[i] = predict(
        beam, A, miller[i],
        seen_D[static_cast<std::size_t>(it - seen_panel.begin())], &ray);
    out.x_mm[i] = ray.x_mm;
    out.y_mm[i] = ray.y_mm;
    out.delpsi[i] = ray.delpsi;
    out.s1[i] = ray.s1;
  }
  return out;
}

Linearisation linearise_crystal(const StillsDetector &detector, const Vec3 &s0,
                                const Mat3 &U, const Mat3 &B,
                                const std::string &hall,
                                const IndexedSpots &spots,
                                const RefineOptions &options) {
  Linearisation out;
  out.used = usable(spots);
  const CellConstraint constraint = CellConstraint::from_hall(hall);
  out.parameter_names = parameter_names(constraint);
  out.parameters = start_parameters(constraint, B);
  if (out.used.empty()) {
    out.status = RefineStatus::no_reflections;
    out.message = to_string(out.status);
    return out;
  }
  try {
    Problem problem(detector, s0, U, constraint, spots, out.used, options);
    std::vector<int> all(out.parameters.size());
    for (std::size_t p = 0; p < all.size(); ++p)
      all[p] = static_cast<int>(p);
    problem.evaluate(out.parameters, all, true);
    out.residual = problem.residual();
    out.weight = problem.weight();
    out.jacobian = problem.jacobian();
  } catch (const Failure &failure) {
    out.status = failure.status;
    out.message = failure.message;
  }
  return out;
}

RefineResult refine_crystal(const StillsDetector &detector, const Vec3 &s0,
                            const Mat3 &U, const Mat3 &B,
                            const std::string &hall, const IndexedSpots &spots,
                            const RefineOptions &options) {
  RefineResult result;
  result.U = U;
  result.B = B;

  const std::vector<std::size_t> used = usable(spots);
  result.n_used = used.size();
  if (used.empty()) {
    result.status = RefineStatus::no_reflections;
    result.message = to_string(result.status);
    return result;
  }

  const CellConstraint constraint = CellConstraint::from_hall(hall);
  const std::size_t n_param = 3 + constraint.size();
  result.parameter_names = parameter_names(constraint);
  std::vector<double> x = start_parameters(constraint, B);

  try {
    Problem problem(detector, s0, U, constraint, spots, used, options);

    // Fix any parameter no spot responds to.
    std::vector<int> all(n_param);
    for (std::size_t p = 0; p < n_param; ++p)
      all[p] = static_cast<int>(p);
    problem.evaluate(x, all, true);
    std::vector<int> free;
    for (std::size_t c = 0; c < n_param; ++c)
      if (problem.column_moves_anything(c, n_param))
        free.push_back(static_cast<int>(c));
    std::vector<std::string> names;
    for (int p : free)
      names.push_back(result.parameter_names[static_cast<std::size_t>(p)]);
    result.parameter_names = names;
    const std::size_t n_free = free.size();

    std::vector<double> normal, rhs, damped, step;
    problem.evaluate(x, free, true);
    problem.normal_equations(n_free, &normal, &rhs);
    for (double v : rhs)
      if (v == 0.0)
        throw Failure{RefineStatus::indeterminate,
                      to_string(RefineStatus::indeterminate)};

    const long dof =
        static_cast<long>(problem.n_equations()) - static_cast<long>(n_free);
    if (dof < 1) {
      result.status = RefineStatus::too_few_degrees_of_freedom;
      result.message = to_string(result.status);
      result.U = problem.U();
      result.B = problem.B();
      record_predictions(problem, &result);
      return result;
    }

    double max_diagonal = normal[0];
    for (std::size_t a = 1; a < n_free; ++a)
      max_diagonal = std::max(max_diagonal, normal[a * n_free + a]);
    double mu = 1e-3 * max_diagonal;
    int nu = 2;
    int n_iterations = 0;
    // The model the refinement ends on: the last evaluated, which after a
    // run of rejected steps is the last trial, not x.
    Mat3 end_U = problem.U(), end_B = problem.B();

    while (true) {
      const double f = problem.objective();
      damped = normal;
      for (std::size_t a = 0; a < n_free; ++a)
        damped[a * n_free + a] += mu;
      step = rhs;
      if (!solve_spd(damped.data(), step.data(), n_free))
        throw Failure{RefineStatus::not_positive_definite,
                      to_string(RefineStatus::not_positive_definite)};

      JournalRow row;
      row.rmsd = problem.rmsd();
      row.objective = f;
      row.mu = mu;
      row.nu = nu;
      row.parameters = free_values(x, free);
      result.journal.push_back(row);

      const std::size_t rows = result.journal.size();
      if (rows >= 2 && rmsds_converged(result.journal[rows - 1].rmsd,
                                       result.journal[rows - 2].rmsd,
                                       options.rmsd_tolerance)) {
        result.status = RefineStatus::rmsd_converged;
        break;
      }
      if (n_iterations == options.max_iterations) {
        result.status = RefineStatus::max_iterations;
        break;
      }

      // The gradient g is -rhs.
      double expected = 0.0;
      for (std::size_t a = 0; a < n_free; ++a)
        expected += step[a] * (mu * step[a] + rhs[a]);
      expected *= 0.5;
      if (expected == 0.0) {
        result.status = RefineStatus::objective_flat;
        break;
      }

      const std::vector<double> x_old = x;
      for (std::size_t c = 0; c < n_free; ++c)
        x[static_cast<std::size_t>(free[c])] += step[c];
      ++n_iterations;
      problem.evaluate(x, free, false);
      end_U = problem.U();
      end_B = problem.B();
      const double rho = (f - problem.objective()) / expected;
      if (rho > 0.0) {
        mu *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * rho - 1.0, 3));
        nu = 2;
      } else {
        x = x_old;
        result.journal.pop_back();
        if (nu >= 8192) {
          result.status = RefineStatus::max_trial_iterations;
          break;
        }
        mu *= nu;
        nu *= 2;
      }
      problem.evaluate(x, free, true);
      end_U = problem.U();
      end_B = problem.B();
      problem.normal_equations(n_free, &normal, &rhs);
    }

    result.message = to_string(result.status);
    result.n_iterations = n_iterations;
    result.U = end_U;
    result.B = end_B;
    record_predictions(problem, &result);
  } catch (const Failure &failure) {
    result.status = failure.status;
    result.message = failure.message;
    result.U = U;
    result.B = B;
  } catch (const GeometryError &e) {
    result.status = RefineStatus::prediction_failed;
    result.message = e.what();
    result.U = U;
    result.B = B;
  }
  return result;
}

} // namespace mxi::stills

// Limited-memory BFGS minimisation (Liu & Nocedal 1989) with the More-Thuente
// line search (More & Thuente 1994), as scitbx implements it and as
// scitbx.lbfgs.run drives it. The Sauter-Poon distribution fit and the Nave
// maximum-likelihood fit both use it; following its arithmetic and its control
// flow exactly is what lets them reproduce dials' results to the last bit.
//
// Only the path those fits reach is here: no user-supplied diagonal and no
// gradient-only line search.
//
// Derived from cctbx_project (cctbx, rstbx, scitbx).
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace mxi::stills::lbfgs {

struct CoreParameters {
  std::size_t m = 5;
  std::size_t maxfev = 20;
  double gtol = 0.9;
  double xtol = 1e-16;
  double stpmin = 1e-20;
  double stpmax = 1e20;
};

struct TerminationParameters {
  bool traditional_convergence_test = true;
  double traditional_convergence_test_eps = 1e-5;
  std::size_t drop_convergence_test_n_test_points = 5;
  double drop_convergence_test_max_drop_eps = 1e-5;
  double drop_convergence_test_iteration_coefficient = 2.0;
  std::size_t min_iterations = 0;
  std::optional<std::size_t> max_iterations;
  std::optional<std::size_t> max_calls;
};

// Which line-search failures end the minimisation quietly. A failure that is
// not ignored is an error of the whole minimisation.
struct ExceptionHandlingParameters {
  bool ignore_line_search_failed_rounding_errors = true;
  bool ignore_line_search_failed_step_at_lower_bound = false;
  bool ignore_line_search_failed_step_at_upper_bound = false;
  bool ignore_line_search_failed_maxfev = false;
  bool ignore_line_search_failed_xtol = false;
  bool ignore_search_direction_not_descent = false;
};

enum class Failure {
  none,
  rounding_errors,
  step_at_upper_bound,
  step_at_lower_bound,
  maxfev,
  xtol,
  search_direction_not_descent,
  internal,
};

const char *to_string(Failure failure);

// The minimiser's state machine. run() takes the function value and gradient
// at x; it returns true when it has moved x and needs them there, and false
// when a line search has ended at x (or, on the first call, when the gradient
// is zero). A failed line search throws LineSearchFailure.
class Minimizer {
public:
  Minimizer(std::size_t n, const CoreParameters &core);

  bool run(double *x, double f, const double *g);

  std::size_t iter() const { return iter_; }
  std::size_t nfun() const { return nfun_; }

private:
  void line_search(double *x, double f, const double *g, std::size_t is0);
  static int step(double &stx, double &fx, double &dx, double &sty, double &fy,
                  double &dy, double &stp, double fp, double dp, bool &brackt,
                  double stpmin, double stpmax);

  const std::size_t n_, m_;
  const CoreParameters core_;
  int iflag_ = 0;
  bool requests_f_and_g_ = false;
  std::size_t iter_ = 0, nfun_ = 0;
  double stp_ = 0.0, stp1_ = 0.0, ys_ = 0.0;
  const double ftol_ = 0.0001;
  std::size_t point_ = 0, npt_ = 0, bound_ = 0, nfev_ = 0;
  const std::size_t ispt_, iypt_;
  int info_ = 0;
  std::vector<double> w_, diag_, wa_;

  // The line search's state between calls.
  int infoc_ = 0;
  double dginit_ = 0.0, finit_ = 0.0, dgtest_ = 0.0, width_ = 0.0,
         width1_ = 0.0;
  bool brackt_ = false, stage1_ = false;
  double stx_ = 0.0, fx_ = 0.0, dgx_ = 0.0, sty_ = 0.0, fy_ = 0.0, dgy_ = 0.0;
  double stmin_ = 0.0, stmax_ = 0.0;
};

struct LineSearchFailure {
  Failure failure;
};

// The objective: its value at x, with the gradient written to g. It may throw
// to end the minimisation; the exception is passed on unchanged.
using Objective =
    std::function<double(const std::vector<double> &x, std::vector<double> &g)>;

struct Result {
  // A line-search failure that ended the minimisation; with `raised` set it was
  // not ignored, and the minimisation is an error.
  Failure failure = Failure::none;
  bool raised = false;
  std::size_t iter = 0, nfun = 0;
};

// scitbx.lbfgs.run: minimise from x, leaving the final point in x. On an
// ignored line-search failure x is the point of lowest f evaluated.
Result run(std::vector<double> &x, const Objective &objective,
           const CoreParameters &core = {},
           const TerminationParameters &termination = {},
           const ExceptionHandlingParameters &exceptions = {});

// The convergence tests run() applies, exposed for the tests.
bool traditional_convergence_test(const std::vector<double> &x,
                                  const std::vector<double> &g, double eps);

class DropConvergenceTest {
public:
  DropConvergenceTest(std::size_t n_test_points, double max_drop_eps,
                      double iteration_coefficient);
  bool operator()(double f);

private:
  std::size_t n_test_points_;
  double max_drop_eps_, iteration_coefficient_;
  std::vector<double> x_, y_;
  double max_drop_ = 0.0, max_f_ = 0.0;
};

} // namespace mxi::stills::lbfgs

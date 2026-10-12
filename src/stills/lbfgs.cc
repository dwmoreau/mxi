// Derived from cctbx_project (cctbx, rstbx, scitbx).
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "lbfgs.hh"

#include <algorithm>
#include <cmath>

namespace mxi::stills::lbfgs {

namespace {

// Dot product and y += a x, unrolled by 5 and 4 as the Fortran ddot and daxpy
// are, so that the sums are formed in the same order.
double ddot(std::size_t n, const double *dx, const double *dy) {
  double t = 0.0;
  const std::size_t m = n % 5;
  std::size_t i = 0;
  for (; i < m; ++i)
    t += dx[i] * dy[i];
  for (; i < n; i += 5) {
    t += dx[i] * dy[i];
    t += dx[i + 1] * dy[i + 1];
    t += dx[i + 2] * dy[i + 2];
    t += dx[i + 3] * dy[i + 3];
    t += dx[i + 4] * dy[i + 4];
  }
  return t;
}

void daxpy(std::size_t n, double a, const double *dx, double *dy) {
  if (n == 0 || a == 0.0)
    return;
  for (std::size_t i = 0; i < n; ++i)
    dy[i] += a * dx[i];
}

double max3(double x, double y, double z) {
  return x < y ? (y < z ? z : y) : (x < z ? z : x);
}

[[noreturn]] void fail(Failure f) { throw LineSearchFailure{f}; }

} // namespace

const char *to_string(Failure failure) {
  switch (failure) {
  case Failure::none:
    return "none";
  case Failure::rounding_errors:
    return "rounding_errors";
  case Failure::step_at_upper_bound:
    return "step_at_upper_bound";
  case Failure::step_at_lower_bound:
    return "step_at_lower_bound";
  case Failure::maxfev:
    return "maxfev";
  case Failure::xtol:
    return "xtol";
  case Failure::search_direction_not_descent:
    return "search_direction_not_descent";
  case Failure::internal:
    return "internal";
  }
  return "unknown";
}

Minimizer::Minimizer(std::size_t n, const CoreParameters &core)
    : n_(n), m_(core.m), core_(core), ispt_(n + 2 * core.m),
      iypt_(n + 2 * core.m + n * core.m) {
  if (n_ == 0 || m_ == 0 || core.maxfev == 0 || core.gtol <= 1e-4 ||
      core.xtol < 0.0 || core.stpmin < 0.0 || core.stpmax < core.stpmin)
    fail(Failure::internal);
  w_.assign(n_ * (2 * m_ + 1) + 2 * m_, 0.0);
  diag_.assign(n_, 0.0);
  wa_.assign(n_, 0.0);
}

bool Minimizer::run(double *x, double f, const double *g) {
  bool whole_iteration = !requests_f_and_g_;
  requests_f_and_g_ = false;
  double *w = w_.data();
  if (iflag_ == 0) {
    nfun_ = 1;
    std::fill(diag_.begin(), diag_.end(), 1.0);
    for (std::size_t i = 0; i < n_; ++i)
      w[ispt_ + i] = -g[i] * diag_[i];
    const double gnorm = std::sqrt(ddot(n_, g, g));
    if (gnorm == 0.0)
      return false;
    stp1_ = 1.0 / gnorm;
    whole_iteration = true;
  }
  if (whole_iteration) {
    bound_ = iter_;
    ++iter_;
    info_ = 0;
    if (iter_ != 1) {
      if (iter_ > m_)
        bound_ = m_;
      ys_ = ddot(n_, w + iypt_ + npt_, w + ispt_ + npt_);
      const double yy = ddot(n_, w + iypt_ + npt_, w + iypt_ + npt_);
      std::fill(diag_.begin(), diag_.end(), ys_ / yy);
    }
    if (iter_ != 1) {
      std::size_t cp = point_;
      if (point_ == 0)
        cp = m_;
      w[n_ + cp - 1] = 1 / ys_;
      for (std::size_t i = 0; i < n_; ++i)
        w[i] = -g[i];
      cp = point_;
      for (std::size_t i = 0; i < bound_; ++i) {
        if (cp == 0)
          cp = m_;
        --cp;
        const double sq = ddot(n_, w + ispt_ + cp * n_, w);
        const std::size_t inmc = n_ + m_ + cp;
        w[inmc] = w[n_ + cp] * sq;
        daxpy(n_, -w[inmc], w + iypt_ + cp * n_, w);
      }
      for (std::size_t i = 0; i < n_; ++i)
        w[i] *= diag_[i];
      for (std::size_t i = 0; i < bound_; ++i) {
        const double yr = ddot(n_, w + iypt_ + cp * n_, w);
        double beta = w[n_ + cp] * yr;
        beta = w[n_ + m_ + cp] - beta;
        daxpy(n_, beta, w + ispt_ + cp * n_, w);
        ++cp;
        if (cp == m_)
          cp = 0;
      }
      std::copy(w, w + n_, w + ispt_ + point_ * n_);
    }
    stp_ = 1.0;
    if (iter_ == 1)
      stp_ = stp1_;
    std::copy(g, g + n_, w);
  }
  line_search(x, f, g, ispt_ + point_ * n_);
  if (info_ == -1) {
    iflag_ = 1;
    requests_f_and_g_ = true;
    return true;
  }
  if (info_ != 1)
    fail(Failure::internal);
  nfun_ += nfev_;
  npt_ = point_ * n_;
  for (std::size_t i = 0; i < n_; ++i) {
    w[ispt_ + npt_ + i] = stp_ * w[ispt_ + npt_ + i];
    w[iypt_ + npt_ + i] = g[i] - w[i];
  }
  ++point_;
  if (point_ == m_)
    point_ = 0;
  return false;
}

void Minimizer::line_search(double *x, double f, const double *g,
                            std::size_t is0) {
  const double *s = w_.data() + is0;
  const double gtol = core_.gtol, xtol = core_.xtol;
  const double stpmin = core_.stpmin, stpmax = core_.stpmax;
  const std::size_t maxfev = core_.maxfev;
  double &stp = stp_;
  if (info_ != -1) {
    infoc_ = 1;
    if (stp <= 0.0 || ftol_ < 0.0)
      fail(Failure::internal);
    dginit_ = 0.0;
    for (std::size_t j = 0; j < n_; ++j)
      dginit_ += g[j] * s[j];
    if (dginit_ >= 0.0)
      fail(Failure::search_direction_not_descent);
    brackt_ = false;
    stage1_ = true;
    nfev_ = 0;
    finit_ = f;
    dgtest_ = ftol_ * dginit_;
    width_ = stpmax - stpmin;
    width1_ = 2.0 * width_;
    std::copy(x, x + n_, wa_.begin());
    stx_ = 0.0;
    fx_ = finit_;
    dgx_ = dginit_;
    sty_ = 0.0;
    fy_ = finit_;
    dgy_ = dginit_;
  }
  for (;;) {
    if (info_ != -1) {
      if (brackt_) {
        stmin_ = std::min(stx_, sty_);
        stmax_ = std::max(stx_, sty_);
      } else {
        stmin_ = stx_;
        stmax_ = stp + 4.0 * (stp - stx_);
      }
      stp = std::max(stp, stpmin);
      stp = std::min(stp, stpmax);
      if ((brackt_ && (stp <= stmin_ || stp >= stmax_)) ||
          nfev_ >= maxfev - 1 || infoc_ == 0 ||
          (brackt_ && stmax_ - stmin_ <= xtol * stmax_))
        stp = stx_;
      for (std::size_t j = 0; j < n_; ++j)
        x[j] = wa_[j] + stp * s[j];
      info_ = -1;
      return;
    }
    info_ = 0;
    ++nfev_;
    double dg = 0.0;
    for (std::size_t j = 0; j < n_; ++j)
      dg += g[j] * s[j];
    const double ftest1 = finit_ + stp * dgtest_;
    if ((brackt_ && (stp <= stmin_ || stp >= stmax_)) || infoc_ == 0)
      fail(Failure::rounding_errors);
    if (stp == stpmax && f <= ftest1 && dg <= dgtest_)
      fail(Failure::step_at_upper_bound);
    if (stp == stpmin && (f > ftest1 || dg >= dgtest_))
      fail(Failure::step_at_lower_bound);
    if (nfev_ >= maxfev)
      fail(Failure::maxfev);
    if (brackt_ && stmax_ - stmin_ <= xtol * stmax_)
      fail(Failure::xtol);
    if (f <= ftest1 && std::abs(dg) <= gtol * (-dginit_)) {
      info_ = 1;
      return;
    }
    if (stage1_ && f <= ftest1 && dg >= std::min(ftol_, gtol) * dginit_)
      stage1_ = false;
    if (stage1_ && f <= fx_ && f > ftest1) {
      const double fm = f - stp * dgtest_;
      double fxm = fx_ - stx_ * dgtest_;
      double fym = fy_ - sty_ * dgtest_;
      const double dgm = dg - dgtest_;
      double dgxm = dgx_ - dgtest_;
      double dgym = dgy_ - dgtest_;
      infoc_ = step(stx_, fxm, dgxm, sty_, fym, dgym, stp, fm, dgm, brackt_,
                    stmin_, stmax_);
      fx_ = fxm + stx_ * dgtest_;
      fy_ = fym + sty_ * dgtest_;
      dgx_ = dgxm + dgtest_;
      dgy_ = dgym + dgtest_;
    } else {
      infoc_ = step(stx_, fx_, dgx_, sty_, fy_, dgy_, stp, f, dg, brackt_,
                    stmin_, stmax_);
    }
    if (brackt_) {
      if (std::abs(sty_ - stx_) >= 0.66 * width1_)
        stp = stx_ + 0.5 * (sty_ - stx_);
      width1_ = width_;
      width_ = std::abs(sty_ - stx_);
    }
  }
}

int Minimizer::step(double &stx, double &fx, double &dx, double &sty,
                    double &fy, double &dy, double &stp, double fp, double dp,
                    bool &brackt, double stpmin, double stpmax) {
  bool bound;
  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  int info = 0;
  if ((brackt && (stp <= std::min(stx, sty) || stp >= std::max(stx, sty))) ||
      dx * (stp - stx) >= 0.0 || stpmax < stpmin)
    return 0;
  sgnd = dp * (dx / std::abs(dx));
  if (fp > fx) {
    info = 1;
    bound = true;
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
    s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
    gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp < stx)
      gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx);
    if (std::abs(stpc - stx) < std::abs(stpq - stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / 2.0;
    brackt = true;
  } else if (sgnd < 0.0) {
    info = 2;
    bound = false;
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
    s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
    gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp > stx)
      gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p / q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (std::abs(stpc - stp) > std::abs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  } else if (std::abs(dp) < std::abs(dx)) {
    info = 3;
    bound = true;
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
    s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
    gamma = s * std::sqrt(std::max(0.0, (theta / s) * (theta / s) -
                                            (dx / s) * (dp / s)));
    if (stp > stx)
      gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p / q;
    if (r < 0.0 && gamma != 0.0)
      stpc = stp + r * (stx - stp);
    else if (stp > stx)
      stpc = stpmax;
    else
      stpc = stpmin;
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (brackt) {
      if (std::abs(stp - stpc) < std::abs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    } else {
      if (std::abs(stp - stpc) > std::abs(stp - stpq))
        stpf = stpc;
      else
        stpf = stpq;
    }
  } else {
    info = 4;
    bound = false;
    if (brackt) {
      theta = 3.0 * (fp - fy) / (sty - stp) + dy + dp;
      s = max3(std::abs(theta), std::abs(dy), std::abs(dp));
      gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
      if (stp > sty)
        gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  stpf = std::min(stpmax, stpf);
  stpf = std::max(stpmin, stpf);
  stp = stpf;
  if (brackt && bound) {
    if (sty > stx)
      stp = std::min(stx + 0.66 * (sty - stx), stp);
    else
      stp = std::max(stx + 0.66 * (sty - stx), stp);
  }
  return info;
}

bool traditional_convergence_test(const std::vector<double> &x,
                                  const std::vector<double> &g, double eps) {
  const double xnorm = std::sqrt(ddot(x.size(), x.data(), x.data()));
  const double gnorm = std::sqrt(ddot(g.size(), g.data(), g.data()));
  return gnorm <= eps * std::max(1.0, xnorm);
}

DropConvergenceTest::DropConvergenceTest(std::size_t n_test_points,
                                         double max_drop_eps,
                                         double iteration_coefficient)
    : n_test_points_(n_test_points), max_drop_eps_(max_drop_eps),
      iteration_coefficient_(iteration_coefficient) {}

bool DropConvergenceTest::operator()(double f) {
  if (!x_.empty())
    max_drop_ = std::max(max_drop_, y_.back() - f);
  max_f_ = std::max(max_f_, std::abs(f));
  x_.push_back(static_cast<double>(x_.size() + 1));
  y_.push_back(f);
  if (x_.size() < n_test_points_)
    return false;
  if (max_f_ == 0.0)
    return true;
  const std::size_t first = y_.size() - n_test_points_;
  // scitbx's linear_regression over the last points, y scaled by max |f|.
  const double *x = x_.data() + first;
  std::vector<double> y(n_test_points_);
  for (std::size_t i = 0; i < n_test_points_; ++i)
    y[i] = y_[first + i] / max_f_;
  double min_x = x[0], max_x = x[0], min_y = y[0], max_y = y[0];
  double sum_x = x[0], sum_x2 = x[0] * x[0], sum_y = y[0], sum_xy = x[0] * y[0];
  for (std::size_t i = 1; i < n_test_points_; ++i) {
    if (min_x > x[i])
      min_x = x[i];
    if (max_x < x[i])
      max_x = x[i];
    if (min_y > y[i])
      min_y = y[i];
    if (max_y < y[i])
      max_y = y[i];
    sum_x += x[i];
    sum_x2 += x[i] * x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
  }
  const double sum_weights = static_cast<double>(n_test_points_);
  double slope = 0.0;
  if (min_y != max_y) {
    const double dx = std::max(std::abs(min_x - sum_x / sum_weights),
                               std::abs(max_x - sum_x / sum_weights));
    const double dy = std::max(std::abs(min_y - sum_y / sum_weights),
                               std::abs(max_y - sum_y / sum_weights));
    if (dy != 0.0) {
      if (dx < dy * 1e-15)
        throw LineSearchFailure{Failure::internal};
      const double d = sum_weights * sum_x2 - sum_x * sum_x;
      if (d != 0.0)
        slope = (sum_weights * sum_xy - sum_x * sum_y) / d;
    }
  }
  const double sliding_tolerance =
      max_drop_ * max_drop_eps_ *
      std::pow(static_cast<double>(x_.size()), iteration_coefficient_);
  return -slope * max_f_ <= sliding_tolerance;
}

Result run(std::vector<double> &x, const Objective &objective,
           const CoreParameters &core, const TerminationParameters &term,
           const ExceptionHandlingParameters &handling) {
  const std::size_t n = x.size();
  Minimizer minimizer(n, core);
  DropConvergenceTest drop(
      std::max(term.drop_convergence_test_n_test_points, term.min_iterations),
      term.drop_convergence_test_max_drop_eps,
      term.drop_convergence_test_iteration_coefficient);
  std::vector<double> g(n), x_min;
  double f_min = 0.0;
  bool have_min = false;
  Result result;
  try {
    for (;;) {
      const double f = objective(x, g);
      if (!have_min) {
        if (!term.traditional_convergence_test)
          drop(f);
        f_min = f;
        x_min = x;
        have_min = true;
      } else if (f_min > f) {
        f_min = f;
        x_min = x;
      }
      if (minimizer.run(x.data(), f, g.data()))
        continue;
      if (term.traditional_convergence_test) {
        if (minimizer.iter() >= term.min_iterations &&
            traditional_convergence_test(x, g,
                                         term.traditional_convergence_test_eps))
          break;
      } else if (drop(f)) {
        break;
      }
      if (term.max_iterations && minimizer.iter() >= *term.max_iterations)
        break;
      if (term.max_calls && minimizer.nfun() > *term.max_calls)
        break;
      if (!minimizer.run(x.data(), f, g.data()))
        break;
    }
  } catch (const LineSearchFailure &e) {
    if (have_min)
      x = x_min;
    result.failure = e.failure;
    // scitbx.lbfgs exception_handling_parameters.filter: 0 and -1 end quietly.
    switch (e.failure) {
    case Failure::rounding_errors:
      result.raised = !handling.ignore_line_search_failed_rounding_errors;
      break;
    case Failure::step_at_lower_bound:
      result.raised = !(traditional_convergence_test(x, g, 1e-5) ||
                        handling.ignore_line_search_failed_step_at_lower_bound);
      break;
    case Failure::step_at_upper_bound:
      result.raised = !handling.ignore_line_search_failed_step_at_upper_bound;
      break;
    case Failure::maxfev:
      result.raised = !handling.ignore_line_search_failed_maxfev;
      break;
    case Failure::xtol:
      result.raised = !handling.ignore_line_search_failed_xtol;
      break;
    case Failure::search_direction_not_descent:
      result.raised = !(traditional_convergence_test(x, g, 1e-5) ||
                        handling.ignore_search_direction_not_descent);
      break;
    default:
      result.raised = true;
    }
  }
  result.iter = minimizer.iter();
  result.nfun = minimizer.nfun();
  return result;
}

} // namespace mxi::stills::lbfgs

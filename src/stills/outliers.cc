// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "outliers.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace mxi::stills {

namespace {

constexpr double kPi = 3.141592653589793;
constexpr double kFitFraction = 0.4;
constexpr std::size_t kCurvePoints = 500;

double rayleigh_cdf(double x, double sigma) {
  const double x_sigma = x / sigma;
  return 1.0 - std::exp(-0.5 * x_sigma * x_sigma);
}

double rayleigh_cdf_derivative(double x, double sigma) {
  const double x_sigma = x / sigma;
  const double pdf = (x_sigma / sigma) * std::exp(-0.5 * x_sigma * x_sigma);
  return -(x / sigma) * pdf;
}

} // namespace

const char *to_string(SauterPoonStatus status) {
  switch (status) {
  case SauterPoonStatus::ok:
    return "ok";
  case SauterPoonStatus::equal_radii:
    return "equal_radii";
  case SauterPoonStatus::too_few:
    return "too_few";
  case SauterPoonStatus::fit_failed:
    return "fit_failed";
  }
  return "unknown";
}

SauterPoonResult sauter_poon(const std::vector<double> &x_obs,
                             const std::vector<double> &y_obs,
                             const std::vector<double> &x_calc,
                             const std::vector<double> &y_calc) {
  const std::size_t n = x_obs.size();
  SauterPoonResult out;

  std::vector<double> radius(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double dx = x_obs[i] - x_calc[i];
    const double dy = y_obs[i] - y_calc[i];
    radius[i] = std::sqrt(dx * dx + dy * dy);
  }
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return radius[a] < radius[b];
  });
  std::vector<double> dr(n);
  for (std::size_t k = 0; k < n; ++k)
    dr[k] = radius[order[k]];
  for (std::size_t k = 1; k < n; ++k)
    if (dr[k] == dr[k - 1]) {
      out.status = SauterPoonStatus::equal_radii;
      return out;
    }

  std::vector<double> cdf(n);
  for (std::size_t k = 0; k < n; ++k)
    cdf[k] = static_cast<double>(k) / static_cast<double>(n);

  const auto limit =
      static_cast<std::size_t>(kFitFraction * static_cast<double>(n));
  out.n_fit = limit;
  if (limit < 4) {
    out.status = SauterPoonStatus::too_few;
    return out;
  }

  // No point of the fitted range has a CDF above one half, so the start is
  // its last point.
  std::vector<double> sigma{dr[limit - 1]};
  out.sigma_start = sigma[0];
  const auto objective = [&](const std::vector<double> &x,
                             std::vector<double> &g) {
    double f = 0.0, g0 = 0.0;
    for (std::size_t k = 0; k < limit; ++k) {
      const double d = rayleigh_cdf(dr[k], x[0]) - cdf[k];
      f += std::sqrt(d * d);
      g0 = g0 + d * rayleigh_cdf_derivative(dr[k], x[0]);
    }
    g[0] = 2.0 * g0;
    return f;
  };
  out.minimiser = lbfgs::run(sigma, objective);
  out.sigma = sigma[0];
  if (out.minimiser.raised) {
    out.status = SauterPoonStatus::fit_failed;
    return out;
  }

  out.sd = std::sqrt((4.0 - kPi) / (2.0) * out.sigma * out.sigma);
  out.first_rejected = n;
  bool found = false;
  for (std::size_t i = 0; i < kCurvePoints && !found; ++i) {
    const double mx =
        static_cast<double>(i) / static_cast<double>(kCurvePoints);
    const double my = rayleigh_cdf(mx, out.sigma);
    for (std::size_t j = 1; j < n; ++j) {
      if (my >= cdf[j - 1] && my < cdf[j] && out.sd <= dr[j] - mx &&
          dr[j - 1] - mx > 0.0) {
        out.first_rejected = j - 1;
        found = true;
        break;
      }
    }
  }

  out.kept.assign(n, true);
  for (std::size_t k = out.first_rejected; k < n; ++k)
    out.kept[order[k]] = false;
  return out;
}

} // namespace mxi::stills

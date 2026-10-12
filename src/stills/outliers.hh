// Sauter-Poon outlier rejection (Sauter & Poon 2010) as dials' stills indexer
// runs it: the radial distances between observed and calculated spot
// positions are taken as a Rayleigh sample, a Rayleigh CDF is fitted to the
// smallest 40 %, and every spot from the first one lying further than the
// distribution's standard deviation beyond the fitted curve outward is
// rejected.
//
// The steps and their arithmetic are those of rstbx's
// find_outliers_from_matches on a first pass, with every spot good:
//   dr      sqrt(dx^2 + dy^2) per spot; two equal values are an error
//   cdf     the sorted dr against k / N, k = 0..N-1
//   fit     sigma of 1 - exp(-r^2 / 2 sigma^2) on the first int(0.4 N) points,
//           from sigma = the last of them, minimising sum |Phi(r_k) - k/N|
//           with the gradient 2 sum (Phi - k/N) dPhi/dsigma, by L-BFGS
//   bar     sd = sqrt((4 - pi) / 2) sigma; scanning the fitted curve at
//           r = i / 500 mm, the first sorted spot j - 1 whose step of the
//           empirical CDF contains the curve's value, and which lies beyond
//           the curve, with spot j at least sd beyond it
//   reject  sorted spots j - 1 .. N - 1
//
// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <cstddef>
#include <vector>

#include "lbfgs.hh"

namespace mxi::stills {

enum class SauterPoonStatus {
  ok,
  equal_radii, // two spots at exactly the same distance from prediction
  too_few,     // fewer than 4 points to fit (N < 10)
  fit_failed,  // the minimiser ended on a line-search error not ignored
};

const char *to_string(SauterPoonStatus status);

struct SauterPoonResult {
  SauterPoonStatus status = SauterPoonStatus::ok;
  std::vector<bool> kept; // input order; empty unless status is ok
  double sigma = 0.0;     // the fitted Rayleigh parameter, mm
  double sigma_start = 0.0;
  double sd = 0.0;
  std::size_t n_fit = 0;
  // The sorted position from which spots are rejected; N when none are.
  std::size_t first_rejected = 0;
  lbfgs::Result minimiser;
};

// Observed and calculated positions in millimetres, in the reflection table's
// order. dials forms them as pixel positions times panel 0's pixel size.
SauterPoonResult sauter_poon(const std::vector<double> &x_obs,
                             const std::vector<double> &y_obs,
                             const std::vector<double> &x_calc,
                             const std::vector<double> &y_calc);

} // namespace mxi::stills

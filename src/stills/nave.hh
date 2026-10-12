// The mosaicity and domain size of a still's crystal, fitted to the
// Delta-Psi of its indexed reflections, as dials' NaveParameters fits them
// (Sauter et al. 2014, eq. 16 and appendix B):
//
//   |Delta-Psi| <= d / D + eta / 2 is the band a reflection of spacing d falls
//   in, for a crystal of mosaic blocks of size D and full angular spread eta.
//
// The steps and their arithmetic are dials':
//   seed      the reflections sorted by 2-theta into min(max(3, N/25), 50)
//             bins of equal count; a straight-line fit of each bin's largest
//             |Delta-Psi| against its mean d gives D = 1/(2 slope), raised to
//             at least 3 V^(1/3), and eta = |2 intercept|
//   ML        the maximum-likelihood alpha = 2/D and eta of serialtbx's
//             max_like.minimizer, from that seed, by L-BFGS on log alpha and
//             log eta, with its gradient as coded there
//   accept    |Delta-Psi| < d alpha / 2 + eta / 2, in degrees
//   crystal   half-mosaicity eta / 2 and domain size 2/alpha, each widened by
//             1.4 (the half-mosaicity multiplied, the domain size divided)
//   reject    the crystal when its Ewald proximal volume exceeds the limit
//
// With safe_likelihood (the default), wherever serialtbx's exponentials would
// pass 1e100 the same likelihood and gradient are evaluated in log-sum-exp
// form, so a line search step far outside the data gets a large finite value
// instead of ending the fit as diverged. Elsewhere the evaluation is
// serialtbx's. Off, the fit is dials' to the last bit, failures included.
//
// d, 2-theta and the cell volume come from the crystal; cell_d_spacings and
// friends compute them from B.
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

#include <array>
#include <cstddef>
#include <vector>

#include "lbfgs.hh"
#include "linalg.hh"

namespace mxi::stills {

struct NaveOptions {
  double ewald_proximal_volume_max = 0.0025;
  double ewald_proximity_resolution_cutoff = 2.0; // Angstrom
  double model_expansion_factor = 1.4;
  bool safe_likelihood = true;
};

enum class NaveStatus {
  ok,
  too_few,          // fewer than three reflections: a bin would be empty
  singular_seed,    // the seed's line fit has no solution or zero slope
  zero_seed_eta,    // the seed's intercept is zero, so log eta is undefined
  outside_domain,   // an exponent beyond the double range
  diverged,         // the model grew past 1e100
  zero_probability, // a reflection with probability zero under the model
  fit_failed,       // the minimiser ended on a line-search error not ignored
  ewald_volume,     // the crystal's Ewald proximal volume exceeds the limit
};

const char *to_string(NaveStatus status);
bool failed(NaveStatus status);

struct NaveResult {
  NaveStatus status = NaveStatus::ok;
  // The seed.
  double eta_seed = 0.0, domain_size_seed = 0.0;
  // The maximum-likelihood alpha = 2/D (per Angstrom) and full mosaicity eta
  // (radians).
  double alpha = 0.0, eta = 0.0;
  std::vector<bool> accepted;
  // The crystal's values, widened by the expansion factor.
  double half_mosaicity_deg = 0.0, domain_size_ang = 0.0;
  double ewald_volume = 0.0;
  lbfgs::Result minimiser;
};

// delpsi in radians, d in Angstrom and two_theta in degrees per reflection;
// the cell volume in cubic Angstrom.
NaveResult nave_parameters(const std::vector<double> &delpsi,
                           const std::vector<double> &d,
                           const std::vector<double> &two_theta,
                           double cell_volume, double wavelength,
                           const NaveOptions &options = {});

// d = 1 / |B h| per reflection.
std::vector<double> cell_d_spacings(const Mat3 &B,
                                    const std::vector<std::array<int, 3>> &h);
// 2 asin(lambda / 2d), in degrees.
std::vector<double> two_theta_deg(const std::vector<double> &d,
                                  double wavelength);
// The direct cell volume, 1 / det B.
double cell_volume(const Mat3 &B);

} // namespace mxi::stills

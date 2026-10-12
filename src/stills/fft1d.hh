// The fft1d basis-vector search as dials runs it (Steller, Bolotovsky &
// Rossmann 1997; Sauter, Grosse-Kunstleve & Adams 2004): the reciprocal
// lattice points projected onto directions over a hemisphere, a histogram of
// each projection Fourier-transformed, and the directions with the strongest
// periodicity, refined and de-duplicated, returned as real-space vectors.
//
// rstbx's version is followed, not the papers', where they differ (the grid,
// the adaptive two-pass search, the peak rule, the clustering): the
// differences are listed in the R2 spec, section 3. Ties between directions of
// equal strength are common, and their order reaches the output, so every
// step follows rstbx's arithmetic and orderings exactly, including the
// unstable std::sort of its final ranking.
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
#include <optional>
#include <stdexcept>
#include <vector>

#include "linalg.hh"

namespace mxi::stills {

// rstbx refuses seven or fewer reciprocal lattice points (an assertion).
struct TooFewSpots : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

// What a search did, for comparing it with rstbx step by step.
struct Fft1dTrace {
  double characteristic_grid = 0, incr = 0;
  bool adaptive = false;
  std::size_t n_quick = 0, n_quick_winners = 0, n_directions = 0;
  // The directions before refinement: kval, uc_length, dvec.
  std::vector<std::array<double, 5>> unrefined;
  // After refinement, before setSolutions.
  std::vector<std::array<double, 5>> refined;
};

// dials' FFT1D.find_basis_vectors: 0 to 30 real-space vectors (Angstrom, lab
// frame), strongest first. characteristic_grid is dials' fft1d parameter of
// that name (radians); unset, it is derived from the data. trace, if given,
// records the intermediate steps.
std::vector<Vec3>
fft1d_basis_vectors(const std::vector<Vec3> &rlp, double max_cell,
                    std::optional<double> characteristic_grid = std::nullopt,
                    Fft1dTrace *trace = nullptr);

// rstbx's SimpleSamplerTool.construct_hemisphere_grid(sampling): the
// directions (as unit vectors) in its order.
std::vector<Vec3> hemisphere_grid(double sampling);

} // namespace mxi::stills

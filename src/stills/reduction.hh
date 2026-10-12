// Niggli reduction as cctbx runs it: Krivy & Gruber (1976) with the
// epsilon comparisons of Grosse-Kunstleve, Sauter & Adams (2004),
// epsilon = V^(1/3) * relative_epsilon, compared against squared lengths.
//
// Derived from cctbx_project (cctbx, rstbx, scitbx) and dxtbx.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.
// dxtbx: Copyright (c) 2019 Diamond Light Source, Lawrence Berkeley National
// Laboratory and the Science and Technology Facilities Council. All rights
// reserved. BSD-3-Clause; see third_party/dxtbx/LICENSE.txt.

#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>

#include "stills/cell.hh"

namespace mxi::stills {

using IntMat3 = std::array<int, 9>;

struct NiggliLimitExceeded : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct NiggliResult {
  // The product of the step matrices: new basis vector j is
  // sum_i r_inv(i, j) old vector i.
  IntMat3 r_inv{};
  std::size_t n_iterations = 0;
};

// Throws NiggliLimitExceeded when the step count reaches iteration_limit, and
// CellError where a step leaves a non-positive squared length (cctbx asserts).
NiggliResult niggli_reduction(const UnitCell &cell,
                              double relative_epsilon = 1e-5,
                              std::size_t iteration_limit = 1000);

// The matrix dxtbx's Crystal::change_basis multiplies A by for the operator
// cctbx returns as the change of basis to the Niggli cell: (r_inv^-1)^T.
Mat3 niggli_change_of_basis(const NiggliResult &r);

} // namespace mxi::stills

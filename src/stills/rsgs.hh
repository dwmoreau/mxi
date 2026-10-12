// dials' real-space grid search (Gildea et al. 2014, section 3, eqs. 10-11):
// each target cell length, laid along each direction of a hemisphere grid, is
// scored by sum_j cos(2 pi r_j . v) over the reciprocal lattice points, and the
// strongest vectors, merged where they agree within 10 % in length and 5
// degrees in direction, are the candidate basis vectors.
//
// dials scores with numpy: the dot products come from OpenBLAS's dgemm, which
// on this machine forms each as fma(a2, b2, fma(a1, b1, a0 * b0)), and each
// row's sum is numpy's pairwise summation. Both are reproduced, so the scores,
// their order and the vectors chosen are dials' to the last bit; a plain
// sequential sum would change the chosen vectors in about 2 % of calls.
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

#include "linalg.hh"

namespace mxi::stills {

class RealSpaceGridSearch {
public:
  // dials' RealSpaceGridSearch(max_cell, target_unit_cell): the search
  // vectors, each distinct cell length (ascending) times each grid direction.
  explicit RealSpaceGridSearch(const std::array<double, 6> &target_cell,
                               double characteristic_grid = 0.02);

  // find_basis_vectors: up to max_vectors, strongest first.
  std::vector<Vec3> find_basis_vectors(const std::vector<Vec3> &rlp,
                                       std::size_t max_vectors = 30) const;

  // score_vectors: one score per search vector, in their order.
  std::vector<double> scores(const std::vector<Vec3> &rlp) const;

  const std::vector<Vec3> &search_vectors() const { return sv_; }

private:
  std::vector<Vec3> sv_;
};

} // namespace mxi::stills

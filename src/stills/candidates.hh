// Candidate crystal models from basis vectors, as dials'
// candidate_orientation_matrices makes them (Gildea et al. 2014, section 2.4),
// and global index assignment, as dials' AssignIndicesGlobal does it.
//
// Derived from dials.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "stills/cell.hh"

namespace mxi::stills {

// dials' candidate_orientation_matrices as the generator it is: each call to
// next() does only the work needed to yield the next candidate, so a Niggli
// reduction that would fail on a later triplet is never reached by a caller
// that stops early.
class CandidateStream {
public:
  explicit CandidateStream(std::vector<Vec3> basis_vectors);
  // The next candidate, or nothing once the triplets are exhausted. Throws
  // NiggliLimitExceeded as dials raises DialsIndexError.
  std::optional<CrystalModel> next();

private:
  std::vector<Vec3> basis_vectors_;
  std::vector<std::array<int, 3>> triplets_;
  std::size_t next_ = 0;
};

// Every triplet (i < j < k) of the first 100 vectors, ordered by
// i^2 + j^2 + k^2 (stable), that passes the 20 degree angle tests, made
// obtuse and right-handed, reduced to its Niggli cell, and kept if
// V > abc / 100. Throws NiggliLimitExceeded as dials raises DialsIndexError.
std::vector<CrystalModel>
candidate_orientation_matrices(const std::vector<Vec3> &basis_vectors);

struct Assignment {
  std::vector<std::array<int, 3>> miller_index;
  // The crystal each row was assigned to, or -1.
  std::vector<int> id;
};

// Rows with id == -1 and 1/|rlp| > d_min (d_min <= 0: no limit) are indexed
// against every crystal: h = A^-1 rlp rounded half away from zero, the crystal
// with the smallest |h - round(h)|^2 kept if that is at most tolerance^2 and
// h is not (0, 0, 0). Of rows given the same h on the same crystal, the
// nearest to its integer keeps it, ties going to the earlier row. Every other
// row comes back as (0, 0, 0) and -1, whatever its id was. All rows are taken
// to be on one image, so the rotation-angle test never applies.
Assignment assign_indices_global(const std::vector<Vec3> &rlp,
                                 const std::vector<int> &id,
                                 const std::vector<Mat3> &A, double d_min,
                                 double tolerance = 0.3);

} // namespace mxi::stills

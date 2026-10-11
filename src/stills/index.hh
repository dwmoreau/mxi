// Index one still: the entry point dials.stills_process calls.

#pragma once

#include <string>

#include "geometry.hh"

namespace mxi::stills {

enum class IndexStatus { indexed, not_indexed };

// A still that does not index is an ordinary outcome, returned with the
// reason. Exceptions are kept for inputs that break the contract.
struct IndexResult {
  IndexStatus status = IndexStatus::not_indexed;
  std::string reason;
  std::size_t n_spots = 0;
};

// Maps the spots to reciprocal space; the search that would follow does not
// exist yet, so every still comes back not indexed.
IndexResult index_still(const StillsDetector &detector, const Vec3 &s0,
                        double wavelength, const Spots &spots);

} // namespace mxi::stills

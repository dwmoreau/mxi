#include "index.hh"

namespace mxi::stills {

IndexResult index_still(const StillsDetector &detector, const Vec3 &s0,
                        double wavelength, const Spots &spots) {
  const MappedSpots mapped = map_spots(detector, s0, wavelength, spots);
  IndexResult result;
  result.status = IndexStatus::not_indexed;
  result.reason = "mxi_stills has no indexer yet";
  result.n_spots = mapped.rlp.size();
  return result;
}

} // namespace mxi::stills

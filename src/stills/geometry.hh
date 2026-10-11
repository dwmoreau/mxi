// The geometry of a still: a fixed multi-panel detector, and the map from a
// spot's pixel centroid to millimetres on its panel, to the diffracted beam
// s1, and to the reciprocal lattice point s1 - s0.
//
// Each step follows the arithmetic dials and dxtbx use, in the same order, so
// that the results agree to the last bit rather than to a tolerance:
//
//   px -> mm   dxtbx ParallaxCorrectedPxMmStrategy::to_millimeter, i.e.
//              parallax_correction_inv evaluated at the uncorrected position;
//   variance   the pixel variance times the pixel size squared, with the
//              parallax Jacobian ignored (dials centroid_px_to_mm_panel);
//   mm -> lab  VirtualPanel::get_lab_coord, d * (x, y, 1);
//   s1, rlp    map_centroids_to_reciprocal_space with no goniometer:
//              s1 = lab / |lab| * (1 / wavelength), rlp = s1 - s0.
//
// Every panel's vectors are in the laboratory frame, with any detector
// hierarchy already composed. dxtbx serialises panels in their local frame, so
// these must come from the live panel, never from an .expt.

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "linalg.hh"

namespace mxi::stills {

struct GeometryError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct StillsPanel {
  Vec3 fast{1.0, 0.0, 0.0};
  Vec3 slow{0.0, -1.0, 0.0};
  Vec3 origin{0.0, 0.0, -100.0};
  double pixel_size[2] = {0.1, 0.1};

  // Linear attenuation coefficient (per mm) and sensor thickness (mm) of the
  // parallax correction. Both zero means a plain pixel-size conversion.
  double mu = 0.0;
  double t0 = 0.0;

  bool parallax() const { return mu != 0.0 || t0 != 0.0; }

  // The mean depth of absorption along the unit ray u: Winter et al. (2018)
  // eqn (7). Throws where dxtbx asserts, for mu <= 0 or a ray that does not
  // meet the sensor from the front.
  double attenuation_length(const Vec3 &u) const;

  std::pair<double, double> px_to_mm(double px_fast, double px_slow) const;

  Vec3 lab_coord(double mm_fast, double mm_slow) const {
    return Vec3{fast.x * mm_fast + slow.x * mm_slow + origin.x,
                fast.y * mm_fast + slow.y * mm_slow + origin.y,
                fast.z * mm_fast + slow.z * mm_slow + origin.z};
  }
};

// The detector is the same for every still a process handles, so it is built
// once and kept.
struct StillsDetector {
  std::vector<StillsPanel> panels;

  std::size_t size() const { return panels.size(); }
  const StillsPanel &panel(std::int64_t index) const;
};

// A still's strong spots, one entry per spot.
struct Spots {
  std::vector<std::int64_t> panel;
  std::vector<double> x_px, y_px;         // centroid
  std::vector<double> var_x_px, var_y_px; // centroid variance

  std::size_t size() const { return panel.size(); }
};

struct MappedSpots {
  std::vector<double> x_mm, y_mm;
  std::vector<double> var_x_mm, var_y_mm;
  std::vector<Vec3> s1;  // per Angstrom
  std::vector<Vec3> rlp; // per Angstrom, laboratory frame
};

// s0 and the wavelength as the beam model gives them; neither is derived from
// the other here.
MappedSpots map_spots(const StillsDetector &detector, const Vec3 &s0,
                      double wavelength, const Spots &spots);

} // namespace mxi::stills

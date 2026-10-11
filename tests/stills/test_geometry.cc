// The stills geometry against dials' own numbers for real spots, on detectors
// of 256, 32, 8 and 1 panels.

#include <cmath>
#include <string>

#include "check.hh"
#include "real_stills.hh"
#include "stills/geometry.hh"
#include "stills/index.hh"

using namespace mxi::stills;

namespace {

// TASKS R5's agreement gate. The corpus check finds every value bit-identical
// with the compiler used for dials; this leaves room for another libm's exp.
constexpr double kGate = 1e-12;

StillsDetector detector_of(const real_stills::RealStill &r, bool local_frame,
                           bool parallax) {
  StillsDetector d;
  for (std::size_t i = 0; i < r.panels.size(); ++i) {
    const auto &p = r.panels[i];
    StillsPanel q;
    q.fast = local_frame ? r.local[i].fast : p.fast;
    q.slow = local_frame ? r.local[i].slow : p.slow;
    q.origin = local_frame ? r.local[i].origin : p.origin;
    q.pixel_size[0] = p.pixel_size[0];
    q.pixel_size[1] = p.pixel_size[1];
    q.mu = parallax ? p.mu : 0.0;
    q.t0 = parallax ? p.t0 : 0.0;
    d.panels.push_back(q);
  }
  return d;
}

Spots spots_of(const real_stills::RealStill &r) {
  Spots s;
  for (const auto &x : r.spots) {
    s.panel.push_back(x.panel);
    s.x_px.push_back(x.x_px);
    s.y_px.push_back(x.y_px);
    s.var_x_px.push_back(x.var_x_px);
    s.var_y_px.push_back(x.var_y_px);
  }
  return s;
}

// The largest difference from dials over every mapped quantity.
double worst_difference(const real_stills::RealStill &r, const MappedSpots &m) {
  double worst = 0.0;
  auto take = [&](double a, double b) {
    worst = std::fmax(worst, std::abs(a - b));
  };
  for (std::size_t i = 0; i < r.spots.size(); ++i) {
    const auto &e = r.spots[i];
    take(m.x_mm[i], e.x_mm);
    take(m.y_mm[i], e.y_mm);
    take(m.var_x_mm[i], e.var_x_mm);
    take(m.var_y_mm[i], e.var_y_mm);
    for (std::size_t k = 0; k < 3; ++k) {
      take(m.s1[i][k], e.s1[k]);
      take(m.rlp[i][k], e.rlp[k]);
    }
  }
  return worst;
}

void check_against_dials(const real_stills::RealStill &r,
                         std::size_t min_panels) {
  check::is_true(r.spots.size() >= 3, r.name + ": fixture has spots");
  check::is_true(r.panels.size() >= min_panels,
                 r.name + ": spots on " + std::to_string(r.panels.size()) +
                     " panels, want " + std::to_string(min_panels));
  const MappedSpots m =
      map_spots(detector_of(r, false, true), r.s0, r.wavelength, spots_of(r));
  check::equal(static_cast<long long>(m.rlp.size()),
               static_cast<long long>(r.spots.size()), r.name + ": mapped");
  check::close(worst_difference(r, m), 0.0, kGate,
               r.name + ": largest difference from dials");
}

bool refuses(const StillsDetector &d, const real_stills::RealStill &r,
             const Spots &s) {
  try {
    map_spots(d, r.s0, r.wavelength, s);
  } catch (const GeometryError &) {
    return true;
  }
  return false;
}

} // namespace

TEST(stills_geometry_cytochrome_matches_dials) {
  check_against_dials(real_stills::real_cytochrome(), 8);
}

TEST(stills_geometry_psii_matches_dials) {
  check_against_dials(real_stills::real_psii(), 7);
}

TEST(stills_geometry_sacla_matches_dials) {
  check_against_dials(real_stills::real_sacla(), 6);
}

TEST(stills_geometry_lysozyme_matches_dials) {
  check_against_dials(real_stills::real_lysozyme(), 1);
}

// The vectors an .expt serialises are in each panel's local frame. Inside a
// hierarchy they are wrong by far more than the gate, or describe a panel the
// ray cannot meet, which the parallax correction refuses.
TEST(stills_geometry_local_frame_vectors_are_wrong) {
  for (const auto &r : {real_stills::real_cytochrome(),
                        real_stills::real_psii(), real_stills::real_sacla()}) {
    const StillsDetector d = detector_of(r, true, true);
    const Spots s = spots_of(r);
    const bool caught =
        refuses(d, r, s) ||
        worst_difference(r, map_spots(d, r.s0, r.wavelength, s)) > 1e-3;
    check::is_true(caught, r.name + ": local-frame vectors not caught");
  }
}

// Leaving the parallax correction out moves the spots by far more than the
// gate.
TEST(stills_geometry_without_parallax_is_wrong) {
  for (const auto &r :
       {real_stills::real_cytochrome(), real_stills::real_psii(),
        real_stills::real_sacla(), real_stills::real_lysozyme()}) {
    const double worst =
        worst_difference(r, map_spots(detector_of(r, false, false), r.s0,
                                      r.wavelength, spots_of(r)));
    check::is_true(worst > 1e-4, r.name + ": dropping parallax changed only " +
                                     std::to_string(worst));
  }
}

TEST(stills_geometry_refuses_bad_input) {
  const auto r = real_stills::real_psii();
  const StillsDetector d = detector_of(r, false, true);

  Spots s = spots_of(r);
  s.panel[0] = static_cast<std::int64_t>(d.size());
  check::is_true(refuses(d, r, s), "a panel past the end is refused");

  s = spots_of(r);
  s.panel[0] = -1;
  check::is_true(refuses(d, r, s), "a negative panel is refused");

  s = spots_of(r);
  s.var_y_px.pop_back();
  check::is_true(refuses(d, r, s), "columns of different lengths are refused");

  // A thickness with no attenuation coefficient, which dxtbx asserts against.
  StillsDetector zero_mu = d;
  for (auto &p : zero_mu.panels)
    p.mu = 0.0;
  check::is_true(refuses(zero_mu, r, spots_of(r)),
                 "parallax with mu = 0 is refused");
}

TEST(stills_index_stub_returns_not_indexed) {
  const auto r = real_stills::real_cytochrome();
  const IndexResult result =
      index_still(detector_of(r, false, true), r.s0, r.wavelength, spots_of(r));
  check::is_true(result.status == IndexStatus::not_indexed, "not indexed");
  check::is_true(!result.reason.empty(), "with a reason");
  check::equal(static_cast<long long>(result.n_spots),
               static_cast<long long>(r.spots.size()), "every spot mapped");
}

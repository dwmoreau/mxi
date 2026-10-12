// The indexer on whole corpus frames, against dials' result: the same refined
// reflections, with dials' indices up to one operator of the lattice
// (candidates that are one lattice in several settings tie on RMSD, and either
// setting may win), the crystal to well within the agreement R6 found for
// refinement, and the Nave model to its fit's flat direction.

#include <algorithm>
#include <cmath>
#include <string>

#include "check.hh"
#include "stills/indexer.hh"
#include "stills/real_frames.hh"

using namespace mxi::stills;

namespace {

StillsDetector detector_of(const real_frames::RealFrame &f) {
  StillsDetector d;
  for (const auto &p : f.panels) {
    StillsPanel q;
    q.fast = p.fast;
    q.slow = p.slow;
    q.origin = p.origin;
    q.pixel_size[0] = p.pixel_size[0];
    q.pixel_size[1] = p.pixel_size[1];
    q.mu = p.mu;
    q.t0 = p.t0;
    d.panels.push_back(q);
  }
  return d;
}

StrongSpots spots_of(const real_frames::RealFrame &f) {
  StrongSpots s;
  for (const auto &p : f.spots) {
    s.spots.panel.push_back(p.panel);
    s.spots.x_px.push_back(p.x);
    s.spots.y_px.push_back(p.y);
    s.spots.var_x_px.push_back(p.var_x);
    s.spots.var_y_px.push_back(p.var_y);
    s.flags.push_back(p.flags);
  }
  return s;
}

IndexerSettings settings_of(const real_frames::RealFrame &f) {
  IndexerSettings s;
  s.unit_cell = f.unit_cell;
  s.space_group_hall = f.hall;
  s.methods = f.methods;
  s.refine_candidates_with_known_symmetry = f.known_symmetry_path;
  s.d_min_start = f.d_min_start;
  s.nave.safe_likelihood = false;
  return s;
}

// The integer matrix M, elements in {-1, 0, 1} and determinant +-1, with
// M h = h_dials for every row, if there is one.
bool setting_operator(const std::vector<std::array<int, 3>> &h,
                      const std::vector<std::array<int, 3>> &h_dials,
                      std::array<int, 9> *out) {
  std::array<int, 9> m{};
  for (int code = 0; code < 19683; ++code) {
    int c = code;
    for (int &e : m) {
      e = c % 3 - 1;
      c /= 3;
    }
    const int det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
                    m[1] * (m[3] * m[8] - m[5] * m[6]) +
                    m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (det != 1 && det != -1)
      continue;
    bool all = true;
    for (std::size_t r = 0; all && r < h.size(); ++r)
      for (int i = 0; all && i < 3; ++i)
        all = m[3 * i] * h[r][0] + m[3 * i + 1] * h[r][1] +
                  m[3 * i + 2] * h[r][2] ==
              h_dials[r][static_cast<std::size_t>(i)];
    if (all) {
      *out = m;
      return true;
    }
  }
  return false;
}

void check_against_dials(const real_frames::RealFrame &f) {
  const StillsIndexer indexer(detector_of(f), settings_of(f));
  const IndexedStill r = indexer.index(f.s0, f.wavelength, spots_of(f));
  check::is_true(r.indexed, f.tag + " indexes: " + r.reason);
  check::is_true(r.row == f.dials_rows, f.tag + ": dials' refined reflections");
  std::array<int, 9> m{};
  check::is_true(setting_operator(r.miller_index, f.dials_miller, &m),
                 f.tag + ": dials' Miller indices, up to one operator");

  // A h = A_dials M h for every h, so A = A_dials M.
  const mxi::Mat3 A = r.crystal.A();
  double largest = 0, worst = 0;
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) {
      double expected = 0;
      for (std::size_t k = 0; k < 3; ++k)
        expected += f.dials_A[3 * i + k] * m[3 * k + j];
      largest = std::max(largest, std::abs(expected));
      worst = std::max(worst, std::abs(A.m[3 * i + j] - expected));
    }
  check::is_true(worst <= 1e-6 * largest,
                 f.tag + ": A within 1e-6 of its largest element");
  check::close(r.half_mosaicity_deg, f.dials_half_mosaicity_deg,
               1e-2 * f.dials_half_mosaicity_deg, f.tag + ": half-mosaicity");
  check::close(r.domain_size_ang, f.dials_domain_size_ang,
               1e-2 * f.dials_domain_size_ang, f.tag + ": domain size");
  check::equal(static_cast<long long>(r.x_cal_px.size()),
               static_cast<long long>(r.row.size()), f.tag + ": predictions");
}

} // namespace

TEST(stills_indexer_matches_dials_on_the_p1_path) {
  check_against_dials(real_frames::real_lysozyme_frame());
}

TEST(stills_indexer_matches_dials_on_the_known_symmetry_path) {
  check_against_dials(real_frames::real_sacla_frame());
}

TEST(stills_indexer_reports_a_still_it_cannot_index) {
  real_frames::RealFrame f = real_frames::real_lysozyme_frame();
  f.spots.resize(5);
  const StillsIndexer indexer(detector_of(f), settings_of(f));
  const IndexedStill r = indexer.index(f.s0, f.wavelength, spots_of(f));
  check::is_true(!r.indexed, "five spots do not index");
  check::equal(static_cast<long long>(r.methods.size()), 2,
               "both methods are tried");
  check::is_true(r.row.empty(), "and nothing is returned");
}

TEST(stills_indexer_with_no_spots_ends_before_any_method) {
  const real_frames::RealFrame f = real_frames::real_lysozyme_frame();
  const StillsIndexer indexer(detector_of(f), settings_of(f));
  const IndexedStill r = indexer.index(f.s0, f.wavelength, StrongSpots{});
  check::is_true(!r.indexed, "not indexed");
  check::is_true(r.reason == "No reflections left to index!", r.reason);
  check::is_true(r.methods.empty(), "no method ran");
}

TEST(stills_indexer_max_cell_is_from_the_primitive_target) {
  const real_frames::RealFrame f = real_frames::real_lysozyme_frame();
  const StillsIndexer indexer(detector_of(f), settings_of(f));
  check::close(indexer.max_cell(), 1.3 * 78.0, 1e-12, "1.3 times the longest");
}

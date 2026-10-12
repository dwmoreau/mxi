// The real-space grid search against dials' own results on real calls, and
// its search vectors' layout.

#include <string>
#include <vector>

#include "check.hh"
#include "real_rsgs.hh"
#include "stills/fft1d.hh"
#include "stills/rsgs.hh"

using namespace mxi::stills;

TEST(stills_rsgs_search_vectors_are_length_major) {
  // Two distinct lengths (78 twice, 37): each length ascending times every
  // direction of the 0.02 rad grid (15955 of them).
  const RealSpaceGridSearch s({78, 78, 37, 90, 90, 90});
  const auto &sv = s.search_vectors();
  const auto grid = hemisphere_grid(0.02);
  check::equal(static_cast<long long>(grid.size()), 15955, "grid size");
  check::equal(static_cast<long long>(sv.size()), 2 * 15955, "two lengths");
  for (std::size_t i : {std::size_t{0}, std::size_t{100}})
    for (std::size_t c = 0; c < 3; ++c) {
      check::is_true(sv[i][c] == grid[i][c] * 37, "37 first");
      check::is_true(sv[15955 + i][c] == grid[i][c] * 78, "then 78");
    }
}

TEST(stills_rsgs_matches_dials) {
  for (const auto &c : real_rsgs::cases) {
    std::vector<mxi::Vec3> rlp;
    for (const auto &v : c.rlp)
      rlp.push_back({v[0], v[1], v[2]});
    const RealSpaceGridSearch s(c.target_cell, c.characteristic_grid);
    const auto got =
        s.find_basis_vectors(rlp, static_cast<std::size_t>(c.max_vectors));
    check::equal(static_cast<long long>(got.size()),
                 static_cast<long long>(c.basis_vectors.size()),
                 c.name + ": number of vectors");
    for (std::size_t i = 0; i < got.size(); ++i)
      for (std::size_t j = 0; j < 3; ++j)
        check::is_true(got[i][j] == c.basis_vectors[i][j],
                       c.name + ": vector " + std::to_string(i));
  }
}

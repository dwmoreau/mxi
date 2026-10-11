// Niggli reduction, candidate generation and global index assignment against
// cctbx's and dials' own results on real cells and calls.

#include <string>
#include <vector>

#include "check.hh"
#include "real_candidates.hh"
#include "stills/candidates.hh"
#include "stills/cell.hh"
#include "stills/reduction.hh"

using namespace mxi::stills;

namespace {

mxi::Vec3 vec3(const std::array<double, 3> &v) { return {v[0], v[1], v[2]}; }

} // namespace

TEST(stills_niggli_matches_cctbx) {
  for (const auto &c : real_candidates::niggli_cases) {
    const auto r = niggli_reduction(UnitCell::from_parameters(c.params));
    for (int i = 0; i < 9; ++i)
      check::equal(r.r_inv[i], c.r_inv[i],
                   c.name + ": r_inv[" + std::to_string(i) + "]");
    check::equal(static_cast<long long>(r.n_iterations), c.n_iterations,
                 c.name + ": iterations");
  }
}

TEST(stills_niggli_iteration_limit) {
  const auto &c = real_candidates::niggli_cases.back();
  check::is_true(c.n_iterations > 1, c.name + " needs more than one step");
  bool thrown = false;
  try {
    niggli_reduction(UnitCell::from_parameters(c.params), 1e-5,
                     static_cast<std::size_t>(c.n_iterations - 1));
  } catch (const NiggliLimitExceeded &) {
    thrown = true;
  }
  check::is_true(thrown, "the limit one short of the step count raises");
  niggli_reduction(UnitCell::from_parameters(c.params), 1e-5,
                   static_cast<std::size_t>(c.n_iterations));
}

TEST(stills_candidates_match_dials) {
  for (const auto &c : real_candidates::candidate_cases) {
    std::vector<mxi::Vec3> bv;
    for (const auto &v : c.basis_vectors)
      bv.push_back(vec3(v));
    const auto models = candidate_orientation_matrices(bv);
    check::equal(static_cast<long long>(models.size()), c.n_candidates,
                 c.name + ": number of candidates");
    for (std::size_t k = 0; k < c.first_A.size(); ++k) {
      const mxi::Mat3 A = models[k].A();
      for (int i = 0; i < 9; ++i)
        check::is_true(A.m[i] == c.first_A[k][i],
                       c.name + ": candidate " + std::to_string(k) + " A[" +
                           std::to_string(i) + "]");
      for (int i = 0; i < 6; ++i)
        check::is_true(models[k].cell.params[i] == c.first_cell[k][i],
                       c.name + ": candidate " + std::to_string(k) + " cell[" +
                           std::to_string(i) + "]");
    }
  }
}

TEST(stills_candidates_need_three_vectors) {
  check::equal(
      static_cast<long long>(
          candidate_orientation_matrices({{10, 0, 0}, {0, 10, 0}}).size()),
      0, "two vectors give no candidates");
}

TEST(stills_assign_matches_dials) {
  for (const auto &c : real_candidates::assign_cases) {
    std::vector<mxi::Vec3> rlp;
    for (const auto &v : c.rlp)
      rlp.push_back(vec3(v));
    mxi::Mat3 A;
    A.m = c.A;
    const auto out = assign_indices_global(rlp, c.id_in, {A}, c.d_min);
    for (std::size_t i = 0; i < rlp.size(); ++i) {
      for (int j = 0; j < 3; ++j)
        check::equal(out.miller_index[i][j], c.miller_index[i][j],
                     c.name + ": row " + std::to_string(i) + " miller");
      check::equal(out.id[i], c.id_out[i],
                   c.name + ": row " + std::to_string(i) + " id");
    }
  }
}

TEST(stills_assign_duplicates_keep_the_nearer) {
  // Two spots near h = (1, 0, 0) on a cubic 10 A cell: the nearer keeps it,
  // and with equal distances the earlier row does.
  const mxi::Mat3 A{0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1};
  auto out =
      assign_indices_global({{0.11, 0, 0}, {0.101, 0, 0}}, {-1, -1}, {A}, 0.0);
  check::equal(out.id[0], -1, "farther row loses");
  check::equal(out.id[1], 0, "nearer row keeps h");
  out = assign_indices_global({{0.1, 0.01, 0}, {0.1, -0.01, 0}}, {-1, -1}, {A},
                              0.0);
  check::equal(out.id[0], 0, "tie: earlier row keeps h");
  check::equal(out.id[1], -1, "tie: later row loses");
  out = assign_indices_global({{0.1, 0, 0}, {0.2, 0, 0}}, {0, -1}, {A}, 0.0);
  check::equal(out.id[0], -1, "a row already indexed is not reassigned");
  check::equal(out.miller_index[1][0], 2, "an unindexed row is");
}

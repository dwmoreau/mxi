// Space-group types, lattice symmetry and dials' known-symmetry filter against
// cctbx's and dials' own results on reference groups, corpus targets and real
// candidate cells.

#include <string>

#include "check.hh"
#include "real_symmetry.hh"
#include "stills/known_symmetry.hh"
#include "stills/sgtbx.hh"

using namespace mxi::stills;
namespace sg = mxi::stills::sgtbx;

namespace {

sg::RtMx rt(const real_symmetry::Rt &v) {
  sg::IntMat3 r{};
  for (int i = 0; i < 9; ++i)
    r[static_cast<std::size_t>(i)] = v[static_cast<std::size_t>(i)];
  return sg::RtMx(sg::RotMx(r, v[9]), sg::TrVec(v[10], v[11], v[12], v[13]));
}

void same_rt(const sg::RtMx &a, const real_symmetry::Rt &b,
             const std::string &what) {
  for (int i = 0; i < 9; ++i)
    check::equal(a.r()[static_cast<std::size_t>(i)],
                 b[static_cast<std::size_t>(i)], what + " r");
  check::equal(a.r().den(), b[9], what + " r_den");
  for (int i = 0; i < 3; ++i)
    check::equal(a.t()[static_cast<std::size_t>(i)],
                 b[static_cast<std::size_t>(10 + i)], what + " t");
  check::equal(a.t().den(), b[13], what + " t_den");
}

} // namespace

TEST(stills_space_group_type_matches_cctbx) {
  for (const auto &c : real_symmetry::type_cases) {
    const std::string name = std::to_string(c.number) + " (" + c.setting + ")";
    const sg::SpaceGroup g = sg::SpaceGroup::from_hall(c.hall).change_basis(
        sg::ChangeOfBasisOp(rt(c.cb_c), rt(c.cb_c_inv)));
    const sg::SpaceGroupType t(g);
    check::equal(t.number(), c.expected_number, name + ": number");
    same_rt(t.cb_op().c(), c.expected_c, name + ": cb_op c");
    same_rt(t.cb_op().c_inv(), c.expected_c_inv, name + ": cb_op c_inv");
  }
}

TEST(stills_cached_space_group_results_match_cctbx) {
  // Each case twice: the second pass reads change_basis, the type and
  // find_affine from the per-thread caches, and must still give cctbx's answers
  // and the same group and matrices.
  for (int pass = 0; pass < 2; ++pass)
    for (const auto &c : real_symmetry::type_cases) {
      const std::string name = std::to_string(c.number) + " (" + c.setting +
                               "), pass " + std::to_string(pass);
      const sg::SpaceGroup h = sg::SpaceGroup::from_hall(c.hall);
      const sg::ChangeOfBasisOp cb(rt(c.cb_c), rt(c.cb_c_inv));
      const sg::SpaceGroup g = h.change_basis(cb);
      sg::Key first, again;
      g.append_key(first);
      h.change_basis(cb).append_key(again);
      check::is_true(first == again, name + ": change_basis");
      const sg::SpaceGroupType t(g);
      check::equal(t.number(), c.expected_number, name + ": number");
      same_rt(t.cb_op().c(), c.expected_c, name + ": cb_op c");
      same_rt(t.cb_op().c_inv(), c.expected_c_inv, name + ": cb_op c_inv");
      if (c.expected_number >= 3 && c.expected_number <= 15) {
        sg::Key affine_first, affine_again;
        for (const sg::RtMx &m : sg::find_affine(g))
          sg::append_key(affine_first, m);
        for (const sg::RtMx &m : sg::find_affine(g))
          sg::append_key(affine_again, m);
        check::is_true(!affine_first.empty() && affine_first == affine_again,
                       name + ": find_affine");
      }
    }
}

TEST(stills_reference_symbols_parse) {
  for (int n = 1; n <= 230; ++n) {
    const sg::SpaceGroup g =
        sg::SpaceGroup::from_hall(sg::hall_symbol_table(n), true);
    const sg::SpaceGroupType t(g);
    check::equal(t.number(), n, "reference setting " + std::to_string(n));
    check::is_true(t.cb_op().is_identity_op(), "reference setting " +
                                                   std::to_string(n) +
                                                   " is its own reference");
  }
}

TEST(stills_known_symmetry_targets_match_dials) {
  for (const auto &c : real_symmetry::target_cases) {
    const KnownSymmetry k = KnownSymmetry::make(
        UnitCell::from_parameters(c.cell), sg::SpaceGroup::from_hall(c.hall),
        c.rel, c.ang, c.max_delta);
    check::is_true(k.target_unit_cell.has_value(), c.dataset + ": target cell");
    for (std::size_t i = 0; i < 6; ++i)
      check::is_true(k.target_unit_cell->params[i] == c.target_unit_cell[i],
                     c.dataset + ": target cell[" + std::to_string(i) + "]");
    check::is_true(k.volume_window.has_value(), c.dataset + ": volume window");
    check::is_true(k.volume_window->first == c.volume_window[0] &&
                       k.volume_window->second == c.volume_window[1],
                   c.dataset + ": volume window");
    check::is_true(k.target_bravais == c.bravais, c.dataset + ": Bravais type");
  }
}

TEST(stills_known_symmetry_filter_matches_dials) {
  for (const auto &t : real_symmetry::target_cases) {
    const KnownSymmetry k = KnownSymmetry::make(
        UnitCell::from_parameters(t.cell), sg::SpaceGroup::from_hall(t.hall),
        t.rel, t.ang, t.max_delta);
    int n = 0;
    for (const auto &c : real_symmetry::model_cases) {
      if (c.dataset != t.dataset)
        continue;
      const std::string name = c.dataset + " model " + std::to_string(n++);
      const auto r =
          filter_known_symmetry(UnitCell::from_parameters(c.cell), k);
      check::equal(static_cast<int>(r.outcome), c.outcome, name + ": outcome");
      check::is_true(r.best.has_value() == (c.max_angular_difference >= 0),
                     name + ": best subgroup");
      if (!r.best)
        continue;
      check::is_true(r.best->max_angular_difference == c.max_angular_difference,
                     name + ": max_angular_difference");
      for (std::size_t i = 0; i < 6; ++i)
        check::is_true(r.best->best_subsym.cell.params[i] == c.best_cell[i],
                       name + ": best cell[" + std::to_string(i) + "]");
    }
  }
}

TEST(stills_lattice_symmetry_of_a_cubic_cell) {
  // A cubic P cell has the 9 two-folds of m-3m, which generate 432 (24
  // operations); a 6-degree shear leaves only the two-folds it keeps.
  const UnitCell cubic = UnitCell::from_parameters({50, 50, 50, 90, 90, 90});
  const sg::SpaceGroup g = lattice_symmetry_group(cubic, 5);
  check::equal(static_cast<long long>(g.n_smx()), 24, "cubic: 432");
  check::is_true(lattice_symmetry_find_max_delta(cubic, g) == 0.0,
                 "cubic: zero obliquity");
  const UnitCell sheared = UnitCell::from_parameters({50, 50, 50, 96, 90, 90});
  const sg::SpaceGroup gs = lattice_symmetry_group(sheared, 5);
  check::is_true(gs.n_smx() < 24, "sheared: lower symmetry");
}

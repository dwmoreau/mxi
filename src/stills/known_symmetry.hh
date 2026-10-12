// dials' known-symmetry filter and the cctbx machinery under it: the fast
// minimum (Buerger) reduction of Grosse-Kunstleve, Sauter & Adams (2004), the
// lattice symmetry from the 81 two-folds of a reduced cell (Le Page 1982,
// Lebedev, Vagin & Murshudov 2006), its subgroups, the best cell, and
// dials' find_matching_symmetry and filter_known_symmetry on top.
//
// Follows cctbx_project cctbx/uctbx/fast_minimum_reduction.h,
// cctbx/sgtbx/lattice_symmetry.cpp, cctbx/sgtbx/subgroups.py,
// cctbx/crystal/__init__.py (symmetry) and cctbx/crystal/find_best_cell.py,
// and dials/algorithms/indexing/symmetry.py and
// basis_vector_search/combinations.py. cctbx Copyright (c) 2006 - 2026, The
// Regents of the University of California, through Lawrence Berkeley National
// Laboratory. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "stills/cell.hh"
#include "stills/sgtbx.hh"

namespace mxi::stills {

// uctbx::unit_cell::change_basis: the metric transformed by the rotation part
// of c_inv, as doubles.
UnitCell change_basis(const UnitCell &cell, const sgtbx::ChangeOfBasisOp &cb);
UnitCell change_basis(const UnitCell &cell, const sgtbx::RotMx &c_inv_r);

// uctbx::unit_cell::is_similar_to, lengths relative and angles absolute.
bool is_similar_to(const UnitCell &a, const UnitCell &b,
                   double relative_length_tolerance = 0.01,
                   double absolute_angle_tolerance = 1.0);

// cctbx's fast_minimum_reduction (iteration limit 100, significant-change
// multiplier 16, 2 steps without significant change).
struct MinimumReductionLimitExceeded : std::runtime_error {
  using std::runtime_error::runtime_error;
};
sgtbx::IntMat3 fast_minimum_reduction(const UnitCell &cell,
                                      std::size_t iteration_limit = 100,
                                      double multiplier = 16,
                                      std::size_t min_n_no_change = 2);

// cctbx's crystal.symmetry: a cell and a space group. Building one tidies the
// group, checks the cell against it (assert_is_compatible_unit_cell, 1% and
// 1 degree) and replaces the cell by its average over the group's operations.
struct CrystalSymmetry {
  UnitCell cell;
  sgtbx::SpaceGroup group;

  static CrystalSymmetry make(const UnitCell &cell, sgtbx::SpaceGroup group,
                              bool assert_is_compatible_unit_cell = true);
  CrystalSymmetry change_basis(const sgtbx::ChangeOfBasisOp &cb) const;
  sgtbx::ChangeOfBasisOp change_of_basis_op_to_minimum_cell() const;
  sgtbx::ChangeOfBasisOp change_of_basis_op_to_niggli_cell() const;
  sgtbx::ChangeOfBasisOp change_of_basis_op_to_reference_setting() const;
  // find_best_cell(self, angular_tolerance=3, best_monoclinic_beta).cb_op()
  sgtbx::ChangeOfBasisOp
  change_of_basis_op_to_best_cell(bool best_monoclinic_beta = true) const;
};

// A cell the space group does not fit (crystal.symmetry's assertion).
struct IncompatibleUnitCell : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// cctbx's lattice_symmetry::group and find_max_delta.
sgtbx::SpaceGroup lattice_symmetry_group(const UnitCell &reduced_cell,
                                         double max_delta = 3,
                                         bool enforce_max_delta = true);
double lattice_symmetry_find_max_delta(const UnitCell &reduced_cell,
                                       const sgtbx::SpaceGroup &group);

// cctbx's subgroups(parent).groups_parent_setting(). parent is tidied first,
// as space_group_info(group=parent) does.
std::vector<sgtbx::SpaceGroup>
subgroups_parent_setting(sgtbx::SpaceGroup parent);

// dials' metric_supergroup: the group and its Euclidean normalizer's
// additional generators, made acentric.
sgtbx::SpaceGroup metric_supergroup(sgtbx::SpaceGroup group);

// dials' calc_acentric_subgroups: the subgroups of the lattice group whose
// reference setting has the target Bravais type, ordered by order_z * 1000 +
// space-group number descending, each with its metric supergroup and its
// change of basis to the reference setting. The subgroups are tidied, as
// dials' call to info() leaves them.
struct AcentricSubgroups {
  std::vector<sgtbx::SpaceGroup> subgroups, supergroups;
  std::vector<sgtbx::ChangeOfBasisOp> cb_ops;
};
AcentricSubgroups calc_acentric_subgroups(sgtbx::SpaceGroup lattice_group,
                                          const std::string &target_bravais);

// dials' find_matching_symmetry(unit_cell, None, max_delta,
// best_monoclinic_beta, target_bravais_str).
struct MatchingSymmetry {
  double max_angular_difference = 0;
  CrystalSymmetry subsym, ref_subsym, best_subsym;
};
std::optional<MatchingSymmetry>
find_matching_symmetry(const UnitCell &unit_cell, double max_delta,
                       const std::string &target_bravais,
                       bool best_monoclinic_beta = true);

// The target side of dials' filter_known_symmetry: the target cell in its
// reference setting's best cell, the volume window, and the target Bravais
// type of the primitive setting's reference group. Built once per target.
struct KnownSymmetry {
  std::optional<UnitCell> target_unit_cell;
  std::optional<std::pair<double, double>> volume_window;
  std::string target_bravais;
  double relative_length_tolerance = 0.1;
  double absolute_angle_tolerance = 5;
  double max_delta = 5;

  // The target symmetry dials passes: its cell as it holds it (already
  // averaged over the group, so not averaged again here), or none, and its
  // space group.
  static KnownSymmetry make(const std::optional<UnitCell> &cell,
                            const sgtbx::SpaceGroup &group,
                            double relative_length_tolerance = 0.1,
                            double absolute_angle_tolerance = 5,
                            double max_delta = 5);
};

enum class KnownSymmetryOutcome {
  accepted,
  outside_volume_window,
  no_matching_subgroup,
  not_similar_to_target,
};

struct KnownSymmetryResult {
  KnownSymmetryOutcome outcome;
  std::optional<MatchingSymmetry> best;
};

// One candidate through filter_known_symmetry.
KnownSymmetryResult filter_known_symmetry(const UnitCell &cell,
                                          const KnownSymmetry &target);

} // namespace mxi::stills

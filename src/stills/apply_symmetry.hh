// dials' SymmetryHandler: the target symmetry a user gives (a space group, a
// cell or both), and apply_symmetry, which puts a candidate crystal into the
// target's lattice symmetry. A candidate's metric subgroups come from rstbx's
// iotbx_converter (cctbx's metric_subgroups on the unreduced cell, sorted by
// rstbx's bestcmp as CPython's list.sort orders them); the best is chosen, the
// orientation is moved to its reference setting and constrained by rstbx's
// s_minimizer (an L-BFGS fit of the nine direct-matrix elements to the crystal
// system's metric relations), and the result is returned in the input setting
// with the change of basis to the user's target.
//
// Follows dials/algorithms/indexing/symmetry.py (SymmetryHandler),
// cctbx_project rstbx/dps_core/lepage.py, rstbx/symmetry/subgroup.py,
// cctbx/sgtbx/lattice_symmetry.py (metric_subgroups),
// rstbx/dps_core/constrainment.py and cctbx/uctbx/crystal_orientation.cpp.
// Derived from dials, cctbx_project (cctbx, rstbx, scitbx) and dxtbx.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.
// dxtbx: Copyright (c) 2019 Diamond Light Source, Lawrence Berkeley National
// Laboratory and the Science and Technology Facilities Council. All rights
// reserved. BSD-3-Clause; see third_party/dxtbx/LICENSE.txt.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "stills/cell.hh"
#include "stills/known_symmetry.hh"
#include "stills/sgtbx.hh"

namespace mxi::stills {

// One metric subgroup of a candidate cell, as iotbx_converter lists them.
struct MetricSubgroup {
  CrystalSymmetry subsym, supersym, ref_subsym, best_subsym;
  sgtbx::ChangeOfBasisOp cb_op_inp_best;
  double max_angular_difference = 0;
  std::string bravais;
  // The crystal system s_minimizer constrains to.
  std::string system;
  // The space-group number of subsym (bestcmp's tie-break).
  int number = 0;
};

// rstbx's iotbx_converter(unit_cell, max_delta, target_bravais=...): the
// metric subgroups of the target Bravais type, in its sorted order.
std::vector<MetricSubgroup> iotbx_converter(const UnitCell &unit_cell,
                                            double max_delta,
                                            const std::string &target_bravais);

// rstbx's crystal_orientation.constrain(system) through s_minimizer: the
// direct matrix (rows a, b, c) fitted to the system's metric relations,
// starting from direct.
Mat3 constrain_orientation(const Mat3 &direct, const std::string &system);

struct AppliedSymmetry {
  CrystalModel model;
  sgtbx::SpaceGroup space_group;
  sgtbx::ChangeOfBasisOp cb_op;
};

class SymmetryHandler {
public:
  // dials' SymmetryHandler(unit_cell, space_group, max_delta).
  SymmetryHandler(const std::optional<UnitCell> &unit_cell,
                  const std::optional<sgtbx::SpaceGroup> &space_group,
                  double max_delta = 5);

  // The target in its reference setting, as filter_known_symmetry gets it.
  const std::optional<UnitCell> &reference_cell() const {
    return reference_cell_;
  }
  const std::optional<sgtbx::SpaceGroup> &reference_group() const {
    return reference_group_;
  }
  const std::string &target_bravais() const { return target_bravais_; }
  // The target cell in its primitive setting (target_symmetry_primitive), set
  // when a cell is given: what max_cell and the grid search are taken from.
  const std::optional<UnitCell> &primitive_cell() const {
    return primitive_cell_;
  }
  KnownSymmetry known_symmetry(double relative_length_tolerance,
                               double absolute_angle_tolerance,
                               double max_delta) const;

  // dials' apply_symmetry: nothing when no metric subgroup has the target's
  // Bravais type. Throws IncompatibleUnitCell where dials asserts that the
  // constrained cell fits the target group.
  std::optional<AppliedSymmetry>
  apply_symmetry(const CrystalModel &crystal) const {
    return apply_symmetry(crystal.A(), crystal.cell);
  }
  // The same from what dials reads of the crystal: get_A() and its cell.
  std::optional<AppliedSymmetry> apply_symmetry(const Mat3 &A,
                                                const UnitCell &cell) const;

private:
  double max_delta_;
  std::optional<UnitCell> inp_cell_;
  std::optional<sgtbx::SpaceGroup> inp_group_;
  std::optional<UnitCell> reference_cell_;
  std::optional<UnitCell> primitive_cell_;
  std::optional<sgtbx::SpaceGroup> reference_group_;
  sgtbx::ChangeOfBasisOp cb_op_inp_ref_, cb_op_ref_inp_;
  std::optional<sgtbx::ChangeOfBasisOp> cb_op_inp_best_;
  sgtbx::SpaceGroup target_sg_ref_;
  std::string target_bravais_;
};

} // namespace mxi::stills

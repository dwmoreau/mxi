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

// Follows dials, rstbx and cctbx (see apply_symmetry.hh). cctbx Copyright (c)
// 2006 - 2026, The Regents of the University of California, through Lawrence
// Berkeley National Laboratory. BSD-3-Clause-LBNL; see
// third_party/cctbx/LICENSE.txt.

#include "stills/apply_symmetry.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "stills/lbfgs.hh"

namespace mxi::stills {

using sgtbx::ChangeOfBasisOp;
using sgtbx::RtMx;
using sgtbx::SpaceGroup;
using sgtbx::SpaceGroupType;

namespace {

// The reference group of a space-group number, tidied, as
// space_group_info(number=n).group() holds it.
SpaceGroup reference_group_of(int number) {
  SpaceGroup g = SpaceGroup::from_hall(sgtbx::hall_symbol_table(number), true);
  g.make_tidy();
  return g;
}

const char *crystal_system_name(sgtbx::CrystalSystem x) {
  switch (x) {
  case sgtbx::CrystalSystem::triclinic:
    return "triclinic";
  case sgtbx::CrystalSystem::monoclinic:
    return "monoclinic";
  case sgtbx::CrystalSystem::orthorhombic:
    return "orthorhombic";
  case sgtbx::CrystalSystem::tetragonal:
    return "tetragonal";
  case sgtbx::CrystalSystem::trigonal:
    return "trigonal";
  case sgtbx::CrystalSystem::hexagonal:
    return "hexagonal";
  case sgtbx::CrystalSystem::cubic:
    return "cubic";
  default:
    throw sgtbx::Error("crystal system undefined");
  }
}

// The reference settings of the 14 centric Bravais holohedries,
// bravais_types.centric.
bool is_bravais_holohedry(int number) {
  static const int numbers[] = {2,   10,  12,  47,  65,  69,  71,
                                123, 139, 166, 191, 221, 225, 229};
  return std::find(std::begin(numbers), std::end(numbers), number) !=
         std::end(numbers);
}

// change_of_basis_op(str(cb_op), r_den=144, t_den=144): the same operator
// with both denominators 144.
ChangeOfBasisOp at_144(const ChangeOfBasisOp &cb) {
  return ChangeOfBasisOp(cb.c().new_denominators(144, 144));
}

// rstbx's bestcmp: obliquities within 0.01 degree are equal and then the
// higher space-group number comes first; otherwise the larger obliquity.
int bestcmp(const MetricSubgroup &a, const MetricSubgroup &b) {
  if (std::fabs(a.max_angular_difference - b.max_angular_difference) <
      1. - 0.99) {
    if (a.number > b.number)
      return -1;
    if (a.number == b.number)
      return 0;
    return 1;
  }
  if (a.max_angular_difference > b.max_angular_difference)
    return -1;
  return 1;
}

// CPython 3.11's list.sort for fewer than 64 items: the first run (reversed if
// strictly descending), then binary insertion of the rest. bestcmp is not a
// consistent order, so the result is the algorithm's, and is reproduced here.
template <typename T, typename Less>
void cpython_list_sort(std::vector<T> &v, Less islt) {
  const std::size_t n = v.size();
  if (n < 2)
    return;
  if (n >= 64)
    throw std::logic_error("cpython_list_sort: only runs below 64 are ported");
  std::size_t run = 2;
  bool descending = false;
  if (islt(v[1], v[0])) {
    descending = true;
    for (std::size_t i = 2; i < n; ++i, ++run)
      if (!islt(v[i], v[i - 1]))
        break;
  } else {
    for (std::size_t i = 2; i < n; ++i, ++run)
      if (islt(v[i], v[i - 1]))
        break;
  }
  if (descending)
    std::reverse(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(run));
  for (std::size_t start = run; start < n; ++start) {
    T pivot = v[start];
    std::size_t l = 0, r = start;
    do {
      const std::size_t p = l + ((r - l) >> 1);
      if (islt(pivot, v[p]))
        r = p;
      else
        l = p + 1;
    } while (l < r);
    for (std::size_t k = start; k > l; --k)
      v[k] = v[k - 1];
    v[l] = pivot;
  }
}

// crystal_orientation(direct, basis_type.direct).unit_cell(): the metric of
// the rows of the direct matrix it gets back by inverting its reciprocal one.
UnitCell orientation_unit_cell(const Mat3 &astar) {
  const Mat3 direct = inverse_scitbx(astar);
  const Vec3 p0 = direct.row(0), p1 = direct.row(1), p2 = direct.row(2);
  return UnitCell::from_metric(
      {p0.dot(p0), p1.dot(p1), p2.dot(p2), p0.dot(p1), p2.dot(p0), p1.dot(p2)});
}

Mat3 mat3_of(const std::vector<double> &x) {
  Mat3 m;
  for (int i = 0; i < 9; ++i)
    m.m[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)];
  return m;
}

// s_minimizer.internals(): f and its gradient for one crystal system, with
// A..F from the metric of Orientation(x, direct).unit_cell().
double s_minimizer_functional(const std::string &system,
                              const std::vector<double> &x,
                              std::vector<double> &g) {
  const Sym6 mm = orientation_unit_cell(inverse_scitbx(mat3_of(x))).metric;
  const double A = mm[0], B = mm[1], C = mm[2], D = mm[5], E = mm[4], F = mm[3];
  const double z = 0;
  const double dA[9] = {2. * x[0], 2. * x[1], 2. * x[2], z, z, z, z, z, z};
  const double dB[9] = {z, z, z, 2. * x[3], 2. * x[4], 2. * x[5], z, z, z};
  const double dC[9] = {z, z, z, z, z, z, 2. * x[6], 2. * x[7], 2. * x[8]};
  const double dD[9] = {z, z, z, x[6], x[7], x[8], x[3], x[4], x[5]};
  const double dE[9] = {x[6], x[7], x[8], z, z, z, x[0], x[1], x[2]};
  const double dF[9] = {x[3], x[4], x[5], x[0], x[1], x[2], z, z, z};
  g.assign(9, 0.0);
  if (system == "monoclinic") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * D * dD[i] + 2. * F * dF[i];
    return D * D + F * F;
  }
  if (system == "orthorhombic") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * D * dD[i] + 2. * E * dE[i] + 2. * F * dF[i];
    return D * D + E * E + F * F;
  }
  if (system == "tetragonal") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * D * dD[i] + 2. * E * dE[i] + 2. * F * dF[i] +
             2. * (A - B) * (dA[i] - dB[i]);
    return D * D + E * E + F * F + (A - B) * (A - B);
  }
  if (system == "cubic") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * D * dD[i] + 2. * E * dE[i] + 2. * F * dF[i] +
             2. * (A - B) * (dA[i] - dB[i]) + 2. * (A - C) * (dA[i] - dC[i]);
    return D * D + E * E + F * F + (A - B) * (A - B) + (A - C) * (A - C);
  }
  if (system == "rhombohedral") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * (A - B) * (dA[i] - dB[i]) + 2. * D * dD[i] + 2. * E * dE[i] +
             2. * (F + A / 2.) * (dF[i] + dA[i] / 2.);
    return (A - B) * (A - B) + D * D + E * E + (F + A / 2.) * (F + A / 2.);
  }
  if (system == "hexagonal") {
    for (int i = 0; i < 9; ++i)
      g[i] = 2. * D * dD[i] + 2. * E * dE[i] + 2. * (A - B) * (dA[i] - dB[i]) +
             2. * (F + A / 2.) * (dF[i] + dA[i] / 2.);
    return D * D + E * E + (A - B) * (A - B) + (F + A / 2.) * (F + A / 2.);
  }
  if (system == "triclinic") {
    g.assign(9, 0.0);
    return 0.0;
  }
  throw std::invalid_argument("s_minimizer: no crystal system " + system);
}

} // namespace

std::vector<MetricSubgroup> iotbx_converter(const UnitCell &unit_cell,
                                            double max_delta,
                                            const std::string &target_bravais) {
  // crystal.symmetry(unit_cell, space_group_symbol="P 1"); force_minimum is
  // False, so the "minimum" symmetry is the input itself.
  const CrystalSymmetry minimum =
      CrystalSymmetry::make(unit_cell, SpaceGroup());
  const ChangeOfBasisOp cb_op_inp_minimum = at_144(ChangeOfBasisOp());
  SpaceGroup lattice_group = lattice_symmetry_group(minimum.cell, max_delta);
  lattice_group.make_tidy();

  // cctbx's metric_subgroups.derive_result_group_list.
  std::vector<SpaceGroup> subgrs = subgroups_parent_setting(lattice_group);
  std::vector<double> sort_values;
  for (const SpaceGroup &g : subgrs)
    sort_values.push_back(static_cast<double>(g.order_z()) * 1000 +
                          SpaceGroupType(g, false).number());
  std::vector<std::size_t> perm(subgrs.size());
  for (std::size_t i = 0; i < perm.size(); ++i)
    perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
    return sort_values[a] > sort_values[b];
  });

  std::vector<MetricSubgroup> items;
  for (std::size_t i : perm) {
    SpaceGroup &acentric_subgroup = subgrs[i];
    const SpaceGroup acentric_supergroup = metric_supergroup(acentric_subgroup);
    acentric_subgroup.make_tidy();
    const int super_number = SpaceGroupType(acentric_supergroup).number();
    const SpaceGroup super_reference = reference_group_of(super_number);
    const std::string bravais = sgtbx::bravais_symbol(super_reference);
    if (bravais != target_bravais)
      continue;
    SpaceGroup centric_group = acentric_subgroup;
    centric_group.expand_inv(sgtbx::TrVec(sgtbx::sg_t_den));
    MetricSubgroup item;
    item.subsym = CrystalSymmetry::make(minimum.cell, centric_group, false);
    item.supersym =
        CrystalSymmetry::make(minimum.cell, acentric_supergroup, false);
    const ChangeOfBasisOp cb_op_minimum_ref =
        item.subsym.change_of_basis_op_to_reference_setting();
    item.ref_subsym = item.subsym.change_basis(cb_op_minimum_ref);
    if (!is_bravais_holohedry(SpaceGroupType(item.ref_subsym.group).number()))
      continue;
    const ChangeOfBasisOp cb_op_best_cell =
        item.ref_subsym.change_of_basis_op_to_best_cell(true);
    item.best_subsym = item.ref_subsym.change_basis(cb_op_best_cell);
    ChangeOfBasisOp cb_op_inp_best =
        at_144(cb_op_best_cell) * at_144(cb_op_minimum_ref) * cb_op_inp_minimum;
    if (is_similar_to(item.best_subsym.cell, minimum.cell)) {
      const ChangeOfBasisOp cb_op_corr = cb_op_inp_best.inverse();
      try {
        const CrystalSymmetry corr = item.best_subsym.change_basis(cb_op_corr);
        if (corr.group == item.best_subsym.group)
          cb_op_inp_best = cb_op_corr * cb_op_inp_best;
      } catch (const sgtbx::RationalError &e) {
        if (std::string(e.what()).find(
                "Unsuitable value for rational rotation matrix.") ==
            std::string::npos)
          throw;
      }
    }
    item.cb_op_inp_best = cb_op_inp_best;
    SpaceGroup supergroup = item.supersym.group;
    item.max_angular_difference =
        lattice_symmetry_find_max_delta(minimum.cell, supergroup);

    // iotbx_converter's additions.
    item.bravais = bravais;
    item.system = crystal_system_name(super_reference.crystal_system());
    if (item.bravais == "hR" && item.system == "trigonal")
      item.system = "rhombohedral";
    if (item.bravais == "hP" && item.system == "trigonal")
      item.system = "hexagonal";
    item.number = SpaceGroupType(item.subsym.group).number();
    items.push_back(std::move(item));
  }
  cpython_list_sort(items,
                    [](const MetricSubgroup &a, const MetricSubgroup &b) {
                      return bestcmp(a, b) < 0;
                    });
  return items;
}

Mat3 constrain_orientation(const Mat3 &direct, const std::string &system) {
  // s_minimizer: x is flex.double(orient.direct_matrix()); scitbx.lbfgs.run
  // with traditional_convergence_test=0, min_iterations=25, max_calls=1000.
  std::vector<double> x(direct.m.begin(), direct.m.end());
  lbfgs::TerminationParameters termination;
  termination.traditional_convergence_test = false;
  termination.min_iterations = 25;
  termination.max_calls = 1000;
  const lbfgs::Result r = lbfgs::run(
      x,
      [&system](const std::vector<double> &xx, std::vector<double> &g) {
        return s_minimizer_functional(system, xx, g);
      },
      lbfgs::CoreParameters{}, termination);
  if (r.raised)
    throw std::runtime_error(std::string("s_minimizer: line search failed: ") +
                             lbfgs::to_string(r.failure));
  // newOrientation().direct_matrix(): through the reciprocal matrix and back.
  return inverse_scitbx(inverse_scitbx(mat3_of(x)));
}

SymmetryHandler::SymmetryHandler(const std::optional<UnitCell> &unit_cell,
                                 const std::optional<SpaceGroup> &space_group,
                                 double max_delta)
    : max_delta_(max_delta) {
  std::optional<SpaceGroup> target_space_group;
  if (space_group)
    target_space_group = space_group->build_derived_group(false, true);
  SpaceGroup primitive_group;
  if (unit_cell) {
    if (!space_group)
      throw std::invalid_argument(
          "space_group must be provided in combination with unit_cell");
    const CrystalSymmetry inp =
        CrystalSymmetry::make(*unit_cell, *target_space_group);
    inp_cell_ = inp.cell;
    inp_group_ = inp.group;
    if (sgtbx::bravais_symbol(inp.group) == "aP")
      cb_op_inp_ref_ = inp.change_of_basis_op_to_minimum_cell();
    else
      cb_op_inp_ref_ = inp.change_of_basis_op_to_reference_setting();
    const CrystalSymmetry ref = inp.change_basis(cb_op_inp_ref_);
    reference_cell_ = ref.cell;
    reference_group_ = ref.group;
    cb_op_inp_best_ = ref.change_of_basis_op_to_best_cell() * cb_op_inp_ref_;
    const ChangeOfBasisOp ref_to_primitive = ref.group.z2p_op();
    const CrystalSymmetry primitive = ref.change_basis(ref_to_primitive);
    primitive_group = primitive.group;
    primitive_cell_ = primitive.cell;
  } else if (target_space_group) {
    SpaceGroup g = *target_space_group;
    g.make_tidy();
    inp_group_ = g;
    cb_op_inp_ref_ = SpaceGroupType(g).cb_op();
    SpaceGroup ref = g.change_basis(cb_op_inp_ref_);
    ref.make_tidy();
    reference_group_ = ref;
    primitive_group = ref.change_basis(ref.z2p_op());
    primitive_group.make_tidy();
  } else {
    throw std::invalid_argument("SymmetryHandler: no target symmetry");
  }
  cb_op_ref_inp_ = cb_op_inp_ref_.inverse();
  primitive_group.make_tidy();
  target_sg_ref_ = reference_group_of(SpaceGroupType(primitive_group).number());
  target_bravais_ = sgtbx::bravais_symbol(target_sg_ref_);
}

KnownSymmetry SymmetryHandler::known_symmetry(double relative_length_tolerance,
                                              double absolute_angle_tolerance,
                                              double max_delta) const {
  return KnownSymmetry::make(reference_cell_, *reference_group_,
                             relative_length_tolerance,
                             absolute_angle_tolerance, max_delta);
}

std::optional<AppliedSymmetry>
SymmetryHandler::apply_symmetry(const Mat3 &A, const UnitCell &cell) const {
  const std::vector<MetricSubgroup> items =
      iotbx_converter(cell, max_delta_, target_bravais_);
  const MetricSubgroup *best = nullptr;
  double best_angular_difference = 1e8;
  for (const MetricSubgroup &item : items) {
    if (item.bravais != target_bravais_)
      continue;
    if (item.max_angular_difference < best_angular_difference) {
      best_angular_difference = item.max_angular_difference;
      best = &item;
    }
  }
  if (best == nullptr)
    return std::nullopt;

  const ChangeOfBasisOp &cb_op_inp_best = best->cb_op_inp_best;
  const ChangeOfBasisOp cb_op_ref_best =
      best->ref_subsym.change_of_basis_op_to_best_cell();
  const ChangeOfBasisOp cb_op_inp_ref =
      cb_op_ref_best.inverse() * cb_op_inp_best;
  const ChangeOfBasisOp cb_op_ref_inp = cb_op_inp_ref.inverse();

  // crystal_orientation(A, True).change_basis(c^T): the reciprocal matrix in
  // the reference setting.
  const Mat3 astar_ref = A * cb_op_inp_ref.c().r().as_double().transpose();
  const Mat3 direct =
      constrain_orientation(inverse_scitbx(astar_ref), best->system);

  AppliedSymmetry out;
  out.model = CrystalModel::from_real_space(direct.row(0), direct.row(1),
                                            direct.row(2));
  if (!is_similar_to(out.model.cell,
                     target_sg_ref_.average_unit_cell(out.model.cell)))
    throw IncompatibleUnitCell(
        "apply_symmetry: the constrained cell does not fit the target group");
  out.model =
      out.model.change_basis(cb_op_ref_inp.c().r().transpose().as_double());
  out.space_group = target_sg_ref_.change_basis(cb_op_ref_inp);

  if (cb_op_inp_best_) {
    out.cb_op = cb_op_inp_best_->inverse() * cb_op_inp_best;
  } else if (!cb_op_ref_inp_.is_identity_op()) {
    if (*inp_group_ == best->best_subsym.group)
      out.cb_op = cb_op_inp_best;
    else
      out.cb_op = cb_op_ref_inp_ * cb_op_inp_ref;
  } else {
    out.cb_op =
        best->ref_subsym.change_of_basis_op_to_best_cell(false) * cb_op_inp_ref;
  }
  return out;
}

} // namespace mxi::stills

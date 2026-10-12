// Follows cctbx and dials (see known_symmetry.hh). cctbx Copyright (c) 2006 -
// 2026, The Regents of the University of California, through Lawrence
// Berkeley National Laboratory. BSD-3-Clause-LBNL; see
// third_party/cctbx/LICENSE.txt.

#include "stills/known_symmetry.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>

#include "stills/reduction.hh"

namespace mxi::stills {

using sgtbx::ChangeOfBasisOp;
using sgtbx::IntMat3;
using sgtbx::RotMx;
using sgtbx::RotMxInfo;
using sgtbx::RtMx;
using sgtbx::SpaceGroup;
using sgtbx::SpaceGroupType;

namespace {

// sym_mat3::tensor_transpose_transform: c^T t c.
Sym6 tensor_transpose_transform(const Sym6 &t, const Mat3 &cm) {
  const auto &c = cm.m;
  const double ctt[9] = {c[0] * t[0] + c[3] * t[3] + c[6] * t[4],
                         c[3] * t[1] + c[0] * t[3] + c[6] * t[5],
                         c[6] * t[2] + c[0] * t[4] + c[3] * t[5],
                         c[1] * t[0] + c[4] * t[3] + c[7] * t[4],
                         c[4] * t[1] + c[1] * t[3] + c[7] * t[5],
                         c[7] * t[2] + c[1] * t[4] + c[4] * t[5],
                         c[2] * t[0] + c[5] * t[3] + c[8] * t[4],
                         c[5] * t[1] + c[2] * t[3] + c[8] * t[5],
                         c[8] * t[2] + c[2] * t[4] + c[5] * t[5]};
  return {ctt[0] * c[0] + ctt[1] * c[3] + ctt[2] * c[6],
          ctt[3] * c[1] + ctt[4] * c[4] + ctt[5] * c[7],
          ctt[6] * c[2] + ctt[7] * c[5] + ctt[8] * c[8],
          ctt[0] * c[1] + ctt[1] * c[4] + ctt[2] * c[7],
          ctt[0] * c[2] + ctt[1] * c[5] + ctt[2] * c[8],
          ctt[3] * c[2] + ctt[4] * c[5] + ctt[5] * c[8]};
}

IntMat3 imatmul(const IntMat3 &a, const IntMat3 &b) {
  IntMat3 r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r[i * 3 + j] =
          a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
  return r;
}

// Kept out of line so that m_new + diff is rounded to double before the
// subtraction, as cctbx's spoil_optimization ensures.
[[gnu::noinline]] double spoil_optimization(double x) { return x; }

class FastMinimumReduction {
public:
  FastMinimumReduction(const UnitCell &cell, std::size_t iteration_limit,
                       double multiplier, std::size_t min_n_no_change)
      : multiplier_(multiplier), min_n_no_change_(min_n_no_change),
        iteration_limit_(iteration_limit) {
    const Sym6 &g = cell.metric;
    a_ = g[0];
    b_ = g[1];
    c_ = g[2];
    d_ = 2 * g[5];
    e_ = 2 * g[4];
    f_ = 2 * g[3];
    last_abc_ = {-a_, -b_, -c_};
    while (step()) {
    }
  }
  const IntMat3 &r_inv() const { return r_inv_; }

private:
  static int entier(double x) {
    int result = static_cast<int>(x);
    if (x - result < 0)
      result--;
    if (!(x - result < 1))
      result++;
    return result;
  }
  std::array<int, 2> def_test() const {
    int n_zero = 0, n_positive = 0;
    if (0 < d_)
      n_positive += 1;
    else if (!(d_ < 0))
      n_zero += 1;
    if (0 < e_)
      n_positive += 1;
    else if (!(e_ < 0))
      n_zero += 1;
    if (0 < f_)
      n_positive += 1;
    else if (!(f_ < 0))
      n_zero += 1;
    return {n_zero, n_positive};
  }
  bool def_gt_0() const {
    const auto nz_np = def_test();
    return nz_np[1] == 3 || (nz_np[0] == 0 && nz_np[1] == 1);
  }
  bool step() {
    if (b_ < a_)
      n1_action();
    if (c_ < b_) {
      n2_action();
      return true;
    }
    if (def_gt_0()) {
      n3_true_action();
    } else {
      n3_false_action();
      if (!significant_change_test())
        return false;
    }
    if (b2_action())
      return true;
    if (b3_action())
      return true;
    if (b4_action())
      return true;
    if (b5_action())
      return true;
    return false;
  }
  bool significant_change_test(double new_value, int i) const {
    const double m_new = multiplier_ * new_value;
    const double diff = new_value - last_abc_[i];
    const double m_new_plus_diff = spoil_optimization(m_new + diff);
    return m_new_plus_diff - m_new != 0;
  }
  bool significant_change_test() {
    if (significant_change_test(a_, 0) || significant_change_test(b_, 1) ||
        significant_change_test(c_, 2)) {
      n_no_change_ = 0;
    } else {
      n_no_change_++;
      if (n_no_change_ == min_n_no_change_)
        return false;
    }
    last_abc_ = {a_, b_, c_};
    return true;
  }
  void cb_update(const IntMat3 &m) {
    if (n_iterations_ == iteration_limit_)
      throw MinimumReductionLimitExceeded("Iteration limit exceeded.");
    r_inv_ = imatmul(r_inv_, m);
    n_iterations_ += 1;
  }
  void n1_action() {
    cb_update({0, -1, 0, -1, 0, 0, 0, 0, -1});
    std::swap(a_, b_);
    std::swap(d_, e_);
  }
  void n2_action() {
    cb_update({-1, 0, 0, 0, 0, -1, 0, -1, 0});
    std::swap(b_, c_);
    std::swap(e_, f_);
  }
  void n3_true_action() {
    IntMat3 m{1, 0, 0, 0, 1, 0, 0, 0, 1};
    if (d_ < 0)
      m[0] = -1;
    if (e_ < 0)
      m[4] = -1;
    if (f_ < 0)
      m[8] = -1;
    cb_update(m);
    d_ = std::fabs(d_);
    e_ = std::fabs(e_);
    f_ = std::fabs(f_);
  }
  void n3_false_action() {
    IntMat3 m{1, 0, 0, 0, 1, 0, 0, 0, 1};
    int z = -1;
    if (0 < d_)
      m[0] = -1;
    else if (!(d_ < 0))
      z = 0;
    if (0 < e_)
      m[4] = -1;
    else if (!(e_ < 0))
      z = 4;
    if (0 < f_)
      m[8] = -1;
    else if (!(f_ < 0))
      z = 8;
    if (m[0] * m[4] * m[8] < 0) {
      if (z == -1)
        throw CellError("fast_minimum_reduction: internal error");
      m[z] = -1;
    }
    cb_update(m);
    d_ = -std::fabs(d_);
    e_ = -std::fabs(e_);
    f_ = -std::fabs(f_);
  }
  static void degenerate() {
    throw CellError("Degenerate unit cell parameters.");
  }
  bool b2_action() {
    if (!(b_ < std::fabs(d_)))
      return false;
    const int j = entier((d_ + b_) / (2 * b_));
    if (j == 0)
      return false;
    cb_update({1, 0, 0, 0, 1, -j, 0, 0, 1});
    c_ += j * j * b_ - j * d_;
    d_ -= 2 * j * b_;
    e_ -= j * f_;
    if (!(0 < c_))
      degenerate();
    return true;
  }
  bool b3_action() {
    if (!(a_ < std::fabs(e_)))
      return false;
    const int j = entier((e_ + a_) / (2 * a_));
    if (j == 0)
      return false;
    cb_update({1, 0, -j, 0, 1, 0, 0, 0, 1});
    c_ += j * j * a_ - j * e_;
    d_ -= j * f_;
    e_ -= 2 * j * a_;
    if (!(0 < c_))
      degenerate();
    return true;
  }
  bool b4_action() {
    if (!(a_ < std::fabs(f_)))
      return false;
    const int j = entier((f_ + a_) / (2 * a_));
    if (j == 0)
      return false;
    cb_update({1, -j, 0, 0, 1, 0, 0, 0, 1});
    b_ += j * j * a_ - j * f_;
    d_ -= j * e_;
    f_ -= 2 * j * a_;
    if (!(0 < b_))
      degenerate();
    return true;
  }
  bool b5_action() {
    const double de = d_ + e_;
    const double fab = f_ + a_ + b_;
    if (!(de + fab < 0))
      return false;
    const int j = entier((de + fab) / (2 * fab));
    if (j == 0)
      return false;
    cb_update({1, 0, -j, 0, 1, -j, 0, 0, 1});
    c_ += j * j * fab - j * de;
    d_ -= j * (2 * b_ + f_);
    e_ -= j * (2 * a_ + f_);
    if (!(0 < c_))
      degenerate();
    return true;
  }

  double multiplier_;
  std::size_t min_n_no_change_, iteration_limit_;
  IntMat3 r_inv_{1, 0, 0, 0, 1, 0, 0, 0, 1};
  std::size_t n_iterations_ = 0, n_no_change_ = 0;
  double a_, b_, c_, d_, e_, f_;
  std::array<double, 3> last_abc_{};
};

int compare_orthorhombic(const UnitCell &lhs, const UnitCell &rhs) {
  for (int i = 0; i < 3; ++i) {
    if (lhs.params[i] < rhs.params[i])
      return -1;
    if (lhs.params[i] > rhs.params[i])
      return 1;
  }
  return 0;
}

int compare_monoclinic(const UnitCell &lhs, const UnitCell &rhs,
                       int unique_axis, double angular_tolerance) {
  const double lhs_ang = lhs.params[unique_axis + 3];
  const double rhs_ang = rhs.params[unique_axis + 3];
  if (std::fabs(lhs_ang - rhs_ang) < angular_tolerance)
    return compare_orthorhombic(lhs, rhs);
  const double lhs_d90 = std::fabs(lhs_ang - 90);
  const double rhs_d90 = std::fabs(rhs_ang - 90);
  if (std::fabs(lhs_d90 - rhs_d90) > angular_tolerance) {
    if (lhs_d90 < rhs_d90)
      return -1;
    if (lhs_d90 > rhs_d90)
      return 1;
  } else {
    if (lhs_ang > 90 && rhs_ang < 90)
      return -1;
    if (lhs_ang < 90 && rhs_ang > 90)
      return 1;
  }
  if (lhs_ang > rhs_ang)
    return -1;
  if (lhs_ang < rhs_ang)
    return 1;
  return 0;
}

// unit_cell::change_of_basis_op_for_best_monoclinic_beta: "a+c,b,c" if it
// brings beta nearer 90 without going below it.
ChangeOfBasisOp best_monoclinic_beta_op(const UnitCell &cell) {
  const RotMx c_inv_r({12, 0, 0, 0, 12, 0, 12, 0, 12}, 12);
  const ChangeOfBasisOp acbc(RtMx(c_inv_r, sgtbx::cb_t_den).inverse(),
                             RtMx(c_inv_r, sgtbx::cb_t_den));
  const UnitCell alt = change_basis(cell, acbc);
  const double beta = cell.params[4];
  const double beta_alt = alt.params[4];
  if (beta_alt >= 90 && beta_alt < beta)
    return acbc;
  return ChangeOfBasisOp();
}

// Python's math.radians.
double py_radians(double x) { return x * (M_PI / 180.0); }

double cell_volume_factor(double ca, double cb, double cg) {
  const double square = 1 - std::pow(ca, 2.0) - std::pow(cb, 2.0) -
                        std::pow(cg, 2.0) + 2 * ca * cb * cg;
  return square > 0 ? std::sqrt(square) : 0;
}

constexpr double kVolumeWindowSlack = 1e-6;
constexpr double kAngleEps = 1e-9;

// dials' _volume_window: the bounds on V / V_target that is_similar_to
// implies.
std::optional<std::pair<double, double>>
volume_ratio_window(const UnitCell &target, double rel, double ang_tol) {
  if (rel >= 1)
    return std::nullopt;
  const double f_target =
      cell_volume_factor(std::cos(py_radians(target.params[3])),
                         std::cos(py_radians(target.params[4])),
                         std::cos(py_radians(target.params[5])));
  if (f_target <= 0)
    return std::nullopt;
  std::array<std::vector<double>, 3> axes;
  for (int k = 0; k < 3; ++k) {
    const double angle = target.params[3 + k];
    const double lo =
        std::cos(py_radians(std::max(angle - ang_tol, kAngleEps)));
    const double hi =
        std::cos(py_radians(std::min(angle + ang_tol, 180 - kAngleEps)));
    auto &c = axes[k];
    c.push_back(lo);
    if (hi != lo)
      c.push_back(hi);
    if (std::min(lo, hi) <= 0 && 0 <= std::max(lo, hi) &&
        std::find(c.begin(), c.end(), 0.0) == c.end())
      c.push_back(0.0);
  }
  double f_lo = std::numeric_limits<double>::infinity();
  double f_hi = -std::numeric_limits<double>::infinity();
  for (double x : axes[0])
    for (double y : axes[1])
      for (double z : axes[2]) {
        const double f = cell_volume_factor(x, y, z);
        if (f > 0) {
          f_lo = std::min(f_lo, f);
          f_hi = std::max(f_hi, f);
        }
      }
  const double length_lo = std::pow(1 - rel, 3.0);
  return std::make_pair(length_lo * f_lo / f_target * (1 - kVolumeWindowSlack),
                        f_hi / (length_lo * f_target) *
                            (1 + kVolumeWindowSlack));
}

} // namespace

UnitCell change_basis(const UnitCell &cell, const ChangeOfBasisOp &cb) {
  return change_basis(cell, cb.c_inv().r());
}

UnitCell change_basis(const UnitCell &cell, const RotMx &c_inv_r) {
  return UnitCell::from_metric(
      tensor_transpose_transform(cell.metric, c_inv_r.as_double()));
}

bool is_similar_to(const UnitCell &a, const UnitCell &b,
                   double relative_length_tolerance,
                   double absolute_angle_tolerance) {
  for (int i = 0; i < 3; ++i)
    if (std::fabs(std::min(a.params[i], b.params[i]) /
                      std::max(a.params[i], b.params[i]) -
                  1) > relative_length_tolerance)
      return false;
  for (int i = 3; i < 6; ++i)
    if (std::fabs(a.params[i] - b.params[i]) > absolute_angle_tolerance)
      return false;
  return true;
}

IntMat3 fast_minimum_reduction(const UnitCell &cell,
                               std::size_t iteration_limit, double multiplier,
                               std::size_t min_n_no_change) {
  return FastMinimumReduction(cell, iteration_limit, multiplier,
                              min_n_no_change)
      .r_inv();
}

CrystalSymmetry CrystalSymmetry::make(const UnitCell &cell, SpaceGroup group,
                                      bool assert_is_compatible_unit_cell) {
  group.make_tidy();
  const UnitCell average = group.average_unit_cell(cell);
  if (assert_is_compatible_unit_cell && !is_similar_to(cell, average))
    throw IncompatibleUnitCell(
        "Space group is incompatible with unit cell parameters.");
  return CrystalSymmetry{average, std::move(group)};
}

CrystalSymmetry CrystalSymmetry::change_basis(const ChangeOfBasisOp &cb) const {
  return make(stills::change_basis(cell, cb), group.change_basis(cb));
}

ChangeOfBasisOp CrystalSymmetry::change_of_basis_op_to_minimum_cell() const {
  const ChangeOfBasisOp z2p = group.z2p_op();
  const UnitCell p_cell = stills::change_basis(cell, z2p.c_inv().r());
  const IntMat3 r_inv = fast_minimum_reduction(p_cell);
  const ChangeOfBasisOp p2n =
      ChangeOfBasisOp(RtMx(RotMx(r_inv, 1), sgtbx::sg_t_den)).inverse();
  return p2n.new_denominators(z2p) * z2p;
}

ChangeOfBasisOp CrystalSymmetry::change_of_basis_op_to_niggli_cell() const {
  const ChangeOfBasisOp z2p = group.z2p_op();
  const UnitCell p_cell = stills::change_basis(cell, z2p.c_inv().r());
  const NiggliResult red = niggli_reduction(p_cell);
  const ChangeOfBasisOp p2n =
      ChangeOfBasisOp(RtMx(RotMx(red.r_inv, 1), sgtbx::sg_t_den)).inverse();
  return p2n.new_denominators(z2p) * z2p;
}

ChangeOfBasisOp
CrystalSymmetry::change_of_basis_op_to_reference_setting() const {
  return SpaceGroupType(group).cb_op();
}

ChangeOfBasisOp CrystalSymmetry::change_of_basis_op_to_best_cell(
    bool best_monoclinic_beta) const {
  const double angular_tolerance = 3;
  const int number = SpaceGroupType(group).number();
  if (number == 1)
    return change_of_basis_op_to_niggli_cell();
  if (number < 3 || number >= 75)
    return ChangeOfBasisOp();
  if (number == 48 || number == 50 || number == 59 || number == 68 ||
      number == 70)
    throw sgtbx::Error("find_best_cell: the A1983 setting of a space group "
                       "with two origin choices is not ported");
  SpaceGroup standard =
      SpaceGroup::from_hall(sgtbx::hall_symbol_table(number), true);
  standard.make_tidy();
  const ChangeOfBasisOp cb_op_inp_ref = SpaceGroupType(group).cb_op();
  const ChangeOfBasisOp cb_op_std_ref = SpaceGroupType(standard).cb_op();
  const ChangeOfBasisOp cb_op_std_inp = cb_op_inp_ref.inverse() * cb_op_std_ref;
  if (standard.change_basis(cb_op_std_inp) != group)
    throw sgtbx::Error("find_best_cell: standard setting does not map onto "
                       "the input group");
  ChangeOfBasisOp best_cb_op;
  CrystalSymmetry best = *this;
  if (number <= 15) {
    const RotMxInfo two_fold(group(1).r());
    if (std::abs(two_fold.type) != 2)
      throw sgtbx::Error("find_best_cell: no two-fold");
    int n_zero = 0;
    for (int x : two_fold.ev)
      if (x == 0)
        n_zero++;
    if (n_zero != 2)
      throw sgtbx::Error("find_best_cell: two-fold not along an axis");
    const int unique_axis =
        static_cast<int>(std::find(two_fold.ev.begin(), two_fold.ev.end(), 1) -
                         two_fold.ev.begin());
    if (unique_axis == 3)
      throw sgtbx::Error("find_best_cell: two-fold axis not positive");
    for (const RtMx &cb_mx : sgtbx::find_affine(group)) {
      ChangeOfBasisOp cb_op =
          ChangeOfBasisOp(cb_mx).new_denominators(best_cb_op);
      CrystalSymmetry alt = change_basis(cb_op);
      if (alt.group == group) {
        if (best_monoclinic_beta && unique_axis == 1) {
          const ChangeOfBasisOp beta_op = best_monoclinic_beta_op(alt.cell);
          if (!beta_op.is_identity_op()) {
            cb_op.update(beta_op);
            alt = change_basis(cb_op);
          }
        }
        if (compare_monoclinic(best.cell, alt.cell, unique_axis,
                               angular_tolerance) > 0) {
          best_cb_op = cb_op;
          best = alt;
        }
      }
    }
  } else {
    const SpaceGroup affine_group =
        SpaceGroup::from_hall("P 4 3*").change_basis(cb_op_std_inp);
    for (std::size_t i = 0; i < affine_group.order_z(); ++i) {
      const ChangeOfBasisOp cb_op =
          ChangeOfBasisOp(affine_group(i)).new_denominators(best_cb_op);
      const CrystalSymmetry alt = change_basis(cb_op);
      if (alt.group == group) {
        if (compare_orthorhombic(best.cell, alt.cell) > 0) {
          best_cb_op = cb_op;
          best = alt;
        }
      }
    }
  }
  return best_cb_op;
}

double lattice_symmetry_find_max_delta(const UnitCell &reduced_cell,
                                       const SpaceGroup &group) {
  if (group.n_ltr() != 1 || group.f_inv() != 1)
    throw sgtbx::Error("find_max_delta: group must be primitive and acentric");
  const Mat3 &frac = reduced_cell.frac;
  const Mat3 &orth = reduced_cell.orth;
  double min_cos_delta = 1;
  for (std::size_t i = 1; i < group.n_smx(); ++i) {
    const RotMx &r = group.smx(i).r();
    const RotMxInfo r_info = r.info();
    if (r_info.type != 2)
      continue;
    const sgtbx::IntVec3 ev = r_info.ev;
    const sgtbx::IntVec3 evt = r.transpose().info().ev;
    const auto &o = orth.m;
    const auto &f = frac.m;
    const double t[3] = {o[0] * ev[0] + o[1] * ev[1] + o[2] * ev[2],
                         o[3] * ev[0] + o[4] * ev[1] + o[5] * ev[2],
                         o[6] * ev[0] + o[7] * ev[1] + o[8] * ev[2]};
    const double tau[3] = {evt[0] * f[0] + evt[1] * f[3] + evt[2] * f[6],
                           evt[0] * f[1] + evt[1] * f[4] + evt[2] * f[7],
                           evt[0] * f[2] + evt[1] * f[5] + evt[2] * f[8]};
    const double numerator =
        std::fabs(t[0] * tau[0] + t[1] * tau[1] + t[2] * tau[2]);
    const double denominator =
        std::sqrt((t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) *
                  (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]));
    if (denominator == 0)
      throw sgtbx::Error("find_max_delta: zero denominator");
    const double cos_delta = numerator / denominator;
    if (cos_delta < min_cos_delta)
      min_cos_delta = cos_delta;
  }
  if (min_cos_delta > 1.0 - std::numeric_limits<double>::epsilon())
    return 0.0;
  return rad_as_deg(std::acos(min_cos_delta));
}

SpaceGroup lattice_symmetry_group(const UnitCell &reduced_cell,
                                  double max_delta, bool enforce_max_delta) {
  const auto &f = reduced_cell.frac.m;
  const auto &o = reduced_cell.orth.m;
  const double min_cos_delta =
      std::min(std::cos(deg_as_rad(max_delta)),
               1.0 - std::numeric_limits<double>::epsilon());
  const auto &two_folds = sgtbx::reduced_cell_two_folds();
  std::vector<std::size_t> kept;
  std::vector<double> cos_deltas;
  for (std::size_t k = 0; k < two_folds.size(); ++k) {
    const auto &u = two_folds[k].u;
    const auto &h = two_folds[k].h;
    const double t[3] = {o[0] * u[0] + o[1] * u[1] + o[2] * u[2],
                         o[3] * u[0] + o[4] * u[1] + o[5] * u[2],
                         o[6] * u[0] + o[7] * u[1] + o[8] * u[2]};
    const double tau[3] = {h[0] * f[0] + h[1] * f[3] + h[2] * f[6],
                           h[0] * f[1] + h[1] * f[4] + h[2] * f[7],
                           h[0] * f[2] + h[1] * f[5] + h[2] * f[8]};
    const double numerator =
        std::fabs(t[0] * tau[0] + t[1] * tau[1] + t[2] * tau[2]);
    const double denominator =
        std::sqrt((t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) *
                  (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]));
    if (denominator == 0)
      throw sgtbx::Error("lattice_symmetry: zero denominator");
    const double cos_delta = numerator / denominator;
    if (cos_delta >= min_cos_delta) {
      kept.push_back(k);
      cos_deltas.push_back(numerator / denominator);
    }
  }
  if (kept.empty())
    return SpaceGroup();
  std::vector<std::size_t> perm(kept.size());
  for (std::size_t i = 0; i < perm.size(); ++i)
    perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
    return cos_deltas[a] > cos_deltas[b];
  });
  SpaceGroup group;
  for (std::size_t p : perm) {
    const RtMx s(RotMx(two_folds[kept[p]].r, 1), sgtbx::sg_t_den);
    SpaceGroup expanded = group;
    try {
      expanded.expand_smx(s);
    } catch (const sgtbx::NonCrystallographic &) {
      break;
    }
    if (enforce_max_delta &&
        lattice_symmetry_find_max_delta(reduced_cell, expanded) > max_delta)
      continue;
    group = expanded;
  }
  return group;
}

std::vector<SpaceGroup> subgroups_parent_setting(SpaceGroup parent) {
  parent.make_tidy();
  const ChangeOfBasisOp z2p = parent.z2p_op();
  SpaceGroup p_parent = parent.change_basis(z2p);
  p_parent.make_tidy();
  if (p_parent.order_p() != p_parent.order_z())
    throw sgtbx::Error("subgroups: parent is not primitive");
  std::vector<SpaceGroup> p_groups;
  for (std::size_t i = 0; i < p_parent.order_p(); ++i) {
    SpaceGroup group_i;
    group_i.expand_smx(p_parent(i));
    for (std::size_t j = i; j < p_parent.order_p(); ++j) {
      SpaceGroup subgroup = group_i;
      subgroup.expand_smx(p_parent(j));
      subgroup.make_tidy();
      bool seen = false;
      for (const SpaceGroup &g : p_groups)
        if (g == subgroup) {
          seen = true;
          break;
        }
      if (!seen)
        p_groups.push_back(subgroup);
    }
  }
  const ChangeOfBasisOp p2z = z2p.inverse();
  std::vector<SpaceGroup> out;
  out.reserve(p_groups.size());
  for (const SpaceGroup &g : p_groups)
    out.push_back(g.change_basis(p2z));
  return out;
}

SpaceGroup metric_supergroup(SpaceGroup group) {
  group.make_tidy();
  return SpaceGroupType(group)
      .expand_addl_generators_of_euclidean_normalizer(true, true)
      .build_derived_acentric_group();
}

namespace {

// The key dials' groups_cache uses: the target Bravais type and the tidied
// group, whose operations are canonical for the group.
std::string acentric_cache_key(const SpaceGroup &tidy_group,
                               const std::string &target_bravais) {
  std::string key = target_bravais;
  for (const auto &t : tidy_group.ltr().elems())
    for (int i = 0; i < 3; ++i)
      key += "," + std::to_string(t[static_cast<std::size_t>(i)]);
  key += "|" + std::to_string(tidy_group.f_inv());
  for (const RtMx &s : tidy_group.smx()) {
    key += ";";
    for (int i = 0; i < 9; ++i)
      key += std::to_string(s.r()[static_cast<std::size_t>(i)]) + ",";
    for (int i = 0; i < 3; ++i)
      key += std::to_string(s.t()[static_cast<std::size_t>(i)]) + ",";
  }
  return key;
}

AcentricSubgroups compute_acentric_subgroups(SpaceGroup lattice_group,
                                             const std::string &target_bravais);

} // namespace

AcentricSubgroups calc_acentric_subgroups(SpaceGroup lattice_group,
                                          const std::string &target_bravais) {
  // dials memoises this per (target Bravais type, group); the result depends
  // only on the group, so a cache hit returns what a recomputation would.
  static std::mutex mutex;
  static std::map<std::string, AcentricSubgroups> cache;
  lattice_group.make_tidy();
  const std::string key = acentric_cache_key(lattice_group, target_bravais);
  {
    const std::lock_guard<std::mutex> lock(mutex);
    const auto it = cache.find(key);
    if (it != cache.end())
      return it->second;
  }
  AcentricSubgroups result =
      compute_acentric_subgroups(std::move(lattice_group), target_bravais);
  const std::lock_guard<std::mutex> lock(mutex);
  cache.emplace(key, result);
  return result;
}

namespace {

AcentricSubgroups
compute_acentric_subgroups(SpaceGroup lattice_group,
                           const std::string &target_bravais) {
  lattice_group.make_tidy();
  std::vector<SpaceGroup> subgrs = subgroups_parent_setting(lattice_group);
  std::vector<double> sort_values;
  sort_values.reserve(subgrs.size());
  for (const SpaceGroup &g : subgrs) {
    const int number = SpaceGroupType(g, false).number();
    sort_values.push_back(static_cast<double>(g.order_z()) * 1000 + number);
  }
  std::vector<std::size_t> perm(subgrs.size());
  for (std::size_t i = 0; i < perm.size(); ++i)
    perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
    return sort_values[a] > sort_values[b];
  });
  AcentricSubgroups out;
  for (std::size_t i : perm) {
    SpaceGroup &subgroup = subgrs[i];
    subgroup.make_tidy();
    const ChangeOfBasisOp cb_op = SpaceGroupType(subgroup).cb_op();
    const SpaceGroup ref = subgroup.change_basis(cb_op);
    if (sgtbx::bravais_symbol(ref) != target_bravais)
      continue;
    out.subgroups.push_back(subgroup);
    out.cb_ops.push_back(cb_op);
    out.supergroups.push_back(metric_supergroup(subgroup));
  }
  return out;
}

} // namespace

std::optional<MatchingSymmetry>
find_matching_symmetry(const UnitCell &unit_cell, double max_delta,
                       const std::string &target_bravais,
                       bool best_monoclinic_beta) {
  const CrystalSymmetry input = CrystalSymmetry::make(unit_cell, SpaceGroup());
  const ChangeOfBasisOp cb_op_inp_minimum =
      input.change_of_basis_op_to_minimum_cell();
  const CrystalSymmetry minimum = input.change_basis(cb_op_inp_minimum);
  const SpaceGroup lattice_group =
      lattice_symmetry_group(minimum.cell, max_delta, true);
  const AcentricSubgroups acentric =
      calc_acentric_subgroups(lattice_group, target_bravais);

  std::optional<MatchingSymmetry> best;
  double best_angular_difference = 1e8;
  for (std::size_t k = 0; k < acentric.subgroups.size(); ++k) {
    const CrystalSymmetry subsym =
        CrystalSymmetry::make(minimum.cell, acentric.subgroups[k], false);
    const CrystalSymmetry ref_subsym = subsym.change_basis(acentric.cb_ops[k]);
    const ChangeOfBasisOp cb_op_best_cell =
        ref_subsym.change_of_basis_op_to_best_cell(best_monoclinic_beta);
    const CrystalSymmetry best_subsym =
        ref_subsym.change_basis(cb_op_best_cell);

    // dials composes cb_op_inp_best at r_den = t_den = 144 and, when the best
    // cell is similar to the input, tries its inverse; that only changes the
    // operator, but its change of basis can raise, so it is kept.
    const ChangeOfBasisOp cb_op_inp_best =
        ChangeOfBasisOp(cb_op_best_cell.c().new_denominators(144, 144)) *
        ChangeOfBasisOp(acentric.cb_ops[k].c().new_denominators(144, 144)) *
        ChangeOfBasisOp(cb_op_inp_minimum.c().new_denominators(144, 144));
    if (is_similar_to(best_subsym.cell, input.cell)) {
      try {
        best_subsym.change_basis(cb_op_inp_best.inverse());
      } catch (const sgtbx::RationalError &e) {
        if (std::string(e.what()).find(
                "Unsuitable value for rational rotation matrix.") ==
            std::string::npos)
          throw;
      }
    }

    const double max_angular_difference =
        lattice_symmetry_find_max_delta(minimum.cell, acentric.supergroups[k]);
    if (max_angular_difference < best_angular_difference) {
      best_angular_difference = max_angular_difference;
      best = MatchingSymmetry{max_angular_difference, subsym, ref_subsym,
                              best_subsym};
    }
  }
  return best;
}

KnownSymmetry KnownSymmetry::make(const std::optional<UnitCell> &cell,
                                  const SpaceGroup &group,
                                  double relative_length_tolerance,
                                  double absolute_angle_tolerance,
                                  double max_delta) {
  KnownSymmetry k;
  k.relative_length_tolerance = relative_length_tolerance;
  k.absolute_angle_tolerance = absolute_angle_tolerance;
  k.max_delta = max_delta;
  SpaceGroup g = group;
  g.make_tidy();
  const ChangeOfBasisOp cb_op_ref_to_primitive = g.z2p_op();
  SpaceGroup primitive_group;
  if (cell) {
    const CrystalSymmetry target{*cell, g};
    const CrystalSymmetry primitive =
        target.change_basis(cb_op_ref_to_primitive);
    const CrystalSymmetry reference =
        target.change_basis(target.change_of_basis_op_to_reference_setting());
    const CrystalSymmetry best =
        reference.change_basis(reference.change_of_basis_op_to_best_cell());
    k.target_unit_cell = best.cell;
    const auto ratio = volume_ratio_window(best.cell, relative_length_tolerance,
                                           absolute_angle_tolerance);
    if (ratio) {
      const double v = primitive.cell.volume;
      k.volume_window = std::make_pair(ratio->first * v, ratio->second * v);
    }
    primitive_group = primitive.group;
  } else {
    primitive_group = g.change_basis(cb_op_ref_to_primitive);
    primitive_group.make_tidy();
  }
  const int number = SpaceGroupType(primitive_group).number();
  SpaceGroup reference =
      SpaceGroup::from_hall(sgtbx::hall_symbol_table(number), true);
  reference.make_tidy();
  k.target_bravais = sgtbx::bravais_symbol(reference);
  return k;
}

KnownSymmetryResult filter_known_symmetry(const UnitCell &cell,
                                          const KnownSymmetry &target) {
  if (target.volume_window && !(target.volume_window->first <= cell.volume &&
                                cell.volume <= target.volume_window->second))
    return {KnownSymmetryOutcome::outside_volume_window, std::nullopt};
  std::optional<MatchingSymmetry> best =
      find_matching_symmetry(cell, target.max_delta, target.target_bravais);
  if (!best)
    return {KnownSymmetryOutcome::no_matching_subgroup, std::nullopt};
  if (target.target_unit_cell &&
      !is_similar_to(best->best_subsym.cell, *target.target_unit_cell,
                     target.relative_length_tolerance,
                     target.absolute_angle_tolerance))
    return {KnownSymmetryOutcome::not_similar_to_target, std::move(best)};
  return {KnownSymmetryOutcome::accepted, std::move(best)};
}

} // namespace mxi::stills

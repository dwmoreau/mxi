// Space-group arithmetic as cctbx's sgtbx does it: rational symmetry
// operations, groups built by cctbx's expansion and kept in its order, the
// Hall symbol parser, and the space-group type with its change of basis to the
// reference setting (Grosse-Kunstleve 1999, made unique by cctbx's search over
// the affine normalizer). The floating results downstream depend on the order
// in which a group holds its operations, so expansion, tidying and change of
// basis follow cctbx statement for statement.
//
// Follows cctbx_project cctbx/sgtbx: rot_mx, tr_vec, rt_mx, tr_group,
// change_of_basis_op, space_group, hall_in, lattice_tr, rot_mx_info,
// row_echelon_solve, select_generators, space_group_type, find_affine and
// reference_settings/normalizer; scitbx/matrix/row_echelon.h.
// cctbx Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "stills/cell.hh"

namespace mxi::stills::sgtbx {

// Denominators cctbx uses by default: symmetry translations in twelfths,
// change-of-basis rotations in twelfths and translations in 144ths.
inline constexpr int sg_t_den = 12;
inline constexpr int cb_r_den = 12;
inline constexpr int cb_t_den = 144;
// The most representative operations a space group holds (m-3m: 48 / 2).
inline constexpr std::size_t n_max_repr_rot_mx = 24;

struct Error : std::runtime_error {
  using std::runtime_error::runtime_error;
};
// A rational matrix or vector that cannot take the requested denominator.
struct RationalError : Error {
  using Error::Error;
};
struct NonCrystallographic : Error {
  NonCrystallographic()
      : Error("Non-crystallographic rotation matrix encountered.") {}
};

using IntMat3 = std::array<int, 9>;
using IntVec3 = std::array<int, 3>;

class TrVec {
public:
  explicit TrVec(int den = sg_t_den) : num_{0, 0, 0}, den_(den) {}
  TrVec(const IntVec3 &v, int den) : num_(v), den_(den) {}
  TrVec(int v0, int v1, int v2, int den) : num_{v0, v1, v2}, den_(den) {}

  const IntVec3 &num() const { return num_; }
  IntVec3 &num() { return num_; }
  int operator[](std::size_t i) const { return num_[i]; }
  int &operator[](std::size_t i) { return num_[i]; }
  int den() const { return den_; }
  int &den() { return den_; }

  bool operator==(const TrVec &o) const {
    return den_ == o.den_ && num_ == o.num_;
  }
  bool operator!=(const TrVec &o) const { return !(*this == o); }
  bool is_valid() const { return den_ != 0; }
  bool is_zero() const { return num_ == IntVec3{0, 0, 0}; }

  TrVec new_denominator(int new_den) const;
  TrVec scale(int factor) const;
  TrVec mod_positive() const;
  TrVec mod_short() const;
  TrVec cancel() const;
  TrVec operator-() const {
    return TrVec({-num_[0], -num_[1], -num_[2]}, den_);
  }
  TrVec operator+(const TrVec &o) const;
  TrVec operator-(const TrVec &o) const;
  TrVec operator*(int f) const;
  TrVec operator/(int d) const;

private:
  IntVec3 num_;
  int den_;
};

inline TrVec operator*(int f, const TrVec &t) { return t * f; }

struct RotMxInfo;

class RotMx {
public:
  explicit RotMx(int den = 1, int diagonal = 1);
  RotMx(const IntMat3 &m, int den) : num_(m), den_(den) {}

  const IntMat3 &num() const { return num_; }
  IntMat3 &num() { return num_; }
  int operator[](std::size_t i) const { return num_[i]; }
  int &operator[](std::size_t i) { return num_[i]; }
  int den() const { return den_; }

  bool operator==(const RotMx &o) const {
    return den_ == o.den_ && num_ == o.num_;
  }
  bool operator!=(const RotMx &o) const { return !(*this == o); }
  bool is_valid() const { return den_ != 0; }
  bool is_unit_mx() const;
  RotMx minus_unit_mx() const;
  RotMx new_denominator(int new_den) const;
  RotMx scale(int factor) const;
  RotMx inverse() const;
  RotMx transpose() const;
  RotMx operator-() const;
  RotMx operator+(const RotMx &o) const;
  RotMx operator-(const RotMx &o) const;
  RotMx operator*(const RotMx &o) const;
  TrVec operator*(const TrVec &t) const;
  IntVec3 operator*(const IntVec3 &v) const;
  RotMx operator*(int f) const;
  RotMx operator/(int d) const;
  RotMx cancel() const;

  // Rotation-part type: 1, 2, 3, 4, 6, -1, -2, -3, -4, -6, or 0.
  int type() const;
  int order(int type = 0) const;
  // I + R + R^2 + ... + R^(order-1); the denominator must be 1.
  RotMx accumulate(int type = 0) const;
  RotMxInfo info() const;
  // The elements divided by the denominator.
  Mat3 as_double() const;

private:
  IntMat3 num_;
  int den_;
};

// The type, axis direction and sense of a rotation part (den 1).
struct RotMxInfo {
  int type = 0;
  IntVec3 ev{0, 0, 0};
  int sense = 0;
  explicit RotMxInfo(const RotMx &r);
  RotMxInfo() = default;
};

class RtMx {
public:
  explicit RtMx(int r_den = 1, int t_den = sg_t_den) : r_(r_den), t_(t_den) {}
  RtMx(const RotMx &r, const TrVec &t) : r_(r), t_(t) {}
  RtMx(const RotMx &r, int t_den) : r_(r), t_(t_den) {}
  RtMx(const TrVec &t, int r_den) : r_(r_den), t_(t) {}

  const RotMx &r() const { return r_; }
  RotMx &r() { return r_; }
  const TrVec &t() const { return t_; }
  TrVec &t() { return t_; }

  bool is_valid() const { return r_.is_valid() && t_.is_valid(); }
  RtMx unit_mx() const { return RtMx(r_.den(), t_.den()); }
  bool is_unit_mx() const { return r_.is_unit_mx() && t_.is_zero(); }
  RtMx new_denominators(int r_den, int t_den = 0) const;
  RtMx new_denominators(const RtMx &other) const {
    return new_denominators(other.r().den(), other.t().den());
  }
  RtMx scale(int factor_r, int factor_t = 0) const;
  RtMx mod_positive() const { return RtMx(r_, t_.mod_positive()); }
  RtMx mod_short() const { return RtMx(r_, t_.mod_short()); }
  RtMx pre_multiply_inv_t(const TrVec &inv_t) const {
    return RtMx(-r_, -t_ + inv_t);
  }
  RtMx inverse() const;
  RtMx operator-() const { return RtMx(-r_, -t_); }
  RtMx operator*(const RtMx &o) const;
  RtMx operator+(const TrVec &t) const { return RtMx(r_, t_ + t); }
  TrVec t_intrinsic_part() const;
  std::string as_xyz() const;

  // cctbx's ordering of symmetry operations, the one make_tidy sorts by.
  bool operator<(const RtMx &o) const;
  bool operator==(const RtMx &o) const { return !(*this < o) && !(o < *this); }
  bool operator!=(const RtMx &o) const { return !(*this == o); }

private:
  RotMx r_;
  TrVec t_;
};

class ChangeOfBasisOp {
public:
  explicit ChangeOfBasisOp(int r_den = cb_r_den, int t_den = cb_t_den)
      : c_(r_den, t_den), c_inv_(r_den, t_den) {}
  explicit ChangeOfBasisOp(const RtMx &c) : c_(c), c_inv_(c.inverse()) {}
  ChangeOfBasisOp(const RtMx &c, const RtMx &c_inv) : c_(c), c_inv_(c_inv) {}

  const RtMx &c() const { return c_; }
  const RtMx &c_inv() const { return c_inv_; }
  const RtMx &select(bool inv) const { return inv ? c_inv_ : c_; }
  bool is_valid() const { return c_.is_valid() && c_inv_.is_valid(); }
  ChangeOfBasisOp identity_op() const {
    return ChangeOfBasisOp(c_.unit_mx(), c_inv_.unit_mx());
  }
  bool is_identity_op() const { return c_.is_unit_mx() && c_inv_.is_unit_mx(); }
  ChangeOfBasisOp new_denominators(int r_den, int t_den) const {
    return ChangeOfBasisOp(c_.new_denominators(r_den, t_den),
                           c_inv_.new_denominators(r_den, t_den));
  }
  ChangeOfBasisOp new_denominators(const ChangeOfBasisOp &other) const {
    return ChangeOfBasisOp(c_.new_denominators(other.c()),
                           c_inv_.new_denominators(other.c_inv()));
  }
  ChangeOfBasisOp inverse() const { return ChangeOfBasisOp(c_inv_, c_); }
  ChangeOfBasisOp mod_short() const {
    return ChangeOfBasisOp(c_.mod_short(), c_inv_.mod_short());
  }
  // cctbx's change_of_basis_op(symbol, "", r_den, t_den): "x+1/2,y,z",
  // "h,k,l" or "a+c,b,c".
  static ChangeOfBasisOp from_string(const std::string &symbol,
                                     int r_den = cb_r_den,
                                     int t_den = cb_t_den);
  RotMx operator()(const RotMx &r) const;
  RtMx operator()(const RtMx &s) const;
  TrVec operator()(const TrVec &t, int sign_identity) const;
  void update(const ChangeOfBasisOp &other);
  void update(const TrVec &shift);
  ChangeOfBasisOp operator*(const ChangeOfBasisOp &rhs) const;
  std::string as_xyz() const { return c_.as_xyz(); }

private:
  RtMx c_, c_inv_;
};

// The lattice translations of a space group, the zero vector first.
class TrGroup {
public:
  explicit TrGroup(int t_den = sg_t_den) { elems_.push_back(TrVec(t_den)); }
  std::size_t size() const { return elems_.size(); }
  const TrVec &operator[](std::size_t i) const { return elems_[i]; }
  const std::vector<TrVec> &elems() const { return elems_; }
  int t_den() const { return elems_[0].den(); }
  bool operator==(const TrGroup &o) const { return elems_ == o.elems_; }
  bool contains(const TrVec &t) const;
  bool add(const TrVec &t);
  bool expand(const TrVec &t);
  TrGroup change_basis(const ChangeOfBasisOp &cb_op) const;
  // P, A, B, C, I, R, Q, H or F; '\0' for anything else.
  char conventional_centring_type_symbol() const;
  ChangeOfBasisOp conventional_z2p_op(int r_den = cb_r_den,
                                      int t_den = cb_t_den) const;
  TrVec tidy(const TrVec &t) const;
  void sort_tail();

private:
  std::vector<TrVec> elems_;
};

enum class CrystalSystem {
  undefined,
  triclinic,
  monoclinic,
  orthorhombic,
  tetragonal,
  trigonal,
  hexagonal,
  cubic
};

// cctbx's matrix_group::code.
struct MatrixGroup {
  CrystalSystem x;
  int l, p, m;
  bool operator==(const MatrixGroup &o) const {
    return x == o.x && l == o.l && p == o.p && m == o.m;
  }
  bool operator!=(const MatrixGroup &o) const { return !(*this == o); }
  MatrixGroup point_group_type() const { return {x, l - p, 0, m + p}; }
  MatrixGroup laue_group_type() const { return {x, 0, 0, m + l}; }
  CrystalSystem crystal_system() const { return x; }
  const char *label() const;
};

namespace matrix_group {
// clang-format off
inline constexpr MatrixGroup unknown   {CrystalSystem::undefined,    0, 0,  1};
inline constexpr MatrixGroup code_1    {CrystalSystem::triclinic,    1, 0,  2};
inline constexpr MatrixGroup code_1b   {CrystalSystem::triclinic,    0, 0,  3};
inline constexpr MatrixGroup code_2    {CrystalSystem::monoclinic,   2, 0,  4};
inline constexpr MatrixGroup code_m    {CrystalSystem::monoclinic,   1, 0,  5};
inline constexpr MatrixGroup code_2_m  {CrystalSystem::monoclinic,   0, 0,  6};
inline constexpr MatrixGroup code_222  {CrystalSystem::orthorhombic, 2, 0,  7};
inline constexpr MatrixGroup code_mm2  {CrystalSystem::orthorhombic, 1, 0,  8};
inline constexpr MatrixGroup code_mmm  {CrystalSystem::orthorhombic, 0, 0,  9};
inline constexpr MatrixGroup code_4    {CrystalSystem::tetragonal,   2, 0, 10};
inline constexpr MatrixGroup code_4b   {CrystalSystem::tetragonal,   1, 0, 11};
inline constexpr MatrixGroup code_4_m  {CrystalSystem::tetragonal,   0, 0, 12};
inline constexpr MatrixGroup code_422  {CrystalSystem::tetragonal,   4, 0, 13};
inline constexpr MatrixGroup code_4mm  {CrystalSystem::tetragonal,   3, 0, 14};
inline constexpr MatrixGroup code_4b2m {CrystalSystem::tetragonal,   2, 1, 15};
inline constexpr MatrixGroup code_4bm2 {CrystalSystem::tetragonal,   1, 0, 16};
inline constexpr MatrixGroup code_4_mmm{CrystalSystem::tetragonal,   0, 0, 17};
inline constexpr MatrixGroup code_3    {CrystalSystem::trigonal,     1, 0, 18};
inline constexpr MatrixGroup code_3b   {CrystalSystem::trigonal,     0, 0, 19};
inline constexpr MatrixGroup code_321  {CrystalSystem::trigonal,     8, 2, 20};
inline constexpr MatrixGroup code_312  {CrystalSystem::trigonal,     7, 1, 21};
inline constexpr MatrixGroup code_32   {CrystalSystem::trigonal,     6, 0, 22};
inline constexpr MatrixGroup code_3m1  {CrystalSystem::trigonal,     5, 2, 23};
inline constexpr MatrixGroup code_31m  {CrystalSystem::trigonal,     4, 1, 24};
inline constexpr MatrixGroup code_3m   {CrystalSystem::trigonal,     3, 0, 25};
inline constexpr MatrixGroup code_3bm1 {CrystalSystem::trigonal,     2, 2, 26};
inline constexpr MatrixGroup code_3b1m {CrystalSystem::trigonal,     1, 1, 27};
inline constexpr MatrixGroup code_3bm  {CrystalSystem::trigonal,     0, 0, 28};
inline constexpr MatrixGroup code_6    {CrystalSystem::hexagonal,    2, 0, 29};
inline constexpr MatrixGroup code_6b   {CrystalSystem::hexagonal,    1, 0, 30};
inline constexpr MatrixGroup code_6_m  {CrystalSystem::hexagonal,    0, 0, 31};
inline constexpr MatrixGroup code_622  {CrystalSystem::hexagonal,    4, 0, 32};
inline constexpr MatrixGroup code_6mm  {CrystalSystem::hexagonal,    3, 0, 33};
inline constexpr MatrixGroup code_6b2m {CrystalSystem::hexagonal,    2, 1, 34};
inline constexpr MatrixGroup code_6bm2 {CrystalSystem::hexagonal,    1, 0, 35};
inline constexpr MatrixGroup code_6_mmm{CrystalSystem::hexagonal,    0, 0, 36};
inline constexpr MatrixGroup code_23   {CrystalSystem::cubic,        1, 0, 37};
inline constexpr MatrixGroup code_m3b  {CrystalSystem::cubic,        0, 0, 38};
inline constexpr MatrixGroup code_432  {CrystalSystem::cubic,        2, 0, 39};
inline constexpr MatrixGroup code_4b3m {CrystalSystem::cubic,        1, 0, 40};
inline constexpr MatrixGroup code_m3bm {CrystalSystem::cubic,        0, 0, 41};
// clang-format on
} // namespace matrix_group

class SpaceGroup {
public:
  explicit SpaceGroup(bool no_expand = false, int t_den = sg_t_den);
  // A group from a Hall symbol (cctbx's space_group(hall_symbol, pedantic,
  // no_centring_type_symbol, no_expand, t_den)).
  static SpaceGroup from_hall(const std::string &hall, bool pedantic = false,
                              bool no_centring_type_symbol = false,
                              bool no_expand = false, int t_den = sg_t_den);

  SpaceGroup &expand_ltr(const TrVec &t);
  SpaceGroup &expand_inv(const TrVec &inv_t);
  SpaceGroup &expand_smx(const RtMx &s);
  std::size_t expand_conventional_centring_type(char symbol);
  std::size_t parse_hall_symbol(const std::string &hall, bool pedantic,
                                bool no_centring_type_symbol);

  int r_den() const { return smx_[0].r().den(); }
  int t_den() const { return smx_[0].t().den(); }
  std::size_t n_smx() const { return smx_.size(); }
  std::size_t n_ltr() const { return ltr_.size(); }
  std::size_t f_inv() const { return f_inv_; }
  std::size_t order_p() const { return f_inv_ * n_smx(); }
  std::size_t order_z() const { return ltr_.size() * f_inv_ * n_smx(); }
  bool is_centric() const { return f_inv_ == 2; }
  const TrGroup &ltr() const { return ltr_; }
  const std::vector<RtMx> &smx() const { return smx_; }
  const RtMx &smx(std::size_t i) const { return smx_[i]; }
  TrVec inv_t(bool tidy = false) const;
  RtMx operator()(std::size_t i_ltr, std::size_t i_inv,
                  std::size_t i_smx) const;
  RtMx operator()(std::size_t i_op) const;
  char conventional_centring_type_symbol() const {
    return ltr_.conventional_centring_type_symbol();
  }

  SpaceGroup &make_tidy();
  bool operator==(const SpaceGroup &o) const;
  bool operator!=(const SpaceGroup &o) const { return !(*this == o); }
  SpaceGroup change_basis(const ChangeOfBasisOp &cb_op) const;
  SpaceGroup build_derived_group(bool discard_z, bool add_inv) const;
  SpaceGroup build_derived_acentric_group() const;
  ChangeOfBasisOp z2p_op(int r_den = cb_r_den, int t_den = cb_t_den) const;
  MatrixGroup point_group_type() const;
  CrystalSystem crystal_system() const {
    return point_group_type().crystal_system();
  }

  // The metric averaged over the rotation parts of the representative
  // operations, sum(R^T G R) / n_smx, as cctbx's average_metrical_matrix.
  Sym6 average_metrical_matrix(const Sym6 &g) const;
  UnitCell average_unit_cell(const UnitCell &cell) const;

private:
  void reset(int t_den);
  void add_inv(const TrVec &t);
  void add_smx(const RtMx &s);
  ChangeOfBasisOp construct_z2p_op(int r_den, int t_den) const;

  bool no_expand_;
  std::size_t n_lsl_ = 1, n_ssl_ = 1, f_inv_ = 1;
  TrGroup ltr_;
  TrVec inv_t_;
  std::vector<RtMx> smx_;
  bool is_tidy_ = false;
};

class SpaceGroupType {
public:
  // cctbx's space_group_type(group, tidy_cb_op, r_den, t_den).
  explicit SpaceGroupType(const SpaceGroup &group, bool tidy_cb_op = true,
                          int r_den = cb_r_den, int t_den = cb_t_den);
  int number() const { return number_; }
  // The change of basis from the group's setting to the reference setting.
  const ChangeOfBasisOp &cb_op() const { return cb_op_; }
  const SpaceGroup &group() const { return group_; }
  std::vector<RtMx> addl_generators_of_euclidean_normalizer(bool k2l,
                                                            bool l2n) const;
  SpaceGroup expand_addl_generators_of_euclidean_normalizer(bool k2l,
                                                            bool l2n) const;

private:
  SpaceGroup group_;
  int number_ = 0;
  ChangeOfBasisOp cb_op_;
};

// The Bravais symbol cctbx's bravais_types.bravais_lattice gives a group, in
// the group's own setting: crystal family letter and centring, e.g. "mP".
std::string bravais_symbol(const SpaceGroup &g);

// cctbx's find_affine(group, range=2): the integer change-of-basis matrices
// with elements in [-range, range] and determinant 1 that map the group onto
// itself, in cctbx's order. Only the general (n_smx > 1) algorithm is here.
std::vector<RtMx> find_affine(const SpaceGroup &group, int range = 2);

// cctbx's reference tables (sgtbx_tables.cc).
const char *hall_symbol_table(int number);
const MatrixGroup &matrix_group_code_table(int number);
struct NormalizerGenerators {
  const char *k2l;
  const char *l2n;
};
const NormalizerGenerators &normalizer_table(int number);
// A two-fold of a reduced cell: its matrix, and its axis in direct and in
// reciprocal space (cctbx's lattice_symmetry::reduced_cell_two_folds).
struct TwoFold {
  IntMat3 r;
  IntVec3 u, h;
};
const std::vector<TwoFold> &reduced_cell_two_folds();
std::vector<RtMx> normalizer_addl_generators(int sg_number, bool affine,
                                             bool k2l, bool l2n);

} // namespace mxi::stills::sgtbx

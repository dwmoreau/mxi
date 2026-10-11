// A unit cell and a crystal model computed as cctbx's uctbx::unit_cell and
// dxtbx's Crystal compute them, operation for operation, so that candidate
// orientations built here equal dials' to the last bit.

#pragma once

#include <array>
#include <stdexcept>

#include "linalg.hh"

namespace mxi::stills {

// A cell that cctbx refuses: a non-positive parameter, an angle of 180 degrees
// or more, a non-positive volume, or a metric that is not positive definite.
struct CellError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

// The six elements of a symmetric 3x3 matrix in cctbx's sym_mat3 order:
// (00, 11, 22, 01, 02, 12).
using Sym6 = std::array<double, 6>;

struct UnitCell {
  // a, b, c in A; alpha, beta, gamma in degrees.
  std::array<double, 6> params{};
  std::array<double, 3> cos_ang{}, sin_ang{};
  double volume = 0.0;
  std::array<double, 6> r_params{};
  std::array<double, 3> r_cos_ang{}, r_sin_ang{};
  // The real-space metric G, in Sym6 order.
  Sym6 metric{};
  // Fractional to Cartesian and back, PDB convention (x along a).
  Mat3 orth, frac;

  static UnitCell from_parameters(const std::array<double, 6> &params);
  static UnitCell from_metric(const Sym6 &g);
  // From an orthogonalisation matrix O, through its metric O^T O.
  static UnitCell from_orthogonalisation(const Mat3 &o);
};

// cctbx's degree conversions: multiplication and division by atan(1)/45.
double deg_as_rad(double deg);
double rad_as_deg(double rad);

// scitbx's mat3::inverse: the transposed cofactors, each divided by the
// determinant. Throws CellError on a zero determinant.
Mat3 inverse_scitbx(const Mat3 &m);

// A crystal as dxtbx holds one: U, B (the transpose of the cell's
// fractionalisation matrix) and the cell. A = U B.
struct CrystalModel {
  Mat3 U, B;
  UnitCell cell;

  Mat3 A() const { return U * B; }

  // dxtbx's Crystal(real_space_a, real_space_b, real_space_c).
  static CrystalModel from_real_space(const Vec3 &a, const Vec3 &b,
                                      const Vec3 &c);
  // dxtbx's Crystal::set_A, which sets B through crystal_orientation.
  static CrystalModel from_A(const Mat3 &A);

  // dxtbx's change_basis: A' = A M, where M is the transpose of the change of
  // basis's rotation part, here passed already transposed.
  CrystalModel change_basis(const Mat3 &M) const;
};

} // namespace mxi::stills

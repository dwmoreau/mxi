// The crystal model of a still as dials refines it: an orientation U, moved by
// three small rotations about the laboratory axes, and a reciprocal-space
// orthogonalisation matrix B, set by the independent elements of the reciprocal
// metric G* = B^T B that the space group allows (Waterman et al. 2016, Table 1
// and eqs. 12-13).

#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "linalg.hh"

namespace mxi::stills {

// A crystal model that cannot be formed from its parameters: a metric that is
// not positive definite, or a cell too small to be real.
struct ComposeError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// The six elements of G*, in the order (a*.a*, b*.b*, c*.c*, a*.b*, a*.c*,
// b*.c*).
using Metric = std::array<double, 6>;

Metric metric_of(const Mat3 &B);

// The constraint a point group puts on G*: G* = sum_k g_k basis[k], with g_k
// the element of G* in column[k]. The independent elements are the free
// columns of the reduced row-echelon form of the equations R G* R^T = G* over
// the group's rotations, taken in column order.
struct CellConstraint {
  std::vector<int> column;
  std::vector<Metric> basis;

  std::size_t size() const { return column.size(); }

  // Rotations as integer matrices on fractional coordinates, row major,
  // x' = R x.
  static CellConstraint
  from_rotations(const std::vector<std::array<int, 9>> &rotations);
  static CellConstraint from_hall(const std::string &hall);

  std::vector<double> independent(const Metric &g) const;
  Metric metric(const std::vector<double> &independent) const;
};

// The lower-triangular B with B^T B = G*, which is dxtbx's B (the transpose of
// the direct cell's fractionalisation matrix). Throws ComposeError where the
// metric has no such factor or the direct cell's volume is at most 40 A^3.
Mat3 b_from_metric(const Metric &g);

// dB for a change dG* in the metric, at B: with M = B^-T dG* B^-1, dB B^-1 is
// the lower triangle of M with its diagonal halved.
Mat3 db_from_metric(const Mat3 &B, const Mat3 &B_inverse, const Metric &dg);

// U = Phi3 Phi2 Phi1 U0, the angles in mrad about the laboratory x, y and z
// axes, and the derivatives of U with respect to each angle (per mrad).
struct Orientation {
  Mat3 U;
  std::array<Mat3, 3> dU;
};
Orientation compose_orientation(const Mat3 &U0,
                                const std::array<double, 3> &mrad);

} // namespace mxi::stills

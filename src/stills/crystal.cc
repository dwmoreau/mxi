// Derived from dials and dxtbx.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// dxtbx: Copyright (c) 2019 Diamond Light Source, Lawrence Berkeley National
// Laboratory and the Science and Technology Facilities Council. All rights
// reserved. BSD-3-Clause; see third_party/dxtbx/LICENSE.txt.

#include "crystal.hh"

#include <cmath>

#include "symmetry.hh"

namespace mxi::stills {

namespace {

// Positions of the six metric elements in a symmetric 3x3 matrix.
constexpr int kRow[6] = {0, 1, 2, 0, 0, 1};
constexpr int kCol[6] = {0, 1, 2, 1, 2, 2};

Mat3 symmetric(const Metric &g) {
  return {g[0], g[3], g[4], g[3], g[1], g[5], g[4], g[5], g[2]};
}

Metric packed(const Mat3 &m) {
  return {m(0, 0), m(1, 1), m(2, 2), m(0, 1), m(0, 2), m(1, 2)};
}

// Rotation by angle (radians) about a unit axis.
Mat3 axis_angle(const Vec3 &u, double angle) {
  const double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
  return {
      t * u.x * u.x + c,       t * u.x * u.y - s * u.z, t * u.x * u.z + s * u.y,
      t * u.x * u.y + s * u.z, t * u.y * u.y + c,       t * u.y * u.z - s * u.x,
      t * u.x * u.z - s * u.y, t * u.y * u.z + s * u.x, t * u.z * u.z + c};
}

// dR/d(angle) for the rotation above.
Mat3 axis_angle_derivative(const Vec3 &u, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {
      s * u.x * u.x - s,       s * u.x * u.y - c * u.z, s * u.x * u.z + c * u.y,
      s * u.y * u.x + c * u.z, s * u.y * u.y - s,       s * u.y * u.z - c * u.x,
      s * u.z * u.x - c * u.y, s * u.z * u.y + c * u.x, s * u.z * u.z - s};
}

Mat3 divided(const Mat3 &m, double d) {
  Mat3 r = m;
  for (double &v : r.m)
    v /= d;
  return r;
}

} // namespace

Metric metric_of(const Mat3 &B) { return packed(B.transpose() * B); }

CellConstraint CellConstraint::from_rotations(
    const std::vector<std::array<int, 9>> &rotations) {
  // One row per (rotation, element): the coefficients of R G* R^T - G* in the
  // six unknowns. Entries are small integers, so elimination in double is
  // exact.
  std::vector<std::array<double, 6>> rows;
  for (const auto &r : rotations) {
    for (int e = 0; e < 6; ++e) {
      std::array<double, 6> row{};
      for (int k = 0; k < 6; ++k) {
        // d(R G R^T)_{ij} / d g_k, with g_k the symmetric pair (p, q).
        const int i = kRow[e], j = kCol[e], p = kRow[k], q = kCol[k];
        double v = r[3 * i + p] * r[3 * j + q];
        if (p != q)
          v += r[3 * i + q] * r[3 * j + p];
        row[k] = v - (k == e ? 1.0 : 0.0);
      }
      rows.push_back(row);
    }
  }

  // Reduced row-echelon form, pivots taken in column order.
  std::vector<int> pivot_col;
  std::size_t lead = 0;
  for (int c = 0; c < 6 && lead < rows.size(); ++c) {
    std::size_t p = lead;
    while (p < rows.size() && std::abs(rows[p][c]) < 1e-9)
      ++p;
    if (p == rows.size())
      continue;
    std::swap(rows[p], rows[lead]);
    const double scale = rows[lead][c];
    for (double &v : rows[lead])
      v /= scale;
    for (std::size_t r = 0; r < rows.size(); ++r) {
      if (r == lead || rows[r][c] == 0.0)
        continue;
      const double f = rows[r][c];
      for (int k = 0; k < 6; ++k)
        rows[r][k] -= f * rows[lead][k];
    }
    pivot_col.push_back(c);
    ++lead;
  }

  CellConstraint out;
  for (int c = 0; c < 6; ++c) {
    bool is_pivot = false;
    for (int p : pivot_col)
      is_pivot = is_pivot || p == c;
    if (is_pivot)
      continue;
    Metric v{};
    v[c] = 1.0;
    for (std::size_t r = 0; r < pivot_col.size(); ++r)
      v[pivot_col[r]] = -rows[r][c];
    out.column.push_back(c);
    out.basis.push_back(v);
  }
  return out;
}

CellConstraint CellConstraint::from_hall(const std::string &hall) {
  return from_rotations(SpaceGroup::from_hall(hall).rotations());
}

std::vector<double> CellConstraint::independent(const Metric &g) const {
  std::vector<double> out(column.size());
  for (std::size_t k = 0; k < column.size(); ++k)
    out[k] = g[column[k]];
  return out;
}

Metric CellConstraint::metric(const std::vector<double> &independent) const {
  Metric g{};
  for (std::size_t k = 0; k < basis.size(); ++k)
    for (int e = 0; e < 6; ++e)
      g[e] += independent[k] * basis[k][e];
  return g;
}

Mat3 b_from_metric(const Metric &g) {
  // B^T B = G* with B lower triangular, solved from the last row up:
  // c* lies along z, b* in the y-z plane.
  if (!(g[2] > 0.0) || !(g[1] - g[5] * g[5] / g[2] > 0.0))
    throw ComposeError("the reciprocal metric is not positive definite");
  const double b22 = std::sqrt(g[2]);
  const double b21 = g[5] / b22;
  const double b20 = g[4] / b22;
  const double b11 = std::sqrt(g[1] - b21 * b21);
  const double b10 = (g[3] - b20 * b21) / b11;
  const double rest = g[0] - b10 * b10 - b20 * b20;
  if (!(rest > 0.0))
    throw ComposeError("the reciprocal metric is not positive definite");
  const double b00 = std::sqrt(rest);
  const Mat3 B{b00, 0.0, 0.0, b10, b11, 0.0, b20, b21, b22};
  // The direct cell's volume is 1 / det(B).
  if (!(1.0 / (b00 * b11 * b22) > 40.0))
    throw ComposeError("the direct cell's volume is 40 A^3 or less");
  return B;
}

Mat3 db_from_metric(const Mat3 &B, const Mat3 &B_inverse, const Metric &dg) {
  const Mat3 M = B_inverse.transpose() * symmetric(dg) * B_inverse;
  const Mat3 L{0.5 * M(0, 0), 0.0,     0.0,     M(1, 0),      0.5 * M(1, 1),
               0.0,           M(2, 0), M(2, 1), 0.5 * M(2, 2)};
  return L * B;
}

Orientation compose_orientation(const Mat3 &U0,
                                const std::array<double, 3> &mrad) {
  const Vec3 x{1.0, 0.0, 0.0}, y{0.0, 1.0, 0.0}, z{0.0, 0.0, 1.0};
  const double p1 = mrad[0] / 1000.0, p2 = mrad[1] / 1000.0,
               p3 = mrad[2] / 1000.0;
  const Mat3 Phi1 = axis_angle(x, p1), dPhi1 = axis_angle_derivative(x, p1);
  const Mat3 Phi2 = axis_angle(y, p2), dPhi2 = axis_angle_derivative(y, p2);
  const Mat3 Phi3 = axis_angle(z, p3), dPhi3 = axis_angle_derivative(z, p3);
  const Mat3 Phi21 = Phi2 * Phi1;
  Orientation o;
  o.U = (Phi3 * Phi21) * U0;
  o.dU[0] = divided(Phi3 * Phi2 * dPhi1 * U0, 1000.0);
  o.dU[1] = divided(Phi3 * dPhi2 * Phi1 * U0, 1000.0);
  o.dU[2] = divided(dPhi3 * Phi21 * U0, 1000.0);
  return o;
}

} // namespace mxi::stills

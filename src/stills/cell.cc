// Derived from cctbx_project (cctbx, rstbx, scitbx) and dxtbx.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.
// dxtbx: Copyright (c) 2019 Diamond Light Source, Lawrence Berkeley National
// Laboratory and the Science and Technology Facilities Council. All rights
// reserved. BSD-3-Clause; see third_party/dxtbx/LICENSE.txt.

#include "stills/cell.hh"

#include <cmath>

namespace mxi::stills {

namespace {

const double kPi180 = std::atan(1.) / 45.;

double acos_deg(double x) { return rad_as_deg(std::acos(x)); }

Sym6 construct_metric(const std::array<double, 6> &p,
                      const std::array<double, 3> &c) {
  return {p[0] * p[0],        p[1] * p[1],        p[2] * p[2],
          p[0] * p[1] * c[2], p[0] * p[2] * c[1], p[1] * p[2] * c[0]};
}

std::array<double, 6> parameters_from_metric(const Sym6 &g) {
  std::array<double, 6> p{};
  for (std::size_t i = 0; i < 3; ++i) {
    if (g[i] <= 0.)
      throw CellError("Corrupt metrical matrix.");
    p[i] = std::sqrt(g[i]);
  }
  p[3] = acos_deg(g[5] / p[1] / p[2]);
  p[4] = acos_deg(g[4] / p[2] / p[0]);
  p[5] = acos_deg(g[3] / p[0] / p[1]);
  return p;
}

void initialize(UnitCell &u) {
  for (std::size_t i = 0; i < 6; ++i)
    if (u.params[i] <= 0.)
      throw CellError("Unit cell parameter is zero or negative.");
  for (std::size_t i = 3; i < 6; ++i) {
    if (u.params[i] >= 180.)
      throw CellError(
          "Unit cell angle is greater than or equal to 180 degrees.");
    const double a_rad = deg_as_rad(u.params[i]);
    u.cos_ang[i - 3] = std::cos(a_rad);
    u.sin_ang[i - 3] = std::sin(a_rad);
    if (u.sin_ang[i - 3] == 0.)
      throw CellError("Unit cell angle is zero or extremely close to zero.");
  }

  const auto &p = u.params;
  const auto &c = u.cos_ang;
  const auto &s = u.sin_ang;
  double d = 1.;
  for (std::size_t i = 0; i < 3; ++i)
    d -= c[i] * c[i];
  d += 2. * c[0] * c[1] * c[2];
  if (d < 0.)
    throw CellError("Square of unit cell volume is negative.");
  u.volume = p[0] * p[1] * p[2] * std::sqrt(d);
  if (u.volume <= 0.)
    throw CellError("Unit cell volume is zero or negative.");

  for (std::size_t i = 0; i < 3; ++i)
    u.r_params[i] = p[(i + 1) % 3] * p[(i + 2) % 3] * s[i] / u.volume;
  for (std::size_t i = 0; i < 3; ++i) {
    const double denom = s[(i + 1) % 3] * s[(i + 2) % 3];
    if (denom == 0)
      throw CellError("Error computing reciprocal unit cell angles.");
    u.r_cos_ang[i] = (c[(i + 1) % 3] * c[(i + 2) % 3] - c[i]) / denom;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    if (u.r_cos_ang[i] < -1 || u.r_cos_ang[i] > 1)
      throw CellError("Error computing reciprocal unit cell angles.");
    const double a_rad = std::acos(u.r_cos_ang[i]);
    u.r_params[i + 3] = rad_as_deg(a_rad);
    u.r_sin_ang[i] = std::sin(a_rad);
    u.r_cos_ang[i] = std::cos(a_rad);
  }

  u.metric = construct_metric(p, c);

  const auto &rc = u.r_cos_ang;
  const double s1rca2 = std::sqrt(1. - rc[0] * rc[0]);
  if (s1rca2 == 0.)
    throw CellError("Reciprocal unit cell alpha angle is zero or extremely "
                    "close to zero.");
  u.orth = {p[0], c[2] * p[1], c[1] * p[2],
            0.,   s[2] * p[1], -s[1] * rc[0] * p[2],
            0.,   0.,          s[1] * p[2] * s1rca2};
  u.frac = {1. / p[0],
            -c[2] / (s[2] * p[0]),
            -(c[2] * s[1] * rc[0] + c[1] * s[2]) /
                (s[1] * s1rca2 * s[2] * p[0]),
            0.,
            1. / (s[2] * p[1]),
            rc[0] / (s1rca2 * s[2] * p[1]),
            0.,
            0.,
            1. / (s[1] * s1rca2 * p[2])};
}

double length(const Vec3 &v) { return std::sqrt(v.dot(v)); }

double angle(const Vec3 &u, const Vec3 &v) {
  return std::acos(u.dot(v) / (length(u) * length(v)));
}

} // namespace

double deg_as_rad(double deg) { return deg * kPi180; }
double rad_as_deg(double rad) { return rad / kPi180; }

Mat3 inverse_scitbx(const Mat3 &d) {
  const auto &m = d.m;
  const double det = d.determinant();
  if (det == 0.)
    throw CellError("Matrix is not invertible.");
  return {(m[4] * m[8] - m[5] * m[7]) / det, (-m[1] * m[8] + m[2] * m[7]) / det,
          (m[1] * m[5] - m[2] * m[4]) / det, (-m[3] * m[8] + m[5] * m[6]) / det,
          (m[0] * m[8] - m[2] * m[6]) / det, (-m[0] * m[5] + m[2] * m[3]) / det,
          (m[3] * m[7] - m[4] * m[6]) / det, (-m[0] * m[7] + m[1] * m[6]) / det,
          (m[0] * m[4] - m[1] * m[3]) / det};
}

UnitCell UnitCell::from_parameters(const std::array<double, 6> &params) {
  UnitCell u;
  u.params = params;
  initialize(u);
  return u;
}

UnitCell UnitCell::from_metric(const Sym6 &g) {
  UnitCell u;
  u.params = parameters_from_metric(g);
  initialize(u);
  return u;
}

UnitCell UnitCell::from_orthogonalisation(const Mat3 &o) {
  const auto &m = o.m;
  return from_metric({m[0] * m[0] + m[3] * m[3] + m[6] * m[6],
                      m[1] * m[1] + m[4] * m[4] + m[7] * m[7],
                      m[2] * m[2] + m[5] * m[5] + m[8] * m[8],
                      m[0] * m[1] + m[3] * m[4] + m[6] * m[7],
                      m[0] * m[2] + m[3] * m[5] + m[6] * m[8],
                      m[1] * m[2] + m[4] * m[5] + m[7] * m[8]});
}

CrystalModel CrystalModel::from_real_space(const Vec3 &a, const Vec3 &b,
                                           const Vec3 &c) {
  const Mat3 A = inverse_scitbx(Mat3::from_rows(a, b, c));
  CrystalModel x;
  x.cell = UnitCell::from_parameters(
      {length(a), length(b), length(c), rad_as_deg(angle(b, c)),
       rad_as_deg(angle(c, a)), rad_as_deg(angle(a, b))});
  x.B = x.cell.frac.transpose();
  x.U = A * inverse_scitbx(x.B);
  return x;
}

CrystalModel CrystalModel::from_A(const Mat3 &A) {
  const UnitCell uc =
      UnitCell::from_orthogonalisation(inverse_scitbx(A.transpose()));
  const Mat3 B = uc.frac.transpose();
  CrystalModel x;
  x.U = A * inverse_scitbx(B);
  // set_B: the cell of crystal_orientation(B, reciprocal), from the rows of
  // the direct matrix B^-1.
  const Mat3 direct = inverse_scitbx(B);
  const Vec3 p0 = direct.row(0), p1 = direct.row(1), p2 = direct.row(2);
  x.cell = UnitCell::from_metric(
      {p0.dot(p0), p1.dot(p1), p2.dot(p2), p0.dot(p1), p2.dot(p0), p1.dot(p2)});
  x.B = x.cell.frac.transpose();
  return x;
}

CrystalModel CrystalModel::change_basis(const Mat3 &M) const {
  return from_A(A() * M);
}

} // namespace mxi::stills

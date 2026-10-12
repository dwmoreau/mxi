// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "predict.hh"

#include <cmath>
#include <math.h>

namespace mxi::stills {

namespace {

Vec3 unit(const Vec3 &v) { return v / v.norm(); }

// The derivative of q rotated about the unit axis u by angle, for a change du
// in the axis.
Vec3 rotate_axis_derivative(const Vec3 &q, const Vec3 &u, double angle,
                            const Vec3 &du) {
  const double c = std::cos(angle), s = std::sin(angle);
  return du * (u.dot(q) * (1.0 - c)) + u * (du.dot(q) * (1.0 - c)) +
         du.cross(q) * s;
}

} // namespace

Beam::Beam(const Vec3 &s0_)
    : s0(s0_), unit_s0(unit(s0_)), s0_length(s0_.norm()),
      wavelength(1.0 / s0_.norm()) {}

Vec3 rotate_about_unit(const Vec3 &q, const Vec3 &u, double angle) {
  const double c = std::cos(angle);
  return q * c + u * (u.dot(q)) * (1.0 - c) + u.cross(q) * std::sin(angle);
}

PredictStatus predict(const Beam &beam, const Mat3 &A,
                      const std::array<int, 3> &h, const Mat3 &D, Ray *out) {
  const Vec3 q = A * Vec3{static_cast<double>(h[0]), static_cast<double>(h[1]),
                          static_cast<double>(h[2])};
  if (!(q.norm() > 0.0))
    return PredictStatus::zero_q;
  const Vec3 e1 = unit(q.cross(beam.unit_s0));
  const Vec3 c0 = unit(beam.unit_s0.cross(e1));

  const double qq = q.norm_squared();
  const double lambda = 1.0 / beam.s0_length;
  const double a = 0.5 * qq * lambda;
  const double tmp = qq - a * a;
  if (!(tmp > 0.0))
    return PredictStatus::beyond_sphere;
  const double b = std::sqrt(tmp);
  const Vec3 r = beam.unit_s0 * (-1.0 * a) + c0 * b;

  const Vec3 q0 = unit(q);
  const Vec3 q1 = unit(q0.cross(e1));
  out->delpsi = -1.0 * std::atan2(r.dot(q1), r.dot(q0));
  out->s1 = unit(beam.s0 + r) * beam.s0_length;
  out->q = q;
  out->r = r;
  out->e1 = e1;

  const Vec3 v = D * out->s1;
  if (v.z <= 0.0)
    return PredictStatus::misses_panel;
  out->pv = v;
  out->x_mm = v.x / v.z;
  out->y_mm = v.y / v.z;
  return PredictStatus::ok;
}

RayFrame ray_frame(const Beam &beam, const Ray &ray) {
  // The derivative code takes r from the predicted s1, and e1 from q0 x s0,
  // as dials' stills prediction parameterisation does.
  RayFrame f;
  f.r = ray.s1 - beam.s0;
  f.q_scalar = ray.q.norm();
  f.qq = f.q_scalar * f.q_scalar;
  f.q0 = unit(ray.q);
  const Vec3 n = f.q0.cross(beam.unit_s0);
  f.e1 = unit(n);
  f.n_norm = n.norm();
  ::sincos(ray.delpsi, &f.sin_delpsi, &f.cos_delpsi);
  f.one_minus_cos = 1.0 - f.cos_delpsi;
  f.e1_r = f.e1.cross(f.r);
  f.e1_r_dot_s0 = f.e1_r.dot(beam.s0);
  f.e1_dot_q_versin = f.e1.dot(ray.q) * f.one_minus_cos;
  f.w_inv = 1.0 / ray.pv.z;
  f.u_w_inv = ray.pv.x * f.w_inv;
  f.v_w_inv = ray.pv.y * f.w_inv;
  return f;
}

RayDerivative ray_derivative(const Beam &beam, const Ray &ray,
                             const RayFrame &f, const Mat3 &D, const Vec3 &dq,
                             AxisDerivative mode) {
  const Vec3 &e1 = f.e1;
  const double c = f.cos_delpsi, s = f.sin_delpsi, versin = f.one_minus_cos;
  const Vec3 dr = dq * c + e1 * (e1.dot(dq)) * versin + e1.cross(dq) * s;
  const double ddelpsi = -1.0 * dr.dot(ray.s1) / f.e1_r_dot_s0;

  const double q_dot_dq = ray.q.dot(dq);
  const Vec3 dq0 = (dq * f.q_scalar - f.q0 * q_dot_dq) / f.qq;
  const Vec3 de1 = dq0.cross(beam.unit_s0);

  Vec3 drde;
  if (mode == AxisDerivative::dials) {
    const double step = 1.0e-8;
    const Vec3 del_e1 = de1 * step;
    const Vec3 fwd =
        rotate_about_unit(ray.q, unit(e1 + del_e1 * 0.5), ray.delpsi);
    const Vec3 rev =
        rotate_about_unit(ray.q, unit(e1 - del_e1 * 0.5), ray.delpsi);
    drde = (fwd - rev) * (1.0 / step);
  } else {
    // e1 = n / |n| with n = q0 x s0_unit, so de1 = (dn - e1 (e1 . dn)) / |n|.
    Vec3 du = de1 - e1 * e1.dot(de1);
    if (mode == AxisDerivative::analytic)
      du = du / f.n_norm;
    drde = du * f.e1_dot_q_versin + e1 * (du.dot(ray.q) * versin) +
           du.cross(ray.q) * s;
  }

  const Vec3 dpv = D * (dr + f.e1_r * ddelpsi + drde);
  RayDerivative out;
  out.dx = f.w_inv * (dpv.x - dpv.z * f.u_w_inv);
  out.dy = f.w_inv * (dpv.y - dpv.z * f.v_w_inv);
  out.ddelpsi = ddelpsi;
  return out;
}

RayDerivative ray_derivative(const Beam &beam, const Ray &ray, const Mat3 &D,
                             const Vec3 &dq, AxisDerivative mode) {
  return ray_derivative(beam, ray, ray_frame(beam, ray), D, dq, mode);
}

} // namespace mxi::stills

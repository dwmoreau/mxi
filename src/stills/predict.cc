#include "predict.hh"

#include <cmath>

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

RayDerivative ray_derivative(const Beam &beam, const Ray &ray, const Mat3 &D,
                             const Vec3 &dq, AxisDerivative mode) {
  // The derivative code takes r from the predicted s1, and e1 from q0 x s0,
  // as dials' stills prediction parameterisation does.
  const Vec3 r = ray.s1 - beam.s0;
  const double q_scalar = ray.q.norm();
  const double qq = q_scalar * q_scalar;
  const Vec3 q0 = unit(ray.q);
  const Vec3 e1 = unit(q0.cross(beam.unit_s0));

  const Vec3 dr = rotate_about_unit(dq, e1, ray.delpsi);
  const Vec3 e1_r = e1.cross(r);
  const double ddelpsi = -1.0 * dr.dot(ray.s1) / e1_r.dot(beam.s0);

  const double q_dot_dq = ray.q.dot(dq);
  const Vec3 dq0 = (dq * q_scalar - q0 * q_dot_dq) / qq;
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
    Vec3 de1_unit = de1 - e1 * e1.dot(de1);
    if (mode == AxisDerivative::analytic)
      de1_unit = de1_unit / q0.cross(beam.unit_s0).norm();
    drde = rotate_axis_derivative(ray.q, e1, ray.delpsi, de1_unit);
  }

  const Vec3 dpv = D * (dr + e1_r * ddelpsi + drde);
  const double w_inv = 1.0 / ray.pv.z;
  const double u_w_inv = ray.pv.x * w_inv;
  const double v_w_inv = ray.pv.y * w_inv;
  RayDerivative out;
  out.dx = w_inv * (dpv.x - dpv.z * u_w_inv);
  out.dy = w_inv * (dpv.y - dpv.z * v_w_inv);
  out.ddelpsi = ddelpsi;
  return out;
}

} // namespace mxi::stills

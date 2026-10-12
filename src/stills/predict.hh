// Stills prediction: the reciprocal lattice point q = A h is rotated about
// e1 = q x s0 onto the Ewald sphere, by the angle Delta-Psi (Sauter et al.
// 2014, eqs. 4-12), and the diffracted ray s1 = s0 + r meets the detector at
// (X, Y). The derivatives of X, Y and Delta-Psi with respect to a change dq in
// the lattice point follow Waterman et al. (2016) eqs. 40-45 with Delta-Psi in
// place of the rotation angle.
//
// The arithmetic follows dials' StillsRayPredictor and its stills prediction
// parameterisation, in the same order.
//
// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <array>

#include "linalg.hh"

namespace mxi::stills {

// How the derivative of r with respect to the rotation axis e1 is formed.
//   dials        as dials does: a central difference, step 1e-8, of q
//                rotated about e1 +- step/2 (dq0 x s0_unit), the axis
//                normalised.
//   dials_limit  the limit that difference approximates, without its
//                rounding error: the axis moves by the part of
//                dq0 x s0_unit perpendicular to e1.
//   analytic     the derivative of q rotated about normalise(q0 x s0_unit),
//                which divides that part by |q0 x s0_unit|.
enum class AxisDerivative { dials, dials_limit, analytic };

// One reflection's prediction, with the vectors its derivatives need.
struct Ray {
  Vec3 q;  // A h
  Vec3 e1; // the rotation axis, unit
  Vec3 r;  // q rotated onto the Ewald sphere
  Vec3 s1; // the diffracted beam
  double delpsi = 0.0;
  Vec3 pv; // D s1
  double x_mm = 0.0, y_mm = 0.0;
};

struct Beam {
  Vec3 s0;
  Vec3 unit_s0;
  double s0_length = 0.0;
  double wavelength = 0.0; // 1 / |s0|

  explicit Beam(const Vec3 &s0_);
};

enum class PredictStatus { ok, zero_q, beyond_sphere, misses_panel };

// Predicts reflection h of the crystal A = U B on the panel with D-matrix D.
PredictStatus predict(const Beam &beam, const Mat3 &A,
                      const std::array<int, 3> &h, const Mat3 &D, Ray *out);

struct RayDerivative {
  double dx = 0.0, dy = 0.0, ddelpsi = 0.0;
};

// The derivatives of X, Y and Delta-Psi for a change dq in q, at a predicted
// ray.
RayDerivative ray_derivative(const Beam &beam, const Ray &ray, const Mat3 &D,
                             const Vec3 &dq, AxisDerivative mode);

// q rotated about the unit axis u by angle, as scitbx's
// unit_rotate_around_origin forms it.
Vec3 rotate_about_unit(const Vec3 &q, const Vec3 &u, double angle);

} // namespace mxi::stills

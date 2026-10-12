// Refinement of one still's crystal against its indexed spots: the orientation
// (three angles) and the cell (the metric elements its space group leaves
// free), with the beam and the detector fixed. The target, weights and
// minimiser are those of dials' stills refinement (Brewster et al. 2018, eq.
// 5):
//
//   L = 1/2 sum_i w_x (X - X_obs)^2 + w_y (Y - Y_obs)^2 + w_psi DeltaPsi^2,
//
// with w_x, w_y the inverse centroid variances and w_psi = 1e6, minimised by
// Levenberg-Marquardt as dials' LevenbergMarquardtIterations does it.
//
// Derived from dials and dxtbx.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// dxtbx: Copyright (c) 2019 Diamond Light Source, Lawrence Berkeley National
// Laboratory and the Science and Technology Facilities Council. All rights
// reserved. BSD-3-Clause; see third_party/dxtbx/LICENSE.txt.

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "crystal.hh"
#include "geometry.hh"
#include "predict.hh"

namespace mxi::stills {

// A still's indexed spots, in millimetres on their panels.
struct IndexedSpots {
  std::vector<std::array<int, 3>> miller;
  std::vector<std::int64_t> panel;
  std::vector<double> x_mm, y_mm;
  std::vector<double> var_x_mm, var_y_mm;
  std::vector<std::uint64_t> flags; // dials' reflection flags

  std::size_t size() const { return miller.size(); }
};

struct RefineOptions {
  AxisDerivative axis_derivative = AxisDerivative::dials;
  int max_iterations = 100;
  // A relative change in every RMSD below this between two accepted steps
  // ends the refinement.
  double rmsd_tolerance = 1e-4;
  double delpsi_weight = 1e6;
};

enum class RefineStatus {
  rmsd_converged,
  max_iterations,
  max_trial_iterations, // nu reached 8192 without an accepted step
  objective_flat,
  too_few_degrees_of_freedom,
  // The failures; the crystal is unchanged.
  no_reflections,
  zero_weights,
  indeterminate, // a parameter's gradient is exactly zero
  not_positive_definite,
  cannot_compose,
  prediction_failed,
};

const char *to_string(RefineStatus status);
bool failed(RefineStatus status);

// One row per accepted state: the start and after each accepted step.
struct JournalRow {
  std::array<double, 3> rmsd{}; // X mm, Y mm, DeltaPsi rad
  double objective = 0.0;
  double mu = 0.0;
  int nu = 0;
  std::vector<double> parameters;
};

struct RefineResult {
  RefineStatus status = RefineStatus::no_reflections;
  std::string message;
  Mat3 U, B;
  std::vector<std::string> parameter_names;
  int n_iterations = 0; // trial steps, accepted or not
  std::size_t n_used = 0;
  std::array<double, 3> rmsd{};
  std::vector<JournalRow> journal;
  // Per used spot (input order): the index into the input, and its final
  // prediction.
  std::vector<std::size_t> used;
  std::vector<double> x_calc, y_calc, delpsi;
  std::vector<Vec3> s1;
};

// Each spot's prediction for the crystal A = U B, in input order; status per
// spot is that of predict().
struct Predictions {
  std::vector<PredictStatus> status;
  std::vector<double> x_mm, y_mm, delpsi;
  std::vector<Vec3> s1;
};
Predictions predict_spots(const StillsDetector &detector, const Vec3 &s0,
                          const Mat3 &U, const Mat3 &B,
                          const std::vector<std::array<int, 3>> &miller,
                          const std::vector<std::int64_t> &panel);

// The weighted least-squares problem at the starting model, before any
// parameter is fixed: the residuals (X, Y, Delta-Psi per used spot), their
// weights, and the Jacobian, row major, one column per parameter.
struct Linearisation {
  RefineStatus status = RefineStatus::rmsd_converged;
  std::string message;
  std::vector<std::string> parameter_names;
  std::vector<double> parameters;
  std::vector<std::size_t> used;
  std::vector<double> residual, weight, jacobian;
};
Linearisation linearise_crystal(const StillsDetector &detector, const Vec3 &s0,
                                const Mat3 &U, const Mat3 &B,
                                const std::string &hall,
                                const IndexedSpots &spots,
                                const RefineOptions &options = {});

// U and B as dxtbx gives them; hall names the space group whose point group
// constrains the cell.
RefineResult refine_crystal(const StillsDetector &detector, const Vec3 &s0,
                            const Mat3 &U, const Mat3 &B,
                            const std::string &hall, const IndexedSpots &spots,
                            const RefineOptions &options = {});
RefineResult refine_crystal(const StillsDetector &detector, const Vec3 &s0,
                            const Mat3 &U, const Mat3 &B,
                            const CellConstraint &constraint,
                            const IndexedSpots &spots,
                            const RefineOptions &options = {});

} // namespace mxi::stills

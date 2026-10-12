// Derived from dials.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.

#include "stills/candidates.hh"

#include <algorithm>
#include <cmath>
#include <map>

#include "stills/reduction.hh"

namespace mxi::stills {

namespace {

double clipped_angle(double dot, double norm_sq_product) {
  return std::acos(std::clamp(dot / std::sqrt(norm_sq_product), -1.0, 1.0));
}

int iround(double x) {
  if (x < 0)
    return static_cast<int>(x - 0.5);
  return static_cast<int>(x + .5);
}

} // namespace

std::vector<CrystalModel>
candidate_orientation_matrices(const std::vector<Vec3> &basis_vectors) {
  std::vector<CrystalModel> out;
  const int n =
      static_cast<int>(std::min<std::size_t>(basis_vectors.size(), 100));
  if (n < 3)
    return out;

  std::vector<std::array<int, 3>> triplets;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k)
        triplets.push_back({i, j, k});
  std::stable_sort(triplets.begin(), triplets.end(),
                   [](const auto &x, const auto &y) {
                     return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] <
                            y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
                   });

  const double half_pi = 0.5 * M_PI;
  const double min_angle = 20. / 180. * M_PI;
  for (const auto &t : triplets) {
    Vec3 a = basis_vectors[t[0]], b = basis_vectors[t[1]],
         c = basis_vectors[t[2]];
    const double a_ns = a.dot(a), b_ns = b.dot(b), c_ns = c.dot(c);

    const double angle_ab = clipped_angle(a.dot(b), a_ns * b_ns);
    if (!(angle_ab >= min_angle && (M_PI - angle_ab) >= min_angle))
      continue;
    if (angle_ab < half_pi)
      b = -b;

    const Vec3 a_cross_b = a.cross(b);
    const double acb_ns = a_cross_b.dot(a_cross_b);
    const double angle_acb_c = clipped_angle(a_cross_b.dot(c), acb_ns * c_ns);
    if (!(std::abs(half_pi - angle_acb_c) >= min_angle))
      continue;

    const double alpha = clipped_angle(b.dot(c), b_ns * c_ns);
    if (alpha < half_pi)
      c = -c;
    if (a_cross_b.dot(c) < 0) {
      a = -a;
      b = -b;
      c = -c;
    }

    const CrystalModel model = CrystalModel::from_real_space(a, b, c);
    const CrystalModel reduced = model.change_basis(
        niggli_change_of_basis(niggli_reduction(model.cell)));
    const auto &p = reduced.cell.params;
    if (reduced.cell.volume > (p[0] * p[1] * p[2] / 100))
      out.push_back(reduced);
  }
  return out;
}

Assignment assign_indices_global(const std::vector<Vec3> &rlp,
                                 const std::vector<int> &id,
                                 const std::vector<Mat3> &A, double d_min,
                                 double tolerance) {
  const std::size_t n_all = rlp.size();
  Assignment out;
  out.miller_index.assign(n_all, {0, 0, 0});
  out.id.assign(n_all, -1);

  std::vector<std::size_t> rows;
  for (std::size_t i = 0; i < n_all; ++i) {
    if (id[i] != -1)
      continue;
    if (d_min > 0 && !(1 / std::sqrt(rlp[i].dot(rlp[i])) > d_min))
      continue;
    rows.push_back(i);
  }
  const std::size_t n = rows.size();

  std::vector<std::vector<std::array<int, 3>>> hkl_ints(A.size());
  std::vector<std::vector<double>> lengths_sq(A.size());
  for (std::size_t l = 0; l < A.size(); ++l) {
    const Mat3 A_inv = inverse_scitbx(A[l]);
    hkl_ints[l].resize(n);
    lengths_sq[l].resize(n);
    for (std::size_t r = 0; r < n; ++r) {
      const Vec3 hkl_f = A_inv * rlp[rows[r]];
      std::array<int, 3> h{};
      for (int j = 0; j < 3; ++j)
        h[j] = iround(hkl_f[j]);
      const Vec3 diff = hkl_f - Vec3(h[0], h[1], h[2]);
      hkl_ints[l][r] = h;
      lengths_sq[l][r] = diff.dot(diff);
    }
  }

  std::vector<std::array<int, 3>> miller(n, {0, 0, 0});
  std::vector<int> crystal(n, -1);
  std::multimap<std::array<int, 3>, std::size_t> by_hkl;
  const double tolerance_sq = tolerance * tolerance;
  for (std::size_t r = 0; r < n; ++r) {
    if (A.empty())
      break;
    std::size_t best = 0;
    for (std::size_t l = 1; l < A.size(); ++l)
      if (lengths_sq[l][r] < lengths_sq[best][r])
        best = l;
    if (lengths_sq[best][r] > tolerance_sq)
      continue;
    const auto &h = hkl_ints[best][r];
    if (h[0] == 0 && h[1] == 0 && h[2] == 0)
      continue;
    miller[r] = h;
    by_hkl.emplace(h, r);
    crystal[r] = static_cast<int>(best);
  }

  for (auto it = by_hkl.begin(); it != by_hkl.end();) {
    auto end = by_hkl.upper_bound(it->first);
    std::vector<std::size_t> same;
    for (auto jt = it; jt != end; ++jt)
      same.push_back(jt->second);
    for (std::size_t i = 0; i < same.size(); ++i)
      for (std::size_t j = i + 1; j < same.size(); ++j) {
        const std::size_t ri = same[i], rj = same[j];
        const int ci = crystal[ri], cj = crystal[rj];
        if (ci != cj || ci == -1)
          continue;
        if (lengths_sq[cj][rj] < lengths_sq[ci][ri]) {
          miller[ri] = {0, 0, 0};
          crystal[ri] = -1;
        } else {
          miller[rj] = {0, 0, 0};
          crystal[rj] = -1;
        }
      }
    it = end;
  }

  for (std::size_t r = 0; r < n; ++r) {
    out.miller_index[rows[r]] = miller[r];
    if (crystal[r] >= 0)
      out.id[rows[r]] = crystal[r];
  }
  return out;
}

} // namespace mxi::stills

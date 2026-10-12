// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "stills/rsgs.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <set>

#include "stills/fft1d.hh"

namespace mxi::stills {

namespace {

// numpy's pairwise summation of a contiguous float64 row: eight accumulators
// over blocks of up to 128, halves above that, starting from -0.0 below 8.
double pairwise_sum(const double *a, std::size_t n) {
  if (n < 8) {
    double res = -0.0;
    for (std::size_t i = 0; i < n; ++i)
      res += a[i];
    return res;
  }
  if (n <= 128) {
    double r[8];
    for (int k = 0; k < 8; ++k)
      r[k] = a[k];
    std::size_t i = 8;
    for (; i < n - (n % 8); i += 8)
      for (int k = 0; k < 8; ++k)
        r[k] += a[i + static_cast<std::size_t>(k)];
    double res =
        ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; ++i)
      res += a[i];
    return res;
  }
  std::size_t n2 = n / 2;
  n2 -= n2 % 8;
  return pairwise_sum(a, n2) + pairwise_sum(a + n2, n - n2);
}

// cos(2 pi x). x is reduced exactly to r = x - q / 4 with q = round(4 x), so
// |r| <= 1/8, and the cosine or sine of 2 pi r, by Taylor series to within
// 1e-17, is taken as q mod 4 requires: cos, -sin, -cos, sin. The two series
// are one Horner evaluation with the coefficients of whichever the quadrant
// needs, the sine's ending a term early with a zero coefficient.
inline double cos_two_pi(double x) {
  constexpr double kShift = 6755399441055744.0; // 1.5 * 2^52: rounds to integer
  const double shifted = 4.0 * x + kShift;
  const double q = shifted - kShift;
  std::uint64_t bits;
  std::memcpy(&bits, &shifted, sizeof bits);
  const double y = (x - 0.25 * q) * (2 * M_PI);
  const double z = y * y;
  const unsigned quadrant = static_cast<unsigned>(bits & 3u);
  const bool odd = (quadrant & 1u) != 0;
  const double a8 = odd ? 0.0 : 1.0 / 20922789888000.0;
  const double a7 = odd ? -1.0 / 1307674368000.0 : -1.0 / 87178291200.0;
  const double a6 = odd ? 1.0 / 6227020800.0 : 1.0 / 479001600;
  const double a5 = odd ? -1.0 / 39916800 : -1.0 / 3628800;
  const double a4 = odd ? 1.0 / 362880 : 1.0 / 40320;
  const double a3 = odd ? -1.0 / 5040 : -1.0 / 720;
  const double a2 = odd ? 1.0 / 120 : 1.0 / 24;
  const double a1 = odd ? -1.0 / 6 : -1.0 / 2;
  const double p =
      1.0 +
      z * (a1 + z * (a2 + z * (a3 + z * (a4 + z * (a5 + z * (a6 + z * (a7 + z * a8)))))));
  const double v = odd ? y * p : p;
  return ((quadrant + 1u) & 2u) ? -v : v;
}

// scitbx.matrix's length, dot and angle (degrees), as Python evaluates them.
double norm_sq(const Vec3 &v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
double length(const Vec3 &v) { return std::sqrt(norm_sq(v)); }

// None (here NaN) where a vector is zero, as scitbx returns it.
double angle_deg(const Vec3 &a, const Vec3 &b) {
  const double a_sq = norm_sq(a), b_sq = norm_sq(b);
  const double d = a_sq * b_sq;
  if (a_sq == 0 || b_sq == 0 || d == 0)
    return std::nan("");
  const double cos_angle = a.dot(b) / std::sqrt(d);
  return std::acos(std::max(-1.0, std::min(1.0, cos_angle))) * (180 / M_PI);
}

struct VectorGroup {
  std::vector<Vec3> vectors;
  std::vector<double> weights;
  // The running mean: sums from 0 in member order, then each divided by the
  // member count.
  Vec3 mean() const {
    double sx = 0, sy = 0, sz = 0;
    for (const Vec3 &v : vectors) {
      sx += v[0];
      sy += v[1];
      sz += v[2];
    }
    const double n = static_cast<double>(vectors.size());
    return {sx / n, sy / n, sz / n};
  }
};

} // namespace

RealSpaceGridSearch::RealSpaceGridSearch(const std::array<double, 6> &cell,
                                         double characteristic_grid) {
  const std::set<double> unique_dims(cell.begin(), cell.begin() + 3);
  const std::vector<Vec3> directions = hemisphere_grid(characteristic_grid);
  sv_.reserve(unique_dims.size() * directions.size());
  for (double l : unique_dims)
    for (const Vec3 &d : directions)
      sv_.push_back(Vec3(d[0] * l, d[1] * l, d[2] * l));
}

std::vector<double>
RealSpaceGridSearch::scores(const std::vector<Vec3> &rlp) const {
  const std::size_t n = rlp.size();
  std::vector<double> bx(n), by(n), bz(n);
  for (std::size_t j = 0; j < n; ++j) {
    bx[j] = rlp[j][0];
    by[j] = rlp[j][1];
    bz[j] = rlp[j][2];
  }
  std::vector<double> out(sv_.size());
  std::vector<double> row(n);
  for (std::size_t i = 0; i < sv_.size(); ++i) {
    const double ax = sv_[i][0], ay = sv_[i][1], az = sv_[i][2];
    for (std::size_t j = 0; j < n; ++j)
      row[j] = cos_two_pi(ax * bx[j] + ay * by[j] + az * bz[j]);
    out[i] = pairwise_sum(row.data(), n);
  }
  return out;
}

std::vector<Vec3>
RealSpaceGridSearch::find_basis_vectors(const std::vector<Vec3> &rlp,
                                        std::size_t max_vectors) const {
  const std::vector<double> w = scores(rlp);
  std::vector<std::size_t> perm(w.size());
  for (std::size_t i = 0; i < perm.size(); ++i)
    perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(),
                   [&](std::size_t a, std::size_t b) { return w[a] > w[b]; });

  // utils.group_vectors(vectors, weights, max_groups=max_vectors).
  const double relative_length_tolerance = 0.1, angular_tolerance = 5;
  std::vector<VectorGroup> groups;
  for (std::size_t p : perm) {
    if (groups.size() == max_vectors)
      break;
    const Vec3 &v = sv_[p];
    const double len = length(v);
    bool matched = false;
    for (VectorGroup &g : groups) {
      const Vec3 mean_v = g.mean();
      const double mean_len = length(mean_v);
      if (std::abs(mean_len - len) / std::max(mean_len, len) <
          relative_length_tolerance) {
        const double angle = angle_deg(mean_v, v);
        if (angle < angular_tolerance) {
          g.vectors.push_back(v);
          g.weights.push_back(w[p]);
          matched = true;
          break;
        }
        if (std::abs(180 - angle) < angular_tolerance) {
          g.vectors.push_back(Vec3(-v[0], -v[1], -v[2]));
          g.weights.push_back(w[p]);
          matched = true;
          break;
        }
      }
    }
    if (!matched) {
      VectorGroup g;
      g.vectors.push_back(v);
      g.weights.push_back(w[p]);
      groups.push_back(std::move(g));
    }
  }
  // The member with the largest weight, the first on a tie.
  std::vector<Vec3> out;
  for (const VectorGroup &g : groups) {
    const auto best = std::max_element(g.weights.begin(), g.weights.end());
    out.push_back(
        g.vectors[static_cast<std::size_t>(best - g.weights.begin())]);
  }
  return out;
}

} // namespace mxi::stills

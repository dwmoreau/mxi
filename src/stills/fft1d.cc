// Derived from dials and cctbx_project (cctbx, rstbx, scitbx).
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "stills/fft1d.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <map>
#include <math.h>
#include <utility>

#include "stills/fftpack.hh"

namespace mxi::stills {

namespace {

const double kPi = 4. * std::atan(1.);
const double kPi_2 = 2. * std::atan(1.);
const double kPi180 = std::atan(1.) / 45.;

// The granularity n of Steller's eq. 6, and the first coefficient rstbx
// considers for the peak.
constexpr double kGranularity = 5;
constexpr std::size_t kF0Cutoff = 3;
// rstbx's kval: coefficients whose squared magnitude is below the running
// best's by more than this fraction are not compared by magnitude.
constexpr double kKvalTieMargin = 1.0e-13;

// libm's pow, called as CPython's math.pow calls it: not folded into x * x.
double py_pow(double x, double y) {
  static double (*volatile pow_fn)(double, double) = ::pow;
  return pow_fn(x, y);
}

// rstbx's compiled code takes each sine and cosine of one angle from a single
// sincos call, which can differ from sin and cos in the last bit; so do these.
std::pair<double, double> sin_cos(double x) {
  double s, c;
  ::sincos(x, &s, &c);
  return {s, c};
}

double length(const Vec3 &v) { return std::sqrt(v.dot(v)); }
Vec3 normalize(const Vec3 &v) { return v / length(v); }

// scitbx's vec3::unit_rotate_around_origin.
Vec3 unit_rotate(const Vec3 &v, const Vec3 &unit, double angle) {
  const auto [sinang, cosang] = sin_cos(angle);
  return v * cosang + unit * unit.dot(v) * (1.0 - cosang) +
         unit.cross(v) * sinang;
}

// rstbx's Direction.
struct Direction {
  double psi = 0, phi = 0;
  Vec3 dvec;
  std::size_t kmax = 1000000;
  double kval = 1.0e300;
  std::size_t m = 1000000;
  double delta_p = 1.0e300, pmin = 1.0e300, uc_length = 1.0e300;

  Direction() = default;
  Direction(double psi_, double phi_) : psi(psi_), phi(phi_) {
    const auto [sin_psi, cos_psi] = sin_cos(psi);
    const auto [sin_phi, cos_phi] = sin_cos(phi);
    dvec = Vec3(sin_psi * cos_phi, sin_psi * sin_phi, cos_psi);
  }
  explicit Direction(const Vec3 &p) {
    dvec = p;
    uc_length = length(dvec);
    dvec = normalize(dvec);
    phi = std::atan2(dvec[1], dvec[0]);
    psi = std::acos(dvec[2]);
  }
};

bool is_nearly_collinear(const Direction &a, const Direction &b) {
  const double dot = a.dvec.dot(b.dvec);
  const double inner_angle = std::acos(std::min(1.0, std::abs(dot)));
  return inner_angle < 10.0 * kPi180 || inner_angle > 170.0 * kPi180;
}

struct KvalCmp {
  bool operator()(const Direction &a, const Direction &b) const {
    if (a.kval == b.kval)
      return false;
    return a.kval > b.kval;
  }
};

std::size_t primecheck(std::size_t input_m) {
  std::size_t residual = input_m;
  std::size_t return_m = input_m;
  while (residual % 2 == 0)
    residual /= 2;
  while (residual % 3 == 0)
    residual /= 3;
  while (residual % 5 == 0)
    residual /= 5;
  if (residual > 1)
    return_m = primecheck(return_m + 1);
  return return_m;
}

// rstbx's dps_core with its directional FFT and workspace.
class DpsCore {
public:
  DpsCore(const std::vector<Vec3> &xyz, double amax, bool single_precision)
      : xyz_(xyz), amax_(amax), single_precision_(single_precision) {}

  std::size_t size() const { return xyz_.size(); }

  // Directional_FFT followed by extract_directional_properties, for every
  // direction in dirs. Directions whose transforms have the same length are
  // transformed kLanes at a time, one in each lane of fftpack::Lanes.
  void evaluate_all(std::vector<Direction> &dirs) {
    const std::size_t n_dirs = dirs.size();
    extents_.resize(n_dirs);
    for (std::size_t i = 0; i < n_dirs; ++i)
      extents_[i] = extent(dirs[i]);
    order_.resize(n_dirs);
    for (std::size_t i = 0; i < n_dirs; ++i)
      order_[i] = i;
    std::stable_sort(order_.begin(), order_.end(),
                     [&](std::size_t a, std::size_t b) {
                       return extents_[a].m < extents_[b].m;
                     });
    const std::size_t width =
        single_precision_ ? fftpack::kFloatLanes : fftpack::kLanes;
    for (std::size_t a = 0; a < n_dirs;) {
      const std::size_t m = extents_[order_[a]].m;
      std::size_t b = a + 1;
      while (b < n_dirs && b - a < width && extents_[order_[b]].m == m)
        ++b;
      if (single_precision_)
        transform(dirs, &order_[a], b - a, freq_float_, scratch_float_);
      else
        transform(dirs, &order_[a], b - a, freq_, scratch_);
      a = b;
    }
  }

  // indexing_api's refine_direction.
  Direction refine_direction(const Direction &candidate, double current_grid,
                             double target_grid) {
    const double incr = current_grid / 4.0;
    if (incr < target_grid)
      return candidate;
    const Direction direction0(candidate.dvec);
    const Direction direction1(candidate.psi + incr, candidate.phi);
    const Vec3 rotaxisi = normalize(direction0.dvec.cross(direction1.dvec));
    const Vec3 rotaxisj = normalize(rotaxisi.cross(direction1.dvec));
    std::vector<Direction> trials;
    trials.reserve(25);
    for (int i = -2; i < 3; ++i) {
      const Vec3 rotate1 = unit_rotate(direction0.dvec, rotaxisi, -i * incr);
      for (int j = -2; j < 3; ++j) {
        const Vec3 rotate2 = unit_rotate(rotate1, rotaxisj, -j * incr);
        if (!(std::abs(1.0 - length(rotate2)) < 0.00001))
          throw std::logic_error("refine_direction: rotation not unit");
        trials.emplace_back(rotate2);
      }
    }
    evaluate_all(trials);
    Direction best = candidate;
    for (const Direction &trial : trials)
      if (trial.kval > best.kval)
        best = trial;
    return refine_direction(best, incr, target_grid);
  }

private:
  // The projections' range and the transform length for one direction.
  struct Extent {
    double pmin = 0, delta_p = 0;
    std::size_t m = 0;
  };

  void project(const Direction &d) {
    const std::size_t n = xyz_.size();
    projections_.resize(n);
    for (std::size_t k = 0; k < n; ++k)
      projections_[k] = xyz_[k][0] * d.dvec[0] + xyz_[k][1] * d.dvec[1] +
                        xyz_[k][2] * d.dvec[2];
  }

  Extent extent(const Direction &d) {
    project(d);
    const double pmin =
        *std::min_element(projections_.begin(), projections_.end());
    const double pmax =
        *std::max_element(projections_.begin(), projections_.end());
    std::size_t m =
        1 + static_cast<std::size_t>((pmax - pmin) * kGranularity * amax_);
    m = primecheck(m);
    m = static_cast<std::size_t>(std::max<int>(2, static_cast<int>(m)));
    const double delta_p = (pmax - pmin) / static_cast<double>(m - 1);
    return {pmin, delta_p, m};
  }

  // The histograms of up to kLanes directions of one transform length, one
  // per lane, transformed together; each direction's peak is then found in
  // its own lane. Unused lanes transform zeros.
  template <class L>
  void transform(std::vector<Direction> &dirs, const std::size_t *index,
                 std::size_t count, std::vector<L> &freq,
                 std::vector<L> &scratch) {
    const std::size_t m = extents_[index[0]].m;
    const fftpack::RealToComplex &fft = plan(m);
    freq.assign(fft.m_real(), L{});
    for (std::size_t l = 0; l < count; ++l) {
      const Extent &e = extents_[index[l]];
      project(dirs[index[l]]);
      for (double p : projections_)
        freq[static_cast<std::size_t>(
            static_cast<int>(((p - e.pmin) / e.delta_p) + 0.5))][l] += 1;
    }
    scratch.resize(std::max(scratch.size(), m));
    fft.forward(freq.data(), scratch.data());
    // The peak of each lane, by rstbx's rule: every lane's squared magnitudes
    // are screened at once against that lane's running bound, and only a
    // coefficient that passes is compared by magnitude.
    constexpr std::size_t kMax = fftpack::kFloatLanes;
    static_assert(fftpack::kLanes <= kMax);
    const std::size_t n_complex = fft.n_complex();
    std::size_t p_kmax[kMax] = {};
    double p_kval[kMax] = {}, p_square[kMax] = {}, bound[kMax] = {};
    for (std::size_t k = kF0Cutoff; k < n_complex; ++k) {
      const L &re = freq[2 * k], &im = freq[2 * k + 1];
      double square[kMax];
      bool any = false;
      for (std::size_t l = 0; l < count; ++l) {
        const double r = re[l], i = im[l];
        square[l] = r * r + i * i;
        any |= square[l] >= bound[l];
      }
      if (!any)
        continue;
      for (std::size_t l = 0; l < count; ++l) {
        if (!(square[l] >= bound[l]))
          continue;
        const double value = std::abs(std::complex<double>(re[l], im[l]));
        if (value > p_kval[l]) {
          p_kval[l] = value;
          p_square[l] = square[l];
          p_kmax[l] = k;
          bound[l] = p_square[l] * (1.0 - kKvalTieMargin);
        }
      }
    }
    for (std::size_t l = 0; l < count; ++l) {
      Direction &d = dirs[index[l]];
      const Extent &e = extents_[index[l]];
      d.kmax = p_kmax[l];
      d.kval = p_kval[l];
      d.m = m;
      d.delta_p = e.delta_p;
      d.pmin = e.pmin;
      d.uc_length =
          static_cast<double>(p_kmax[l]) / (static_cast<double>(m) * e.delta_p);
    }
  }

  const fftpack::RealToComplex &plan(std::size_t m) {
    auto it = plans_.find(m);
    if (it == plans_.end())
      it = plans_.emplace(m, fftpack::RealToComplex(m)).first;
    return it->second;
  }

  const std::vector<Vec3> &xyz_;
  double amax_;
  std::map<std::size_t, fftpack::RealToComplex> plans_;
  std::vector<double> projections_;
  std::vector<Extent> extents_;
  std::vector<std::size_t> order_;
  bool single_precision_;
  std::vector<fftpack::Lanes> freq_, scratch_;
  std::vector<fftpack::FloatLanes> freq_float_, scratch_float_;
};

std::vector<Direction> construct_hemisphere_grid(double sampling) {
  const int psi_index_range = static_cast<int>(0.5 + kPi_2 / sampling);
  const double adjusted_psi_incr = kPi_2 / psi_index_range;
  std::vector<Direction> angles;
  angles.reserve(
      static_cast<std::size_t>(4 * psi_index_range * psi_index_range));
  for (int x = 0; x <= psi_index_range; ++x) {
    double psi = x * adjusted_psi_incr;
    if (psi > kPi) {
      const double eps = 1E-4;
      psi = kPi - eps;
    }
    if (psi == 0) {
      angles.emplace_back(psi, 0.);
    } else {
      const int phi_index_range =
          static_cast<int>(0.5 + 2. * kPi * std::sin(psi) / sampling);
      const double adjusted_phi_incr = 2. * kPi / phi_index_range;
      for (int y = 0; y < phi_index_range; ++y)
        angles.emplace_back(psi, y * adjusted_phi_incr);
    }
  }
  return angles;
}

// flex.sort_permutation(kvals, reverse=True): stable, descending.
std::vector<std::size_t> sort_by_kval(const std::vector<Direction> &d) {
  std::vector<std::size_t> perm(d.size());
  for (std::size_t i = 0; i < perm.size(); ++i)
    perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
    return d[a].kval > d[b].kval;
  });
  return perm;
}

// SimpleSampler.get_top_solutions: the directions above cutoff, the strongest
// size of them if there are three or more.
std::vector<Direction> simple_top_solutions(DpsCore &ai,
                                            const std::vector<Direction> &dirs,
                                            std::size_t size,
                                            double cutoff_divisor) {
  const double kval_cutoff = static_cast<double>(ai.size()) / cutoff_divisor;
  std::vector<Direction> evaluated = dirs;
  ai.evaluate_all(evaluated);
  std::vector<Direction> kept;
  for (const Direction &d : evaluated)
    if (d.kval > kval_cutoff)
      kept.push_back(d);
  if (kept.size() < 3)
    return kept;
  const std::vector<std::size_t> perm = sort_by_kval(kept);
  std::vector<Direction> out;
  for (std::size_t i = 0; i < std::min(size, perm.size()); ++i)
    out.push_back(kept[perm[i]]);
  return out;
}

// HemisphereSampler.get_top_solutions: the directions shorter than max_cell
// and above cutoff, then the strongest of each 0.087-chord cluster, up to
// size.
std::vector<Direction>
hemisphere_top_solutions(DpsCore &ai, const std::vector<Direction> &dirs,
                         std::size_t size, double cutoff_divisor,
                         double max_cell) {
  const double kval_cutoff = static_cast<double>(ai.size()) / cutoff_divisor;
  std::vector<Direction> evaluated = dirs;
  ai.evaluate_all(evaluated);
  std::vector<Direction> kept;
  for (const Direction &d : evaluated)
    if (d.uc_length < max_cell && d.kval > kval_cutoff)
      kept.push_back(d);
  if (kept.size() < 3)
    return kept;
  const std::vector<std::size_t> perm = sort_by_kval(kept);
  std::vector<Direction> out;
  std::size_t perm_idx = 0, unique_clusters = 0;
  while (perm_idx < perm.size() && unique_clusters < size) {
    const Direction &test = kept[perm[perm_idx]];
    bool ok = true;
    for (const Direction &item : out) {
      const double distance = std::sqrt(py_pow(item.dvec[0] - test.dvec[0], 2) +
                                        py_pow(item.dvec[1] - test.dvec[1], 2) +
                                        py_pow(item.dvec[2] - test.dvec[2], 2));
      if (distance < 0.087) {
        ok = false;
        break;
      }
    }
    if (ok) {
      unique_clusters++;
      out.push_back(test);
    }
    perm_idx++;
  }
  return out;
}

// restrict_target_angles_to_quick_winners: the grid directions among the
// kn nearest (squared distance summed over x, y, z, as ANN sums it) of a quick
// winner, or of its opposite for a winner on the equator, that lie closer than
// quick_grid. cctbx builds ANN without self-matches, so a grid direction at
// distance zero from the query is never one of the nearest.
//
// Only directions closer than quick_grid can be accepted, and among them the
// kn nearest are the ones the kn nearest overall would contribute, so only
// those need finding: the grid is hashed into cubes of side quick_grid and a
// query looks in the 27 around it. ANN's exact kd-tree search can differ from
// this only where two candidates tie exactly at the kn-th place.
std::vector<Direction>
restrict_to_quick_winners(const std::vector<Direction> &angles,
                          const std::vector<Direction> &quick_sampling,
                          double second_round_sampling, double quick_grid) {
  const std::size_t kn =
      static_cast<std::size_t>(static_cast<int>(2. * second_round_sampling));
  std::vector<Vec3> queries;
  for (const Direction &q : quick_sampling) {
    queries.push_back(q.dvec);
    if (std::fabs(kPi_2 - q.psi) < 0.0001)
      queries.push_back(Vec3(-q.dvec[0], -q.dvec[1], -q.dvec[2]));
  }
  auto cell = [quick_grid](double x) {
    return static_cast<long>(std::floor((x + 1.0) / quick_grid));
  };
  auto key = [](long ix, long iy, long iz) {
    return (ix * 4096 + iy) * 4096 + iz;
  };
  std::map<long, std::vector<std::size_t>> buckets;
  for (std::size_t i = 0; i < angles.size(); ++i) {
    const Vec3 &p = angles[i].dvec;
    buckets[key(cell(p[0]), cell(p[1]), cell(p[2]))].push_back(i);
  }
  std::vector<bool> accept(angles.size(), false);
  std::vector<std::pair<double, std::size_t>> near;
  for (const Vec3 &q : queries) {
    near.clear();
    const long cx = cell(q[0]), cy = cell(q[1]), cz = cell(q[2]);
    for (long dx = -1; dx <= 1; ++dx)
      for (long dy = -1; dy <= 1; ++dy)
        for (long dz = -1; dz <= 1; ++dz) {
          const auto it = buckets.find(key(cx + dx, cy + dy, cz + dz));
          if (it == buckets.end())
            continue;
          for (std::size_t i : it->second) {
            const Vec3 &p = angles[i].dvec;
            double dist = 0;
            for (int c = 0; c < 3; ++c) {
              const double t = q[c] - p[c];
              dist = dist + t * t;
            }
            if (dist != 0 && std::sqrt(dist) < quick_grid)
              near.emplace_back(dist, i);
          }
        }
    const std::size_t k = std::min(kn, near.size());
    std::partial_sort(near.begin(),
                      near.begin() + static_cast<std::ptrdiff_t>(k),
                      near.end());
    for (std::size_t r = 0; r < k; ++r)
      accept[near[r].second] = true;
  }
  std::vector<Direction> restricted;
  for (std::size_t i = 0; i < angles.size(); ++i)
    if (accept[i])
      restricted.push_back(angles[i]);
  return restricted;
}

// dps_core::setSolutions: rank by kval with libstdc++'s std::sort (not
// stable), zero any direction within 10 degrees of a stronger one, rank again
// and drop the zeros.
std::vector<Direction> set_solutions(std::vector<Direction> s) {
  std::sort(s.begin(), s.end(), KvalCmp());
  if (s.empty())
    return s;
  for (std::size_t i = s.size() - 1; i > 0; --i)
    for (std::size_t j = 0; j < i; ++j)
      if (is_nearly_collinear(s[i], s[j])) {
        s[i].kval = 0.0;
        break;
      }
  std::sort(s.begin(), s.end(), KvalCmp());
  while (s[s.size() - 1].kval == 0.0)
    s.pop_back();
  return s;
}

std::array<double, 5> trace_row(const Direction &d) {
  return {d.kval, d.uc_length, d.dvec[0], d.dvec[1], d.dvec[2]};
}

} // namespace

std::vector<Vec3> hemisphere_grid(double sampling) {
  std::vector<Vec3> out;
  for (const Direction &d : construct_hemisphere_grid(sampling))
    out.push_back(d.dvec);
  return out;
}

std::vector<Vec3> fft1d_basis_vectors(const std::vector<Vec3> &rlp,
                                      double max_cell,
                                      std::optional<double> characteristic_grid,
                                      Fft1dTrace *trace,
                                      const Fft1dGrid &settings) {
  if (!(rlp.size() > 7))
    throw TooFewSpots("fft1d needs more than 7 reciprocal lattice points");

  // lattice.py's index(): the grid from the highest-resolution spot.
  double cg;
  if (characteristic_grid) {
    cg = *characteristic_grid;
  } else {
    const double rossmann_suggestion = settings.max_grid;
    double max_norm = length(rlp[0]);
    for (const Vec3 &v : rlp)
      max_norm = std::max(max_norm, length(v));
    const double naive_obs_highest_resolution = 1. / max_norm;
    const double grid = naive_obs_highest_resolution / max_cell;
    const double new_suggestion = 2. * grid;
    cg = std::min(rossmann_suggestion, new_suggestion);
  }
  DpsCore ai(rlp, max_cell, settings.single_precision);

  // hemisphere_shortcut: HemisphereSampler(max_grid=0.029,
  // characteristic_grid=cg, quick_grid=0.016).hemisphere(ai, max_cell, 30, 4).
  using Clock = std::chrono::steady_clock;
  const auto seconds = [](Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  const auto t0 = Clock::now();
  const double max_grid = settings.max_grid, quick_grid = settings.quick_grid;
  const double incr = std::min(max_grid, cg);
  std::vector<Direction> angles;
  bool adaptive = false;
  std::size_t n_quick = 0, n_quick_winners = 0;
  if (quick_grid >= incr) {
    adaptive = true;
    const double estimate_number_of_quick_grids =
        2. * kPi / (quick_grid * quick_grid);
    const double second_round_sampling =
        (kPi * quick_grid * quick_grid) / (incr * incr);
    const std::size_t target_size = static_cast<std::size_t>(static_cast<int>(
        estimate_number_of_quick_grids / 2. / second_round_sampling));
    const std::vector<Direction> quick_directions =
        construct_hemisphere_grid(quick_grid);
    n_quick = quick_directions.size();
    const std::vector<Direction> quick_sampling =
        simple_top_solutions(ai, quick_directions, target_size, 10.);
    n_quick_winners = quick_sampling.size();
    angles = restrict_to_quick_winners(construct_hemisphere_grid(incr),
                                       quick_sampling, second_round_sampling,
                                       quick_grid);
  } else {
    angles = construct_hemisphere_grid(incr);
  }
  const auto t1 = Clock::now();
  const std::vector<Direction> unrefined =
      hemisphere_top_solutions(ai, angles, 30, 4., max_cell);
  const auto t2 = Clock::now();
  std::vector<Direction> refined;
  const double final_target_grid = cg / settings.refine_divisor;
  for (const Direction &item : unrefined)
    refined.push_back(ai.refine_direction(item, incr, final_target_grid));
  const auto t3 = Clock::now();
  if (trace) {
    trace->seconds_quick = seconds(t0, t1);
    trace->seconds_scan = seconds(t1, t2);
    trace->seconds_refine = seconds(t2, t3);
    trace->characteristic_grid = cg;
    trace->incr = incr;
    trace->adaptive = adaptive;
    trace->n_quick = n_quick;
    trace->n_quick_winners = n_quick_winners;
    trace->n_directions = angles.size();
    trace->unrefined.clear();
    for (const Direction &d : unrefined)
      trace->unrefined.push_back(trace_row(d));
    trace->refined.clear();
    for (const Direction &d : refined)
      trace->refined.push_back(trace_row(d));
  }
  const std::vector<Direction> solutions = set_solutions(std::move(refined));
  std::vector<Vec3> out;
  out.reserve(solutions.size());
  for (const Direction &s : solutions)
    out.push_back(Vec3(s.uc_length * s.dvec[0], s.uc_length * s.dvec[1],
                       s.uc_length * s.dvec[2]));
  return out;
}

} // namespace mxi::stills

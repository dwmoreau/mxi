#include "nave.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numeric>

namespace mxi::stills {

namespace {

constexpr double kPi = 3.141592653589793;
// The likelihood's switch between inside and outside the band: 10 is soft.
constexpr double kEscalate = 10.0;

// Ends the likelihood evaluation with a failure.
struct ObjectiveFailure {
  NaveStatus status;
};

// flex.mean: the sum times the reciprocal of the count.
double mean_of(const std::vector<double> &v) {
  double s = v[0];
  for (std::size_t i = 1; i < v.size(); ++i)
    s += v[i];
  return s * (1. / static_cast<double>(v.size()));
}

double sum_of(const std::vector<double> &v) {
  double s = 0.0;
  for (double x : v)
    s += x;
  return s;
}

// serialtbx max_like.minimizer: the negative log-likelihood of the
// Delta-Psi values and its gradient in (log alpha, log eta).
class Likelihood {
public:
  Likelihood(const std::vector<double> &d, const std::vector<double> &psi)
      : d_(d), psi_(psi), safelog_(-1. + std::log(DBL_MAX)) {}

  double operator()(const std::vector<double> &x,
                    std::vector<double> &g) const {
    if (x[0] > safelog_ || x[1] > safelog_)
      throw ObjectiveFailure{NaveStatus::outside_domain};
    const double alpha = std::exp(x[0]);
    const double eta = std::exp(x[1]);
    const double f = functional(alpha, eta);
    double g0 = 0., g1 = 0.;
    const double expnu = 1. + std::exp(kEscalate);
    for (std::size_t i = 0; i < d_.size(); ++i) {
      const double psi_model = (d_[i] * alpha + eta) / 2.;
      const double pm0 = 0.5 * d_[i] * alpha;
      const double pm1 = 0.5 * eta;
      const double psi_i = psi_[i];
      const double B = kEscalate / psi_model;
      if (psi_model > 1e100 || pm0 > 1e100)
        throw ObjectiveFailure{NaveStatus::diverged};
      const double expB = std::exp(B * (psi_i + psi_model));
      const double expBneg = std::exp(-B * (psi_i - psi_model));
      if (expB > 1e100 || expBneg < -1e100 || expBneg > 1e100)
        throw ObjectiveFailure{NaveStatus::diverged};
      const double Spos = 1. + expB;
      const double Sneg = 1. + expBneg;
      const double pB0 = (-kEscalate / (psi_model * psi_model)) * pm0;
      const double pB1 = (-kEscalate / (psi_model * psi_model)) * pm1;
      const double pSpos0 = expB * ((psi_i + psi_model) * pB0 + B * pm0);
      const double pSpos1 = expB * ((psi_i + psi_model) * pB1 + B * pm1);
      const double pSneg0 = expBneg * ((-psi_i + psi_model) * pB0 + B * pm0);
      const double pSneg1 = expBneg * ((-psi_i + psi_model) * pB1 + B * pm1);
      const double pG0 = -expnu * std::pow(Sneg, -2.) * pSneg0;
      const double pG1 = -expnu * std::pow(Sneg, -2.) * pSneg1;
      const double Sfac = 2. * psi_model * Spos;
      const double pF0 =
          -expnu * std::pow(Sfac, -2.) * 2 * (psi_model * pSpos0 + Spos * pm0);
      const double pF1 =
          -expnu * std::pow(Sfac, -2.) * 2 * (psi_model * pSpos1 + Spos * pm1);
      const double fx = (0.5 / psi_model) / (Spos)*expnu;
      const double gx = (1. / Sneg) * expnu;
      const double prob = fx * gx;
      const double pprob0 = fx * pG0 + gx * pF0;
      const double pprob1 = fx * pG1 + gx * pF1;
      g0 -= (1. / prob) * pprob0;
      g1 -= (1. / prob) * pprob1;
    }
    g[0] = g0;
    g[1] = g1;
    return f;
  }

private:
  double functional(double alpha, double eta) const {
    const double expnu = 1 + std::exp(kEscalate);
    double f = 0.;
    for (std::size_t i = 0; i < d_.size(); ++i) {
      const double psi_model = (d_[i] * alpha + eta) / 2.;
      const double psi_i = psi_[i];
      const double B = kEscalate / psi_model;
      const double arg = B * (psi_i + psi_model);
      const double negarg = -B * (psi_i - psi_model);
      if (std::abs(arg) > safelog_ || std::abs(negarg) > safelog_)
        throw ObjectiveFailure{NaveStatus::outside_domain};
      const double fx = (0.5 / psi_model) / (1 + std::exp(arg)) * expnu;
      const double gx = 1. / (1 + std::exp(negarg)) * expnu;
      const double prob = fx * gx;
      if (prob <= 0.)
        throw ObjectiveFailure{NaveStatus::zero_probability};
      f -= std::log(prob);
    }
    return f;
  }

  const std::vector<double> &d_, &psi_;
  const double safelog_;
};

// NaveParameters.ewald_proximal_volume: half the volume of reciprocal space
// within the model's band, out to the resolution cutoff.
double ewald_proximal_volume(double wavelength, double cutoff,
                             double domain_size_ang, double full_mosaicity) {
  const double R_L = 1.0 / wavelength;
  const double TT = 2.0 * std::asin(wavelength / (2.0 * cutoff));
  const double part_vol = kPi * (2.0 / 3.0) * (1.0 - std::cos(TT));
  const double ewald_sphere_volume = part_vol * std::pow(R_L, 3.0);
  const double R_prime = R_L + 1.0 / domain_size_ang;
  const double domain_size_volume = part_vol * std::pow(R_prime, 3.0);
  constexpr int kTerms = 100;
  double summation = 0.0;
  for (int x = 0; x < kTerms; ++x) {
    const double phi = (static_cast<double>(x) / kTerms) * TT;
    const double integral =
        std::pow(R_prime + (full_mosaicity * R_L * std::sin(phi) / 2.0), 3.0) -
        std::pow(R_prime, 3.0);
    summation += (integral * std::sin(phi)) * (TT / kTerms);
  }
  const double mosaicity_volume = (2.0 / 3.0) * kPi * summation;
  return (domain_size_volume - ewald_sphere_volume) + mosaicity_volume;
}

} // namespace

const char *to_string(NaveStatus status) {
  switch (status) {
  case NaveStatus::ok:
    return "ok";
  case NaveStatus::too_few:
    return "too_few";
  case NaveStatus::singular_seed:
    return "singular_seed";
  case NaveStatus::zero_seed_eta:
    return "zero_seed_eta";
  case NaveStatus::outside_domain:
    return "outside_domain";
  case NaveStatus::diverged:
    return "diverged";
  case NaveStatus::zero_probability:
    return "zero_probability";
  case NaveStatus::fit_failed:
    return "fit_failed";
  case NaveStatus::ewald_volume:
    return "ewald_volume";
  }
  return "unknown";
}

bool failed(NaveStatus status) { return status != NaveStatus::ok; }

NaveResult nave_parameters(const std::vector<double> &delpsi,
                           const std::vector<double> &d,
                           const std::vector<double> &two_theta,
                           double cell_volume, double wavelength,
                           const NaveOptions &options) {
  NaveResult out;
  const std::size_t n = delpsi.size();
  const std::size_t n_bins =
      std::min<std::size_t>(std::max<std::size_t>(3, n / 25), 50);
  const std::size_t bin_size = n / n_bins;
  if (bin_size == 0) {
    out.status = NaveStatus::too_few;
    return out;
  }

  // The seed.
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t a, std::size_t b) {
                     return two_theta[a] < two_theta[b];
                   });
  std::vector<double> d_env(n_bins), psi_env(n_bins), subset(bin_size);
  for (std::size_t b = 0; b < n_bins; ++b) {
    for (std::size_t k = 0; k < bin_size; ++k)
      subset[k] = d[order[b * bin_size + k]];
    d_env[b] = mean_of(subset);
    double largest = std::abs(delpsi[order[b * bin_size]]);
    for (std::size_t k = 1; k < bin_size; ++k) {
      const double v = std::abs(delpsi[order[b * bin_size + k]]);
      if (largest < v)
        largest = v;
    }
    psi_env[b] = largest;
  }
  std::vector<double> dd(n_bins), dpsi(n_bins);
  for (std::size_t b = 0; b < n_bins; ++b) {
    dd[b] = d_env[b] * d_env[b];
    dpsi[b] = d_env[b] * psi_env[b];
  }
  const double s_dd = sum_of(dd), s_d = sum_of(d_env), s_dpsi = sum_of(dpsi),
               s_psi = sum_of(psi_env);
  // scitbx.matrix: the 2x2 inverse as the transposed cofactors over the
  // determinant, then the product with the right-hand side.
  const double count = static_cast<double>(n_bins);
  const double det = s_dd * count - s_d * s_d;
  if (det == 0.0) {
    out.status = NaveStatus::singular_seed;
    return out;
  }
  const double i00 = count / det, i01 = -s_d / det, i10 = -s_d / det,
               i11 = s_dd / det;
  const double slope = i00 * s_dpsi + i01 * s_psi;
  const double intercept = i10 * s_dpsi + i11 * s_psi;
  if (slope == 0.0) {
    out.status = NaveStatus::singular_seed;
    return out;
  }
  const double scherrer = 1.0 / (2 * slope);
  const double lower_limit = std::pow(cell_volume, 1.0 / 3.0) * 3;
  out.domain_size_seed = lower_limit > scherrer ? lower_limit : scherrer;
  out.eta_seed = std::abs(2.0 * intercept);
  if (out.eta_seed == 0.0) {
    out.status = NaveStatus::zero_seed_eta;
    return out;
  }

  // The maximum-likelihood fit.
  std::vector<double> x{std::log(2. / out.domain_size_seed),
                        std::log(out.eta_seed)};
  const Likelihood likelihood(d, delpsi);
  lbfgs::CoreParameters core;
  core.gtol = 0.1;
  lbfgs::TerminationParameters termination;
  termination.traditional_convergence_test = false;
  termination.drop_convergence_test_max_drop_eps = 1.e-5;
  termination.min_iterations = 0;
  termination.max_iterations = 100;
  termination.max_calls = 200;
  lbfgs::ExceptionHandlingParameters handling;
  handling.ignore_line_search_failed_step_at_lower_bound = true;
  try {
    out.minimiser = lbfgs::run(x, likelihood, core, termination, handling);
  } catch (const ObjectiveFailure &e) {
    out.status = e.status;
    return out;
  }
  if (out.minimiser.raised) {
    out.status = NaveStatus::fit_failed;
    return out;
  }
  out.alpha = std::exp(x[0]);
  out.eta = std::exp(x[1]);

  // The acceptance band, in degrees.
  const double outer = 0.5 * out.eta * 180.0 / kPi;
  out.accepted.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double band = d[i] / (2.0 / out.alpha) * 180.0 / kPi + outer;
    out.accepted[i] = std::abs(delpsi[i] * 180.0 / kPi) < band;
  }

  const double ml_domain_size = 2.0 / out.alpha;
  out.half_mosaicity_deg =
      out.eta * 180.0 / (2.0 * kPi) * options.model_expansion_factor;
  out.domain_size_ang = ml_domain_size / options.model_expansion_factor;
  out.ewald_volume = ewald_proximal_volume(
      wavelength, options.ewald_proximity_resolution_cutoff, ml_domain_size,
      out.eta);
  if (out.ewald_volume > options.ewald_proximal_volume_max)
    out.status = NaveStatus::ewald_volume;
  return out;
}

std::vector<double> cell_d_spacings(const Mat3 &B,
                                    const std::vector<std::array<int, 3>> &h) {
  std::vector<double> d(h.size());
  for (std::size_t i = 0; i < h.size(); ++i) {
    const Vec3 q =
        B * Vec3{static_cast<double>(h[i][0]), static_cast<double>(h[i][1]),
                 static_cast<double>(h[i][2])};
    d[i] = 1.0 / std::sqrt(q.dot(q));
  }
  return d;
}

std::vector<double> two_theta_deg(const std::vector<double> &d,
                                  double wavelength) {
  std::vector<double> t(d.size());
  for (std::size_t i = 0; i < d.size(); ++i)
    t[i] = 2.0 * std::asin(wavelength / (2.0 * d[i])) * 180.0 / kPi;
  return t;
}

double cell_volume(const Mat3 &B) { return 1.0 / B.determinant(); }

} // namespace mxi::stills

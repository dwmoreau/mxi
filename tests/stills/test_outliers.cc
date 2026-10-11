// The L-BFGS minimiser against trajectories recorded from scitbx, and
// Sauter-Poon outlier rejection and the Nave fit against dials' own calls on
// real stills.

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "check.hh"
#include "real_outliers.hh"
#include "stills/lbfgs.hh"
#include "stills/nave.hh"
#include "stills/outliers.hh"

using namespace mxi::stills;

namespace {

// The extended Rosenbrock function, as the reference script evaluated it.
double rosenbrock(const std::vector<double> &x, std::vector<double> &g) {
  const std::size_t n = x.size();
  double f = 0.0;
  std::fill(g.begin(), g.end(), 0.0);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const double a = x[i + 1] - x[i] * x[i];
    const double b = 1.0 - x[i];
    f += 100.0 * a * a + b * b;
    g[i] += -400.0 * x[i] * a - 2.0 * b;
    g[i + 1] += 200.0 * a;
  }
  return f;
}

void same_bits(const std::vector<double> &a, const std::vector<double> &b,
               const std::string &what) {
  check::equal(static_cast<long long>(a.size()),
               static_cast<long long>(b.size()), what + ": size");
  for (std::size_t i = 0; i < a.size(); ++i)
    check::is_true(a[i] == b[i], what + "[" + std::to_string(i) +
                                     "]: " + std::to_string(a[i]) + " vs " +
                                     std::to_string(b[i]));
}

void same_flags(const std::vector<bool> &a, const std::vector<bool> &b,
                const std::string &what) {
  check::equal(static_cast<long long>(a.size()),
               static_cast<long long>(b.size()), what + ": size");
  for (std::size_t i = 0; i < a.size(); ++i)
    check::is_true(a[i] == b[i], what + " flag " + std::to_string(i));
}

} // namespace

// Final points, iteration and call counts of scitbx.lbfgs.run on the same
// function (claude/rewrite/runs/r7-outliers-nave/PROVENANCE.txt).
TEST(stills_lbfgs_matches_scitbx) {
  struct Case {
    std::vector<double> x0, x1;
    std::size_t iter, nfun;
    bool drop;
    std::optional<std::size_t> max_iterations;
  };
  const std::vector<Case> cases = {
      {{-1.2, 1.0},
       {0x1.000000114f0ecp+0, 0x1.0000002128eb9p+0},
       39,
       48,
       false,
       {}},
      {{-1.2, 1.0, -1.2, 1.0, -1.2, 1.0, 0.5},
       {0x1.0000002835fb7p+0, 0x1.ffffffe133c02p-1, 0x1.ffffffa4c9badp-1,
        0x1.ffffffc3ff793p-1, 0x1.ffffff595a3d6p-1, 0x1.ffffff06cb8cep-1,
        0x1.fffffe0c6b59bp-1},
       66,
       82,
       false,
       {}},
      {{-1.2, 1.0, -1.2, 1.0, -1.2, 1.0, 0.5},
       {0x1.26e16071c2333p-1, 0x1.46b4574bb8f44p-2, 0x1.35c5dcd04ec13p-4,
        0x1.b98ddfb3e6a66p-6, 0x1.98ee1b8e9ef42p-8, 0x1.d53de6170b708p-9,
        0x1.9137b555b73f2p-8},
       8,
       21,
       true,
       100},
      {{-1.2, 1.0},
       {-0x1.05797f72f6e66p+0, 0x1.0a0f326ff4d90p+0},
       5,
       7,
       false,
       5},
  };
  for (const Case &c : cases) {
    std::vector<double> x = c.x0;
    lbfgs::CoreParameters core;
    lbfgs::TerminationParameters term;
    lbfgs::ExceptionHandlingParameters handling;
    term.max_iterations = c.max_iterations;
    if (c.drop) {
      core.gtol = 0.1;
      term.traditional_convergence_test = false;
      term.max_calls = 200;
      handling.ignore_line_search_failed_step_at_lower_bound = true;
    }
    const lbfgs::Result r = lbfgs::run(x, rosenbrock, core, term, handling);
    check::is_true(!r.raised, "no error");
    same_bits(x, c.x1, "x");
    check::equal(static_cast<long long>(r.iter), static_cast<long long>(c.iter),
                 "iterations");
    check::equal(static_cast<long long>(r.nfun), static_cast<long long>(c.nfun),
                 "function calls");
  }
}

TEST(stills_drop_test_on_a_flat_objective) {
  lbfgs::DropConvergenceTest test(5, 1e-5, 2.0);
  for (int i = 0; i < 4; ++i)
    check::is_true(!test(3.0), "fewer than five points never converge");
  check::is_true(test(3.0), "a constant objective has converged");
}

TEST(stills_sauter_poon_matches_dials) {
  for (const auto &c : real_outliers::sauter_poon_cases()) {
    const SauterPoonResult r =
        sauter_poon(c.x_obs, c.y_obs, c.x_calc, c.y_calc);
    check::is_true(std::string(to_string(r.status)) == c.outcome,
                   c.name + ": outcome " + to_string(r.status));
    if (c.outcome != "ok")
      continue;
    check::is_true(r.sigma == c.sigma, c.name + ": sigma");
    same_flags(r.kept, c.kept, c.name);
  }
}

TEST(stills_sauter_poon_refuses_equal_radii) {
  std::vector<double> x_obs(20), y_obs(20), x_calc(20, 0.0), y_calc(20, 0.0);
  for (std::size_t i = 0; i < 20; ++i)
    x_obs[i] = 0.01 * static_cast<double>(i + 1);
  // Spot 3 moved to spot 7's distance, along the other axis.
  x_obs[3] = 0.0;
  y_obs[3] = x_obs[7];
  const SauterPoonResult r = sauter_poon(x_obs, y_obs, x_calc, y_calc);
  check::is_true(r.status == SauterPoonStatus::equal_radii,
                 std::string("outcome ") + to_string(r.status));
}

TEST(stills_sauter_poon_rejects_a_far_spot) {
  // Radii spread as a Rayleigh sample of sigma 0.05 mm, plus one 1 mm off.
  const std::size_t n = 60;
  std::vector<double> x_obs(n), y_obs(n, 0.0), x_calc(n, 0.0), y_calc(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) {
    const double p = (static_cast<double>(k) + 0.5) / static_cast<double>(n);
    x_obs[k] = 0.05 * std::sqrt(-2.0 * std::log(1.0 - p));
  }
  x_obs[n - 1] = 1.0;
  const SauterPoonResult r = sauter_poon(x_obs, y_obs, x_calc, y_calc);
  check::is_true(r.status == SauterPoonStatus::ok, "outcome");
  check::close(r.sigma, 0.05, 0.005, "sigma");
  check::is_true(!r.kept[n - 1], "the far spot is rejected");
  check::is_true(r.kept[0] && r.kept[n / 2], "close spots are kept");
}

TEST(stills_nave_matches_dials) {
  for (const auto &c : real_outliers::nave_cases()) {
    const NaveResult r =
        nave_parameters(c.delpsi, c.d, c.two_theta, c.volume, c.wavelength);
    check::is_true(std::string(to_string(r.status)) == c.outcome,
                   c.name + ": outcome " + to_string(r.status));
    check::is_true(r.eta_seed == c.eta_seed, c.name + ": seed eta");
    check::is_true(r.domain_size_seed == c.domain_size_seed,
                   c.name + ": seed domain size");
    if (c.accepted.empty())
      continue;
    check::is_true(r.alpha == c.alpha, c.name + ": alpha");
    check::is_true(r.eta == c.eta, c.name + ": eta");
    same_flags(r.accepted, c.accepted, c.name);
    check::is_true(r.half_mosaicity_deg == c.half_mosaicity_deg,
                   c.name + ": half-mosaicity");
    check::is_true(r.domain_size_ang == c.domain_size_ang,
                   c.name + ": domain size");
    check::is_true(r.ewald_volume == c.ewald_volume, c.name + ": Ewald volume");
  }
}

TEST(stills_nave_needs_three_reflections) {
  const NaveResult r =
      nave_parameters({1e-3, -2e-3}, {3.0, 2.5}, {10.0, 12.0}, 1e5, 1.0);
  check::is_true(r.status == NaveStatus::too_few, "outcome");
}

TEST(stills_cell_spacings_and_volume) {
  const mxi::Mat3 B{1.0 / 50.0, 0.0, 0.0, 0.0,       1.0 / 60.0,
                    0.0,        0.0, 0.0, 1.0 / 70.0};
  const auto d = cell_d_spacings(B, {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}});
  check::close(d[0], 50.0, 1e-12, "d(100)");
  check::close(d[1], 30.0, 1e-12, "d(020)");
  check::close(d[2], 70.0 / 3.0, 1e-12, "d(003)");
  check::close(cell_volume(B), 50.0 * 60.0 * 70.0, 1e-7, "volume");
  const auto t = two_theta_deg({50.0}, 1.0);
  check::close(t[0], 2.0 * std::asin(0.01) * 180.0 / M_PI, 1e-12, "2-theta");
}

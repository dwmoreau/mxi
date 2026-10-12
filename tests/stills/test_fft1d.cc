// The fft1d basis-vector search: the FFT against a direct transform, the
// hemisphere grid's sizes, and real searches against dials' own results.

#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include "check.hh"
#include "real_fft1d.hh"
#include "stills/fft1d.hh"
#include "stills/fftpack.hh"

using namespace mxi::stills;

TEST(stills_fftpack_matches_a_direct_transform) {
  // Lengths with each factor rfftf handles: 4, 2, 3, 5 and their products.
  for (std::size_t n : {2, 3, 4, 5, 6, 8, 12, 30, 45, 64, 100, 360, 750}) {
    const fftpack::RealToComplex fft(n);
    std::vector<double> x(fft.m_real(), 0.0), scratch(n);
    for (std::size_t j = 0; j < n; ++j)
      x[j] =
          std::sin(0.7 * static_cast<double>(j)) + static_cast<double>(j % 3);
    const std::vector<double> input(x.begin(),
                                    x.begin() + static_cast<long>(n));
    fft.forward(x.data(), scratch.data());
    for (std::size_t k = 0; k < fft.n_complex(); ++k) {
      std::complex<double> sum = 0;
      for (std::size_t j = 0; j < n; ++j)
        sum +=
            input[j] * std::polar(1.0, -2 * M_PI * static_cast<double>(j * k) /
                                           static_cast<double>(n));
      check::close(x[2 * k], sum.real(), 1e-9,
                   "n " + std::to_string(n) + " re");
      check::close(x[2 * k + 1], sum.imag(), 1e-9,
                   "n " + std::to_string(n) + " im");
    }
  }
}

TEST(stills_hemisphere_grid_sizes) {
  // The R2 spec's counts from the built rstbx extension.
  const std::vector<std::pair<double, long long>> sizes = {{0.029, 7557},
                                                           {0.022, 13055},
                                                           {0.016, 24700},
                                                           {0.0116, 46816},
                                                           {0.010, 63109}};
  for (const auto &[s, n] : sizes)
    check::equal(static_cast<long long>(hemisphere_grid(s).size()), n,
                 "grid " + std::to_string(s));
}

TEST(stills_fft1d_matches_dials) {
  for (const auto &c : real_fft1d::cases) {
    std::vector<mxi::Vec3> rlp;
    for (const auto &v : c.rlp)
      rlp.push_back({v[0], v[1], v[2]});
    const auto got = fft1d_basis_vectors(rlp, c.max_cell);
    check::equal(static_cast<long long>(got.size()),
                 static_cast<long long>(c.basis_vectors.size()),
                 c.frame + ": number of vectors");
    for (std::size_t i = 0; i < got.size(); ++i)
      for (std::size_t j = 0; j < 3; ++j)
        check::is_true(got[i][j] == c.basis_vectors[i][j],
                       c.frame + ": vector " + std::to_string(i));
  }
}

TEST(stills_fft1d_refuses_seven_spots) {
  std::vector<mxi::Vec3> rlp;
  for (const auto &v : real_fft1d::refused.rlp)
    rlp.push_back({v[0], v[1], v[2]});
  bool thrown = false;
  try {
    fft1d_basis_vectors(rlp, real_fft1d::refused.max_cell);
  } catch (const TooFewSpots &) {
    thrown = true;
  }
  check::is_true(thrown, real_fft1d::refused.frame + ": refused");
}

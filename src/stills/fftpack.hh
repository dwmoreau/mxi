// The forward real-to-complex FFT as scitbx's fftpack computes it: FFTPACK's
// rfftf (Swarztrauber), factorised and with its twiddle factors computed as
// scitbx does, so that every coefficient equals scitbx's to the last bit. The
// fft1d basis-vector search ranks directions by these magnitudes, and its
// ties are common, so nothing less than that will do.
//
// Derived from cctbx_project (cctbx, rstbx, scitbx): scitbx/fftpack.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#pragma once

#include <cstddef>
#include <type_traits>
#include <vector>

namespace mxi::stills::fftpack {

// As many doubles as a vector register holds: four with AVX, otherwise the two
// of SSE2. Arithmetic on Lanes is lane by lane, each lane exactly as on a
// double.
#if defined(__AVX__)
inline constexpr std::size_t kLanes = 4;
#else
inline constexpr std::size_t kLanes = 2;
#endif
using Lanes = double __attribute__((vector_size(kLanes * sizeof(double))));
// The same register as floats: 2 * kLanes of them.
inline constexpr std::size_t kFloatLanes = 2 * kLanes;
using FloatLanes =
    float __attribute__((vector_size(kFloatLanes * sizeof(float))));

class RealToComplex {
public:
  explicit RealToComplex(std::size_t n_real);

  std::size_t n_real() const { return n_; }
  // n_real / 2 + 1 complex values.
  std::size_t n_complex() const { return n_ / 2 + 1; }
  // The doubles the transform reads and writes: 2 * n_complex().
  std::size_t m_real() const { return 2 * n_complex(); }

  // In place: the first n_real() values of seq are the input; on return the
  // first m_real() hold n_complex() complex values (re, im). scratch holds at
  // least n_real() values. T is double, or Lanes for kLanes transforms at once,
  // one in each lane, each computed exactly as the double transform computes
  // it; float and FloatLanes are the same transform in single precision, with
  // the twiddles rounded to float.
  template <class T> void forward(T *seq, T *scratch) const;

private:
  template <class T> void forward_compressed(T *c, T *ch) const;
  void compute_twiddles();
  template <class W> const W *twiddles() const {
    if constexpr (std::is_same_v<W, float>)
      return wa_float_.data();
    else
      return wa_.data();
  }

  std::size_t n_;
  std::vector<int> factors_;
  std::vector<double> wa_;
  std::vector<float> wa_float_;
};

} // namespace mxi::stills::fftpack

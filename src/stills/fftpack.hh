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
#include <vector>

namespace mxi::stills::fftpack {

class RealToComplex {
public:
  explicit RealToComplex(std::size_t n_real);

  std::size_t n_real() const { return n_; }
  // n_real / 2 + 1 complex values.
  std::size_t n_complex() const { return n_ / 2 + 1; }
  // The doubles the transform reads and writes: 2 * n_complex().
  std::size_t m_real() const { return 2 * n_complex(); }

  // In place: the first n_real() doubles of seq are the input; on return the
  // first m_real() hold n_complex() complex values (re, im). scratch holds at
  // least n_real() doubles.
  void forward(double *seq, double *scratch) const;

private:
  void forward_compressed(double *c, double *ch) const;

  std::size_t n_;
  std::vector<int> factors_;
  std::vector<double> wa_;
};

} // namespace mxi::stills::fftpack

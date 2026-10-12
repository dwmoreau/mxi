// Derived from cctbx_project (cctbx, rstbx, scitbx): scitbx/fftpack.
// cctbx: Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to receipt
// of any required approvals from the U.S. Dept. of Energy). All rights
// reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "stills/fftpack.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <math.h>

namespace mxi::stills::fftpack {

// scitbx's compiled transform takes each sine and cosine of one angle from a
// single sincos call, which can differ from sin and cos in the last bit; so do
// these.

namespace {

// The scalar a value type is made of: double for double and Lanes, float for
// float and FloatLanes.
template <class T> struct Element {
  using type = double;
};
template <> struct Element<float> {
  using type = float;
};
template <> struct Element<FloatLanes> {
  using type = float;
};
template <class T> using E = typename Element<T>::type;

// scitbx's ref_3d_tp: (ix, iy, iz) -> start[(iz * ny + iy) * nx + ix].
template <class T> struct Ref3 {
  T *p;
  std::size_t nx, ny;
  T &operator()(std::size_t ix, std::size_t iy, std::size_t iz) const {
    return p[(iz * ny + iy) * nx + ix];
  }
};

template <class T> struct Ref2 {
  T *p;
  std::size_t nx;
  T &operator()(std::size_t ix, std::size_t iy) const {
    return p[iy * nx + ix];
  }
};

int count_reduce(int &red_n, int factor) {
  int result = 0;
  while (red_n % factor == 0) {
    red_n /= factor;
    result++;
  }
  return result;
}

// scitbx's factorization(n, real_to_complex=true): 4s, then 2s, then 3s,
// then odd factors from 5 up.
std::vector<int> factorize(std::size_t n) {
  const std::array<int, 3> opt_factors = {3, 4, 2};
  const std::array<int, 3> perm = {2, 1, 0};
  std::array<int, 3> count = {0, 0, 0};
  int red_n = static_cast<int>(n);
  for (std::size_t i = 0; red_n > 1 && i < opt_factors.size(); i++)
    count[i] = count_reduce(red_n, opt_factors[i]);
  std::vector<int> factors;
  for (std::size_t i = 0; i < opt_factors.size(); i++)
    factors.insert(factors.end(), static_cast<std::size_t>(count[perm[i]]),
                   opt_factors[perm[i]]);
  for (int factor = 5; red_n > 1; factor += 2)
    factors.insert(factors.end(),
                   static_cast<std::size_t>(count_reduce(red_n, factor)),
                   factor);
  return factors;
}

template <class T>
void passf2(std::size_t ido, std::size_t l1, T *cc_begin, T *ch_begin,
            const E<T> *wa1) {
  const Ref3<T> cc{cc_begin, ido, l1};
  const Ref3<T> ch{ch_begin, ido, 2};
  for (std::size_t k = 0; k < l1; k++) {
    ch(0, 0, k) = cc(0, k, 0) + cc(0, k, 1);
    ch(ido - 1, 1, k) = cc(0, k, 0) - cc(0, k, 1);
  }
  if (ido < 2)
    return;
  if (ido > 2) {
    for (std::size_t k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
        const std::size_t i1 = i0 + 1;
        const std::size_t ic = ido - i1;
        const T tr2 = wa1[i0 - 1] * cc(i0, k, 1) + wa1[i0] * cc(i1, k, 1);
        const T ti2 = wa1[i0 - 1] * cc(i1, k, 1) - wa1[i0] * cc(i0, k, 1);
        ch(i1, 0, k) = cc(i1, k, 0) + ti2;
        ch(ic, 1, k) = ti2 - cc(i1, k, 0);
        ch(i0, 0, k) = cc(i0, k, 0) + tr2;
        ch(ic - 1, 1, k) = cc(i0, k, 0) - tr2;
      }
    }
    if (ido % 2 != 0)
      return;
  }
  for (std::size_t k = 0; k < l1; k++) {
    ch(0, 1, k) = -cc(ido - 1, k, 1);
    ch(ido - 1, 0, k) = cc(ido - 1, k, 0);
  }
}

template <class T>
void passf3(std::size_t ido, std::size_t l1, T *cc_begin, T *ch_begin,
            const E<T> *wa1, const E<T> *wa2) {
  const Ref3<T> cc{cc_begin, ido, l1};
  const Ref3<T> ch{ch_begin, ido, 3};
  const double taur_d = -.5;
  const E<T> taur = static_cast<E<T>>(taur_d);
  const E<T> taui = static_cast<E<T>>(-taur_d * std::sqrt(3.0));
  for (std::size_t k = 0; k < l1; k++) {
    const T cr2 = cc(0, k, 1) + cc(0, k, 2);
    ch(0, 0, k) = cc(0, k, 0) + cr2;
    ch(0, 2, k) = taui * (cc(0, k, 2) - cc(0, k, 1));
    ch(ido - 1, 1, k) = cc(0, k, 0) + taur * cr2;
  }
  if (ido == 1)
    return;
  for (std::size_t k = 0; k < l1; k++) {
    for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
      const std::size_t i1 = i0 + 1;
      const std::size_t ic = ido - i1;
      const T dr2 = wa1[i0 - 1] * cc(i0, k, 1) + wa1[i0] * cc(i1, k, 1);
      const T di2 = wa1[i0 - 1] * cc(i1, k, 1) - wa1[i0] * cc(i0, k, 1);
      const T dr3 = wa2[i0 - 1] * cc(i0, k, 2) + wa2[i0] * cc(i1, k, 2);
      const T di3 = wa2[i0 - 1] * cc(i1, k, 2) - wa2[i0] * cc(i0, k, 2);
      const T cr2 = dr2 + dr3;
      const T ci2 = di2 + di3;
      ch(i0, 0, k) = cc(i0, k, 0) + cr2;
      ch(i1, 0, k) = cc(i1, k, 0) + ci2;
      const T tr2 = cc(i0, k, 0) + taur * cr2;
      const T ti2 = cc(i1, k, 0) + taur * ci2;
      const T tr3 = taui * (di2 - di3);
      const T ti3 = taui * (dr3 - dr2);
      ch(i0, 2, k) = tr2 + tr3;
      ch(ic - 1, 1, k) = tr2 - tr3;
      ch(i1, 2, k) = ti2 + ti3;
      ch(ic, 1, k) = ti3 - ti2;
    }
  }
}

template <class T>
void passf4(std::size_t ido, std::size_t l1, T *cc_begin, T *ch_begin,
            const E<T> *wa1, const E<T> *wa2, const E<T> *wa3) {
  const Ref3<T> cc{cc_begin, ido, l1};
  const Ref3<T> ch{ch_begin, ido, 4};
  const E<T> hsqt2 = static_cast<E<T>>(.5 * std::sqrt(2.0));
  for (std::size_t k = 0; k < l1; k++) {
    const T tr1 = cc(0, k, 1) + cc(0, k, 3);
    const T tr2 = cc(0, k, 0) + cc(0, k, 2);
    ch(0, 0, k) = tr1 + tr2;
    ch(ido - 1, 3, k) = tr2 - tr1;
    ch(ido - 1, 1, k) = cc(0, k, 0) - cc(0, k, 2);
    ch(0, 2, k) = cc(0, k, 3) - cc(0, k, 1);
  }
  if (ido < 2)
    return;
  if (ido > 2) {
    for (std::size_t k = 0; k < l1; k++) {
      for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
        const std::size_t i1 = i0 + 1;
        const std::size_t ic = ido - i1;
        const T cr2 = wa1[i0 - 1] * cc(i0, k, 1) + wa1[i0] * cc(i1, k, 1);
        const T ci2 = wa1[i0 - 1] * cc(i1, k, 1) - wa1[i0] * cc(i0, k, 1);
        const T cr3 = wa2[i0 - 1] * cc(i0, k, 2) + wa2[i0] * cc(i1, k, 2);
        const T ci3 = wa2[i0 - 1] * cc(i1, k, 2) - wa2[i0] * cc(i0, k, 2);
        const T cr4 = wa3[i0 - 1] * cc(i0, k, 3) + wa3[i0] * cc(i1, k, 3);
        const T ci4 = wa3[i0 - 1] * cc(i1, k, 3) - wa3[i0] * cc(i0, k, 3);
        const T tr1 = cr2 + cr4;
        const T tr4 = cr4 - cr2;
        const T ti1 = ci2 + ci4;
        const T ti4 = ci2 - ci4;
        const T ti2 = cc(i1, k, 0) + ci3;
        const T ti3 = cc(i1, k, 0) - ci3;
        const T tr2 = cc(i0, k, 0) + cr3;
        const T tr3 = cc(i0, k, 0) - cr3;
        ch(i0, 0, k) = tr1 + tr2;
        ch(ic - 1, 3, k) = tr2 - tr1;
        ch(i1, 0, k) = ti1 + ti2;
        ch(ic, 3, k) = ti1 - ti2;
        ch(i0, 2, k) = ti4 + tr3;
        ch(ic - 1, 1, k) = tr3 - ti4;
        ch(i1, 2, k) = tr4 + ti3;
        ch(ic, 1, k) = tr4 - ti3;
      }
    }
    if (ido % 2 != 0)
      return;
  }
  for (std::size_t k = 0; k < l1; k++) {
    const T ti1 = -hsqt2 * (cc(ido - 1, k, 1) + cc(ido - 1, k, 3));
    const T tr1 = hsqt2 * (cc(ido - 1, k, 1) - cc(ido - 1, k, 3));
    ch(ido - 1, 0, k) = tr1 + cc(ido - 1, k, 0);
    ch(ido - 1, 2, k) = cc(ido - 1, k, 0) - tr1;
    ch(0, 1, k) = ti1 - cc(ido - 1, k, 2);
    ch(0, 3, k) = ti1 + cc(ido - 1, k, 2);
  }
}

template <class T>
void passf5(std::size_t ido, std::size_t l1, T *cc_begin, T *ch_begin,
            const E<T> *wa1, const E<T> *wa2, const E<T> *wa3,
            const E<T> *wa4) {
  const Ref3<T> cc{cc_begin, ido, l1};
  const Ref3<T> ch{ch_begin, ido, 5};
  const E<T> tr11 =
      static_cast<E<T>>(.309016994374947424102293417182819058860154589903);
  const E<T> ti11 =
      static_cast<E<T>>(.951056516295153572116439333379382143405698634126);
  const E<T> tr12 =
      static_cast<E<T>>(-.809016994374947424102293417182819058860154589903);
  const E<T> ti12 =
      static_cast<E<T>>(.587785252292473129168705954639072768597652437643);
  for (std::size_t k = 0; k < l1; k++) {
    const T cr2 = cc(0, k, 4) + cc(0, k, 1);
    const T ci5 = cc(0, k, 4) - cc(0, k, 1);
    const T cr3 = cc(0, k, 3) + cc(0, k, 2);
    const T ci4 = cc(0, k, 3) - cc(0, k, 2);
    ch(0, 0, k) = cc(0, k, 0) + cr2 + cr3;
    ch(ido - 1, 1, k) = cc(0, k, 0) + tr11 * cr2 + tr12 * cr3;
    ch(0, 2, k) = ti11 * ci5 + ti12 * ci4;
    ch(ido - 1, 3, k) = cc(0, k, 0) + tr12 * cr2 + tr11 * cr3;
    ch(0, 4, k) = ti12 * ci5 - ti11 * ci4;
  }
  if (ido == 1)
    return;
  for (std::size_t k = 0; k < l1; k++) {
    for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
      const std::size_t i1 = i0 + 1;
      const std::size_t ic = ido - i1;
      const T dr2 = wa1[i0 - 1] * cc(i0, k, 1) + wa1[i0] * cc(i1, k, 1);
      const T di2 = wa1[i0 - 1] * cc(i1, k, 1) - wa1[i0] * cc(i0, k, 1);
      const T dr3 = wa2[i0 - 1] * cc(i0, k, 2) + wa2[i0] * cc(i1, k, 2);
      const T di3 = wa2[i0 - 1] * cc(i1, k, 2) - wa2[i0] * cc(i0, k, 2);
      const T dr4 = wa3[i0 - 1] * cc(i0, k, 3) + wa3[i0] * cc(i1, k, 3);
      const T di4 = wa3[i0 - 1] * cc(i1, k, 3) - wa3[i0] * cc(i0, k, 3);
      const T dr5 = wa4[i0 - 1] * cc(i0, k, 4) + wa4[i0] * cc(i1, k, 4);
      const T di5 = wa4[i0 - 1] * cc(i1, k, 4) - wa4[i0] * cc(i0, k, 4);
      const T cr2 = dr2 + dr5;
      const T ci5 = dr5 - dr2;
      const T cr5 = di2 - di5;
      const T ci2 = di2 + di5;
      const T cr3 = dr3 + dr4;
      const T ci4 = dr4 - dr3;
      const T cr4 = di3 - di4;
      const T ci3 = di3 + di4;
      ch(i0, 0, k) = cc(i0, k, 0) + cr2 + cr3;
      ch(i1, 0, k) = cc(i1, k, 0) + ci2 + ci3;
      const T tr2 = cc(i0, k, 0) + tr11 * cr2 + tr12 * cr3;
      const T ti2 = cc(i1, k, 0) + tr11 * ci2 + tr12 * ci3;
      const T tr3 = cc(i0, k, 0) + tr12 * cr2 + tr11 * cr3;
      const T ti3 = cc(i1, k, 0) + tr12 * ci2 + tr11 * ci3;
      const T tr5 = ti11 * cr5 + ti12 * cr4;
      const T ti5 = ti11 * ci5 + ti12 * ci4;
      const T tr4 = ti12 * cr5 - ti11 * cr4;
      const T ti4 = ti12 * ci5 - ti11 * ci4;
      ch(i0, 2, k) = tr2 + tr5;
      ch(ic - 1, 1, k) = tr2 - tr5;
      ch(i1, 2, k) = ti2 + ti5;
      ch(ic, 1, k) = ti5 - ti2;
      ch(i0, 4, k) = tr3 + tr4;
      ch(ic - 1, 3, k) = tr3 - tr4;
      ch(i1, 4, k) = ti3 + ti4;
      ch(ic, 3, k) = ti4 - ti3;
    }
  }
}

template <class T>
void passfg(std::size_t ido, std::size_t ip, std::size_t l1, std::size_t idl1,
            T *cc_begin, T *c1_begin, T *c2_begin, T *ch_begin, T *ch2_begin,
            const E<T> *wa) {
  const Ref3<T> cc{cc_begin, ido, ip};
  const Ref3<T> c1{c1_begin, ido, l1};
  const Ref2<T> c2{c2_begin, idl1};
  const Ref3<T> ch{ch_begin, ido, l1};
  const Ref2<T> ch2{ch2_begin, idl1};
  const double tpi = 8.0 * std::atan(1.0);
  const double arg = tpi / static_cast<double>(ip);
  double dsp_d, dcp_d;
  ::sincos(arg, &dsp_d, &dcp_d);
  const E<T> dsp = static_cast<E<T>>(dsp_d), dcp = static_cast<E<T>>(dcp_d);
  const std::size_t ipph = (ip + 1) / 2;
  const std::size_t nbd = (ido - 1) / 2;
  if (ido == 1) {
    for (std::size_t ik = 0; ik < idl1; ik++)
      c2(ik, 0) = ch2(ik, 0);
  } else {
    for (std::size_t ik = 0; ik < idl1; ik++)
      ch2(ik, 0) = c2(ik, 0);
    for (std::size_t j = 1; j < ip; j++)
      for (std::size_t k = 0; k < l1; k++)
        ch(0, k, j) = c1(0, k, j);
    if (nbd <= l1) {
      std::size_t is = 0;
      for (std::size_t j = 1; j < ip; j++) {
        std::size_t idij = is;
        for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
          const std::size_t i1 = i0 + 1;
          for (std::size_t k = 0; k < l1; k++) {
            ch(i0, k, j) =
                wa[idij] * c1(i0, k, j) + wa[idij + 1] * c1(i1, k, j);
            ch(i1, k, j) =
                wa[idij] * c1(i1, k, j) - wa[idij + 1] * c1(i0, k, j);
          }
          idij += 2;
        }
        is += ido;
      }
    } else {
      std::size_t is = 0;
      for (std::size_t j = 1; j < ip; j++) {
        for (std::size_t k = 0; k < l1; k++) {
          std::size_t idij = is;
          for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
            const std::size_t i1 = i0 + 1;
            ch(i0, k, j) =
                wa[idij] * c1(i0, k, j) + wa[idij + 1] * c1(i1, k, j);
            ch(i1, k, j) =
                wa[idij] * c1(i1, k, j) - wa[idij + 1] * c1(i0, k, j);
            idij = idij + 2;
          }
        }
        is += ido;
      }
    }
    if (nbd >= l1) {
      for (std::size_t j = 1; j < ipph; j++) {
        const std::size_t jc = ip - j;
        for (std::size_t k = 0; k < l1; k++)
          for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
            const std::size_t i1 = i0 + 1;
            c1(i0, k, j) = ch(i0, k, j) + ch(i0, k, jc);
            c1(i0, k, jc) = ch(i1, k, j) - ch(i1, k, jc);
            c1(i1, k, j) = ch(i1, k, j) + ch(i1, k, jc);
            c1(i1, k, jc) = ch(i0, k, jc) - ch(i0, k, j);
          }
      }
    } else {
      for (std::size_t j = 1; j < ipph; j++) {
        const std::size_t jc = ip - j;
        for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
          const std::size_t i1 = i0 + 1;
          for (std::size_t k = 0; k < l1; k++) {
            c1(i0, k, j) = ch(i0, k, j) + ch(i0, k, jc);
            c1(i0, k, jc) = ch(i1, k, j) - ch(i1, k, jc);
            c1(i1, k, j) = ch(i1, k, j) + ch(i1, k, jc);
            c1(i1, k, jc) = ch(i0, k, jc) - ch(i0, k, j);
          }
        }
      }
    }
  }
  for (std::size_t j = 1; j < ipph; j++) {
    const std::size_t jc = ip - j;
    for (std::size_t k = 0; k < l1; k++) {
      c1(0, k, j) = ch(0, k, j) + ch(0, k, jc);
      c1(0, k, jc) = ch(0, k, jc) - ch(0, k, j);
    }
  }
  E<T> ar1 = 1.;
  E<T> ai1 = 0.;
  for (std::size_t l = 1; l < ipph; l++) {
    const std::size_t lc = ip - l;
    const E<T> ar1h = dcp * ar1 - dsp * ai1;
    ai1 = dcp * ai1 + dsp * ar1;
    ar1 = ar1h;
    for (std::size_t ik = 0; ik < idl1; ik++) {
      ch2(ik, l) = c2(ik, 0) + ar1 * c2(ik, 1);
      ch2(ik, lc) = ai1 * c2(ik, ip - 1);
    }
    const E<T> dc2 = ar1;
    const E<T> ds2 = ai1;
    E<T> ar2 = ar1;
    E<T> ai2 = ai1;
    for (std::size_t j = 2; j < ipph; j++) {
      const std::size_t jc = ip - j;
      const E<T> ar2h = dc2 * ar2 - ds2 * ai2;
      ai2 = dc2 * ai2 + ds2 * ar2;
      ar2 = ar2h;
      for (std::size_t ik = 0; ik < idl1; ik++) {
        ch2(ik, l) = ch2(ik, l) + ar2 * c2(ik, j);
        ch2(ik, lc) = ch2(ik, lc) + ai2 * c2(ik, jc);
      }
    }
  }
  for (std::size_t j = 1; j < ipph; j++)
    for (std::size_t ik = 0; ik < idl1; ik++)
      ch2(ik, 0) = ch2(ik, 0) + c2(ik, j);
  if (ido >= l1) {
    for (std::size_t k = 0; k < l1; k++)
      for (std::size_t i = 0; i < ido; i++)
        cc(i, 0, k) = ch(i, k, 0);
  } else {
    for (std::size_t i = 0; i < ido; i++)
      for (std::size_t k = 0; k < l1; k++)
        cc(i, 0, k) = ch(i, k, 0);
  }
  for (std::size_t j = 1; j < ipph; j++) {
    const std::size_t jc = ip - j;
    const std::size_t j2 = j + j;
    for (std::size_t k = 0; k < l1; k++) {
      cc(ido - 1, j2 - 1, k) = ch(0, k, j);
      cc(0, j2, k) = ch(0, k, jc);
    }
  }
  if (ido == 1)
    return;
  if (nbd >= l1) {
    for (std::size_t j = 1; j < ipph; j++) {
      const std::size_t jc = ip - j;
      const std::size_t j2 = j + j;
      for (std::size_t k = 0; k < l1; k++)
        for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
          const std::size_t i1 = i0 + 1;
          const std::size_t ic = ido - i1;
          cc(i0, j2, k) = ch(i0, k, j) + ch(i0, k, jc);
          cc(ic - 1, j2 - 1, k) = ch(i0, k, j) - ch(i0, k, jc);
          cc(i1, j2, k) = ch(i1, k, j) + ch(i1, k, jc);
          cc(ic, j2 - 1, k) = ch(i1, k, jc) - ch(i1, k, j);
        }
    }
  } else {
    for (std::size_t j = 1; j < ipph; j++) {
      const std::size_t jc = ip - j;
      const std::size_t j2 = j + j;
      for (std::size_t i0 = 1; i0 < ido - 1; i0 += 2) {
        const std::size_t i1 = i0 + 1;
        const std::size_t ic = ido - i1;
        for (std::size_t k = 0; k < l1; k++) {
          cc(i0, j2, k) = ch(i0, k, j) + ch(i0, k, jc);
          cc(ic - 1, j2 - 1, k) = ch(i0, k, j) - ch(i0, k, jc);
          cc(i1, j2, k) = ch(i1, k, j) + ch(i1, k, jc);
          cc(ic, j2 - 1, k) = ch(i1, k, jc) - ch(i1, k, j);
        }
      }
    }
  }
}

} // namespace

RealToComplex::RealToComplex(std::size_t n_real)
    : n_(n_real), factors_(factorize(n_real)), wa_(n_real) {
  compute_twiddles();
  wa_float_.assign(wa_.begin(), wa_.end());
}

void RealToComplex::compute_twiddles() {
  if (n_ < 2)
    return;
  double *wa = wa_.data();
  const double tpi = 8.0 * std::atan(1.0);
  const double argh = tpi / static_cast<double>(n_);
  std::size_t is = 0;
  const std::size_t nfm1 = factors_.size() - 1;
  std::size_t l1 = 1;
  if (nfm1 == 0)
    return;
  for (std::size_t k1 = 0; k1 < nfm1; k1++) {
    const std::size_t ip = static_cast<std::size_t>(factors_[k1]);
    std::size_t ld = 0;
    const std::size_t l2 = l1 * ip;
    const std::size_t ido = n_ / l2;
    const std::size_t ipm = ip - 1;
    for (std::size_t j = 1; j <= ipm; j++) {
      ld = ld + l1;
      std::size_t i = is;
      const double argld = static_cast<double>(ld) * argh;
      double fi = 0.;
      for (std::size_t ii = 3; ii <= ido; ii += 2) {
        i = i + 2;
        fi = fi + 1.;
        const double arg = fi * argld;
        ::sincos(arg, &wa[i - 1], &wa[i - 1 - 1]);
      }
      is = is + ido;
    }
    l1 = l2;
  }
}

template <class T> void RealToComplex::forward_compressed(T *c, T *ch) const {
  if (n_ < 2)
    return;
  const E<T> *wa = twiddles<E<T>>();
  std::size_t na = 1;
  std::size_t l2 = n_;
  std::size_t iw = n_;
  for (std::size_t k1 = 1; k1 <= factors_.size(); k1++) {
    const std::size_t kh = factors_.size() - k1;
    const std::size_t ip = static_cast<std::size_t>(factors_[kh + 3 - 2 - 1]);
    const std::size_t l1 = l2 / ip;
    const std::size_t ido = n_ / l2;
    const std::size_t idl1 = ido * l1;
    iw = iw - (ip - 1) * ido;
    na = 1 - na;
    if (ip == 4) {
      const std::size_t ix2 = iw + ido;
      const std::size_t ix3 = ix2 + ido;
      if (na == 0)
        passf4(ido, l1, c, ch, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);
      else
        passf4(ido, l1, ch, c, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);
    } else if (ip == 2) {
      if (na == 0)
        passf2(ido, l1, c, ch, wa + iw - 1);
      else
        passf2(ido, l1, ch, c, wa + iw - 1);
    } else if (ip == 3) {
      const std::size_t ix2 = iw + ido;
      if (na == 0)
        passf3(ido, l1, c, ch, wa + iw - 1, wa + ix2 - 1);
      else
        passf3(ido, l1, ch, c, wa + iw - 1, wa + ix2 - 1);
    } else if (ip == 5) {
      const std::size_t ix2 = iw + ido;
      const std::size_t ix3 = ix2 + ido;
      const std::size_t ix4 = ix3 + ido;
      if (na == 0)
        passf5(ido, l1, c, ch, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1,
               wa + ix4 - 1);
      else
        passf5(ido, l1, ch, c, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1,
               wa + ix4 - 1);
    } else {
      if (ido == 1)
        na = 1 - na;
      if (na == 0) {
        passfg(ido, ip, l1, idl1, c, c, c, ch, ch, wa + iw - 1);
        na = 1;
      } else {
        passfg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa + iw - 1);
        na = 0;
      }
    }
    l2 = l1;
  }
  if (na == 1)
    return;
  for (std::size_t i = 1; i <= n_; i++)
    c[i - 1] = ch[i - 1];
}

template <class T> void RealToComplex::forward(T *seq, T *scratch) const {
  forward_compressed(seq, scratch);
  std::copy_backward(seq + 1, seq + n_, seq + n_ + 1);
  seq[1] = T{};
  if ((n_ & 1) == 0)
    seq[n_ + 1] = T{};
}

template void RealToComplex::forward(double *, double *) const;
template void RealToComplex::forward(Lanes *, Lanes *) const;
template void RealToComplex::forward(float *, float *) const;
template void RealToComplex::forward(FloatLanes *, FloatLanes *) const;

} // namespace mxi::stills::fftpack

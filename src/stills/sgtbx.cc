// Follows cctbx_project cctbx/sgtbx (see sgtbx.hh). cctbx Copyright (c) 2006 -
// 2026, The Regents of the University of California, through Lawrence
// Berkeley National Laboratory. BSD-3-Clause-LBNL; see
// third_party/cctbx/LICENSE.txt.

#include "stills/sgtbx.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <map>
#include <numeric>
#include <unordered_map>

namespace mxi::stills::sgtbx {

namespace {

[[noreturn]] void fail(const char *what) { throw Error(what); }

void check(bool ok, const char *what = "sgtbx internal error") {
  if (!ok)
    fail(what);
}

int iabs(int x) { return x < 0 ? -x : x; }

int gcd_int(int a, int b) {
  for (;;) {
    if (b == 0)
      return a < 0 ? -a : a;
    const int next_b = a % b;
    a = b;
    b = next_b;
  }
}

int lcm_int(int a, int b) {
  if (a == 0 || b == 0)
    return 0;
  return iabs(a / gcd_int(a, b) * b);
}

int mod_positive(int x, int y) {
  if (y > 0) {
    x %= y;
    if (x < 0)
      x += y;
  }
  return x;
}

int mod_short(int x, int y) {
  x = mod_positive(x, y);
  if (x > y / 2)
    x -= y;
  return x;
}

int determinant(const IntMat3 &m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) -
         m[1] * (m[3] * m[8] - m[5] * m[6]) +
         m[2] * (m[3] * m[7] - m[4] * m[6]);
}

IntMat3 co_factor_matrix_transposed(const IntMat3 &m) {
  return {m[4] * m[8] - m[5] * m[7], -m[1] * m[8] + m[2] * m[7],
          m[1] * m[5] - m[2] * m[4], -m[3] * m[8] + m[5] * m[6],
          m[0] * m[8] - m[2] * m[6], -m[0] * m[5] + m[2] * m[3],
          m[3] * m[7] - m[4] * m[6], -m[0] * m[7] + m[1] * m[6],
          m[0] * m[4] - m[1] * m[3]};
}

IntMat3 matmul(const IntMat3 &a, const IntMat3 &b) {
  IntMat3 r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r[i * 3 + j] =
          a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
  return r;
}

IntVec3 matvec(const IntMat3 &a, const IntVec3 &v) {
  return {a[0] * v[0] + a[1] * v[1] + a[2] * v[2],
          a[3] * v[0] + a[4] * v[1] + a[5] * v[2],
          a[6] * v[0] + a[7] * v[1] + a[8] * v[2]};
}

void set_column(IntMat3 &m, int c, const IntVec3 &v) {
  m[c] = v[0];
  m[3 + c] = v[1];
  m[6 + c] = v[2];
}

// utils::change_denominator.
template <std::size_t N>
bool change_denominator(const std::array<int, N> &old_num, int old_den,
                        std::array<int, N> &new_num, int new_den) {
  for (std::size_t i = 0; i < N; ++i) {
    new_num[i] = old_num[i] * new_den;
    if (new_num[i] % old_den)
      return false;
    new_num[i] /= old_den;
  }
  return true;
}

// utils::cmp_i_vec: fewer zeros first, then by pattern of zeros, then by
// magnitude, then larger first.
bool cmp_i_vec(const int *a, const int *b, std::size_t n) {
  std::size_t n0a = 0, n0b = 0;
  for (std::size_t i = 0; i < n; ++i)
    if (a[i] == 0)
      n0a++;
  for (std::size_t i = 0; i < n; ++i)
    if (b[i] == 0)
      n0b++;
  if (n0a > n0b)
    return true;
  if (n0a < n0b)
    return false;
  for (std::size_t i = 0; i < n; ++i) {
    if (a[i] != 0 && b[i] == 0)
      return true;
    if (a[i] == 0 && b[i] != 0)
      return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (iabs(a[i]) < iabs(b[i]))
      return true;
    if (iabs(a[i]) > iabs(b[i]))
      return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (a[i] > b[i])
      return true;
    if (a[i] < b[i])
      return false;
  }
  return false;
}

// A view of a row-major integer matrix whose shape changes in place, as
// scitbx's af::ref<int, af::mat_grid> does.
struct IntRef {
  int *p = nullptr;
  std::size_t nr = 0, nc = 0;
  int &operator()(std::size_t i, std::size_t j) const { return p[i * nc + j]; }
  std::size_t size() const { return nr * nc; }
  bool is_square() const { return nr == nc; }
  bool is_diagonal() const {
    if (!is_square())
      return false;
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = 0; j < nc; ++j)
        if (i != j && (*this)(i, j))
          return false;
    return true;
  }
  void set_identity() const {
    check(is_square());
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = 0; j < nc; ++j)
        (*this)(i, j) = (i == j) ? 1 : 0;
  }
  void transpose_square_in_place() const {
    check(is_square());
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = i + 1; j < nc; ++j)
        std::swap((*this)(i, j), (*this)(j, i));
  }
  void transpose_in_place() {
    if (is_square()) {
      transpose_square_in_place();
      return;
    }
    std::vector<int> mt(size());
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = 0; j < nc; ++j)
        mt[j * nr + i] = (*this)(i, j);
    std::copy(mt.begin(), mt.end(), p);
    std::swap(nr, nc);
  }
};

// scitbx::matrix::row_echelon::form_t. m shrinks to its rank.
std::size_t row_echelon_form_t(IntRef &m, const IntRef &t) {
  const std::size_t mr = m.nr, mc = m.nc, tc = t.nc;
  if (tc)
    check(t.p != nullptr && t.nr >= mr);
  auto swap_rows = [](const IntRef &x, std::size_t a, std::size_t b,
                      std::size_t n) {
    for (std::size_t c = 0; c < n; ++c)
      std::swap(x(a, c), x(b, c));
  };
  std::size_t i = 0, j = 0;
  while (i < mr && j < mc) {
    std::size_t k = i;
    while (k < mr && m(k, j) == 0)
      k++;
    if (k == mr) {
      j++;
      continue;
    }
    if (i != k) {
      swap_rows(m, i, k, mc);
      if (tc)
        swap_rows(t, i, k, tc);
    }
    for (k++; k < mr; k++) {
      const int a = iabs(m(k, j));
      if (a != 0 && a < iabs(m(i, j))) {
        swap_rows(m, i, k, mc);
        if (tc)
          swap_rows(t, i, k, tc);
      }
    }
    if (m(i, j) < 0) {
      for (std::size_t ic = 0; ic < mc; ++ic)
        m(i, ic) *= -1;
      if (tc)
        for (std::size_t ic = 0; ic < tc; ++ic)
          t(i, ic) *= -1;
    }
    bool cleared = true;
    for (k = i + 1; k < mr; k++) {
      const int a = m(k, j) / m(i, j);
      if (a != 0) {
        for (std::size_t ic = 0; ic < mc; ++ic)
          m(k, ic) -= a * m(i, ic);
        if (tc)
          for (std::size_t ic = 0; ic < tc; ++ic)
            t(k, ic) -= a * t(i, ic);
      }
      if (m(k, j) != 0)
        cleared = false;
    }
    if (cleared) {
      i++;
      j++;
    }
  }
  m.nr = i;
  return i;
}

std::size_t row_echelon_form(IntRef &m) { return row_echelon_form_t(m, {}); }

// scitbx::matrix::row_echelon::back_substitution_int.
int back_substitution_int(const IntRef &re, const int *v, int *sol,
                          bool *flag_indep) {
  const std::size_t nr = re.nr, nc = re.nc;
  if (flag_indep)
    for (std::size_t ic = 0; ic < nc; ++ic)
      flag_indep[ic] = true;
  int d = 1;
  for (std::size_t ir = nr; ir > 0;) {
    ir--;
    std::size_t ic = 0;
    for (; ic < nc; ++ic)
      if (re(ir, ic))
        break;
    if (ic == nc) {
      if (v && v[ir] != 0)
        return 0;
      continue;
    }
    if (flag_indep)
      flag_indep[ic] = false;
    if (sol) {
      const std::size_t icp = ic + 1;
      const std::size_t nv = nc - icp;
      if (nv) {
        int s = 0;
        for (std::size_t k = 0; k < nv; ++k)
          s += re(ir, icp + k) * sol[icp + k];
        sol[ic] = -s;
      } else {
        sol[ic] = 0;
      }
      if (v)
        sol[ic] += d * v[ir];
      const int mrc = re(ir, ic);
      int f = gcd_int(sol[ic], mrc);
      if (mrc < 0)
        f *= -1;
      sol[ic] /= f;
      f = mrc / f;
      if (f != 1) {
        for (std::size_t jc = 0; jc < nc; ++jc)
          if (jc != ic)
            sol[jc] *= f;
        d *= f;
      }
    }
  }
  return d;
}

std::vector<std::size_t> independent_indices(const IntRef &re) {
  bool flag_indep[9];
  check(re.nc <= 9);
  check(back_substitution_int(re, nullptr, nullptr, flag_indep) >= 1);
  std::vector<std::size_t> out;
  for (std::size_t ic = 0; ic < re.nc; ++ic)
    if (flag_indep[ic])
      out.push_back(ic);
  return out;
}

int sign_hemisphere(const IntVec3 &v) {
  if (v[2] > 0)
    return 1;
  if (v[2] == 0) {
    if (v[1] > 0)
      return 1;
    if (v[1] == 0) {
      if (v[0] > 0)
        return 1;
      if (v[0] == 0)
        return 0;
    }
  }
  return -1;
}

IntVec3 homog_rank_2(const IntRef &re) {
  check(re.nr == 2);
  const auto indep = independent_indices(re);
  check(indep.size() == 1);
  IntVec3 ev{0, 0, 0};
  ev[indep[0]] = 1;
  check(back_substitution_int(re, nullptr, ev.data(), nullptr) >= 1);
  if (sign_hemisphere(ev) < 0)
    for (int &x : ev)
      x *= -1;
  return ev;
}

int sense_of_rotation(const RotMx &r, int type, const IntVec3 &ev) {
  const int f = type < 0 ? -1 : 1;
  const IntMat3 &m = r.num();
  const int trace = f * (m[0] + m[4] + m[8]);
  if (trace == 3 || trace == -1)
    return 0;
  if (ev[1] == 0 && ev[2] == 0) {
    if (ev[0] * f * r[7] > 0)
      return 1;
  } else {
    if (f * (r[3] * ev[2] - r[6] * ev[1]) > 0)
      return 1;
  }
  return -1;
}

TrVec tr_vec_12(int a, int b, int c) {
  return TrVec(a, b, c, 12) * (sg_t_den / 12);
}

struct CentringType {
  char symbol;
  std::vector<TrVec> t;
};

const std::vector<CentringType> &centring_types() {
  static const std::vector<CentringType> table = {
      {'P', {tr_vec_12(0, 0, 0)}},
      {'A', {tr_vec_12(0, 0, 0), tr_vec_12(0, 6, 6)}},
      {'B', {tr_vec_12(0, 0, 0), tr_vec_12(6, 0, 6)}},
      {'C', {tr_vec_12(0, 0, 0), tr_vec_12(6, 6, 0)}},
      {'I', {tr_vec_12(0, 0, 0), tr_vec_12(6, 6, 6)}},
      {'R', {tr_vec_12(0, 0, 0), tr_vec_12(8, 4, 4), tr_vec_12(4, 8, 8)}},
      {'Q', {tr_vec_12(0, 0, 0), tr_vec_12(4, 8, 4), tr_vec_12(8, 4, 8)}},
      {'H', {tr_vec_12(0, 0, 0), tr_vec_12(8, 4, 0), tr_vec_12(4, 8, 0)}},
      {'F',
       {tr_vec_12(0, 0, 0), tr_vec_12(0, 6, 6), tr_vec_12(6, 0, 6),
        tr_vec_12(6, 6, 0)}},
  };
  return table;
}

const CentringType *get_centring_type(char symbol) {
  symbol = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol)));
  if (symbol == 'Q')
    return nullptr;
  for (const auto &e : centring_types())
    if (e.symbol == symbol)
      return &e;
  return nullptr;
}

RotMx conventional_z2p_matrix(char symbol) {
  switch (symbol) {
  case 'P':
    return RotMx({1, 0, 0, 0, 1, 0, 0, 0, 1}, 1);
  case 'A':
    return RotMx({-1, 0, 0, 0, -1, 1, 0, 1, 1}, 1);
  case 'B':
    return RotMx({-1, 0, 1, 0, -1, 0, 1, 0, 1}, 1);
  case 'C':
    return RotMx({1, 1, 0, 1, -1, 0, 0, 0, -1}, 1);
  case 'I':
    return RotMx({0, 1, 1, 1, 0, 1, 1, 1, 0}, 1);
  case 'R':
    return RotMx({1, 0, 1, -1, 1, 1, 0, -1, 1}, 1);
  case 'H':
    return RotMx({1, 1, 0, -1, 2, 0, 0, 0, 1}, 1);
  case 'F':
    return RotMx({-1, 1, 1, 1, -1, 1, 1, 1, -1}, 1);
  default:
    return RotMx(0);
  }
}

// cmp_smx, the ordering make_tidy sorts the representative operations by.
int compare_rot_mx(const RotMx &a, const RotMx &b) {
  if (a.den() == b.den()) {
    for (int i = 0; i < 9; ++i) {
      if (a[i] < b[i])
        return -1;
      if (a[i] > b[i])
        return 1;
    }
  } else {
    const int l = lcm_int(a.den(), b.den());
    const RotMx al = a.new_denominator(l), bl = b.new_denominator(l);
    for (int i = 0; i < 9; ++i) {
      if (al[i] < bl[i])
        return -1;
      if (al[i] > bl[i])
        return 1;
    }
  }
  return 0;
}

bool cmp_smx(const RtMx &a, const RtMx &b) {
  bool cmp_as_symop = a.r().den() == 1 && b.r().den() == 1 &&
                      a.r().type() != 0 && b.r().type() != 0;
  RotMxInfo ri_a, ri_b;
  if (cmp_as_symop) {
    try {
      ri_a = RotMxInfo(a.r());
      ri_b = RotMxInfo(b.r());
    } catch (const Error &) {
      cmp_as_symop = false;
    }
    if (cmp_as_symop) {
      if (iabs(ri_a.type) > iabs(ri_b.type))
        return true;
      if (iabs(ri_a.type) < iabs(ri_b.type))
        return false;
      if (ri_a.type > ri_b.type)
        return true;
      if (ri_a.type < ri_b.type)
        return false;
      if (cmp_i_vec(ri_a.ev.data(), ri_b.ev.data(), 3))
        return true;
      if (cmp_i_vec(ri_b.ev.data(), ri_a.ev.data(), 3))
        return false;
      if (ri_a.sense > ri_b.sense)
        return true;
      if (ri_a.sense < ri_b.sense)
        return false;
    }
  }
  if (!cmp_as_symop) {
    const int o = compare_rot_mx(a.r(), b.r());
    if (o < 0)
      return true;
    if (o > 0)
      return false;
  }
  TrVec at = a.t(), bt = b.t();
  if (at.den() != bt.den()) {
    const int l = lcm_int(at.den(), bt.den());
    at = at.new_denominator(l);
    bt = bt.new_denominator(l);
  }
  if (cmp_i_vec(at.num().data(), bt.num().data(), 3))
    return true;
  if (cmp_i_vec(bt.num().data(), at.num().data(), 3))
    return false;
  if (cmp_as_symop) {
    const int o = compare_rot_mx(a.r(), b.r());
    if (o < 0)
      return true;
    if (o > 0)
      return false;
  }
  for (int i = 0; i < 3; ++i) {
    if (at[i] < bt[i])
      return true;
    if (at[i] > bt[i])
      return false;
  }
  return false;
}

bool cmp_ltr(const TrVec &a, const TrVec &b) {
  for (int i = 0; i < 3; ++i) {
    if (a[i] < b[i])
      return true;
    if (a[i] > b[i])
      return false;
  }
  return false;
}

std::string format_rational(int num, int den) {
  const int g = gcd_int(num, den);
  if (g) {
    num /= g;
    den /= g;
  }
  if (den < 0) {
    num = -num;
    den = -den;
  }
  if (den == 1)
    return std::to_string(num);
  return std::to_string(num) + "/" + std::to_string(den);
}

} // namespace

// ---- TrVec ----

TrVec TrVec::new_denominator(int new_den) const {
  TrVec r(new_den);
  if (!change_denominator(num_, den_, r.num_, new_den))
    throw RationalError("Unsuitable value for rational translation vector.");
  return r;
}

TrVec TrVec::scale(int factor) const {
  if (factor == 1)
    return *this;
  return TrVec({num_[0] * factor, num_[1] * factor, num_[2] * factor},
               den_ * factor);
}

TrVec TrVec::mod_positive() const {
  TrVec r(*this);
  for (int i = 0; i < 3; ++i)
    r.num_[i] = sgtbx::mod_positive(r.num_[i], den_);
  return r;
}

TrVec TrVec::mod_short() const {
  TrVec r(*this);
  for (int i = 0; i < 3; ++i)
    r.num_[i] = sgtbx::mod_short(r.num_[i], den_);
  return r;
}

TrVec TrVec::cancel() const {
  int g = den_;
  for (int i = 0; i < 3; ++i)
    g = gcd_int(g, num_[i]);
  if (g == 0)
    return *this;
  return TrVec({num_[0] / g, num_[1] / g, num_[2] / g}, den_ / g);
}

TrVec TrVec::operator+(const TrVec &o) const {
  check(den_ == o.den_);
  return TrVec({num_[0] + o.num_[0], num_[1] + o.num_[1], num_[2] + o.num_[2]},
               den_);
}

TrVec TrVec::operator-(const TrVec &o) const {
  check(den_ == o.den_);
  return TrVec({num_[0] - o.num_[0], num_[1] - o.num_[1], num_[2] - o.num_[2]},
               den_);
}

TrVec TrVec::operator*(int f) const {
  return TrVec({num_[0] * f, num_[1] * f, num_[2] * f}, den_);
}

TrVec TrVec::operator/(int d) const {
  IntVec3 n{};
  for (int i = 0; i < 3; ++i) {
    if (num_[i] % d)
      throw RationalError("Unsuitable value for rational translation vector.");
    n[i] = num_[i] / d;
  }
  return TrVec(n, den_);
}

// ---- RotMx ----

RotMx::RotMx(int den, int diagonal) : num_{}, den_(den) {
  num_[0] = num_[4] = num_[8] = diagonal * den;
}

bool RotMx::is_unit_mx() const {
  return num_ == IntMat3{den_, 0, 0, 0, den_, 0, 0, 0, den_};
}

RotMx RotMx::minus_unit_mx() const {
  RotMx r(*this);
  for (int i = 0; i < 9; i += 4)
    r.num_[i] -= den_;
  return r;
}

RotMx RotMx::new_denominator(int new_den) const {
  RotMx r(new_den);
  if (!change_denominator(num_, den_, r.num_, new_den))
    throw RationalError("Unsuitable value for rational rotation matrix.");
  return r;
}

RotMx RotMx::scale(int factor) const {
  if (factor == 1)
    return *this;
  IntMat3 n = num_;
  for (int &x : n)
    x *= factor;
  return RotMx(n, den_ * factor);
}

RotMx RotMx::inverse() const {
  const int det_den3 = determinant(num_);
  if (det_den3 == 0)
    fail("Rotation matrix is not invertible.");
  IntMat3 n = co_factor_matrix_transposed(num_);
  for (int &x : n)
    x *= den_ * den_;
  return RotMx(n, den_) / det_den3;
}

RotMx RotMx::transpose() const {
  const IntMat3 &m = num_;
  return RotMx({m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]}, den_);
}

RotMx RotMx::operator-() const {
  IntMat3 n = num_;
  for (int &x : n)
    x = -x;
  return RotMx(n, den_);
}

RotMx RotMx::operator+(const RotMx &o) const {
  check(den_ == o.den_);
  IntMat3 n{};
  for (int i = 0; i < 9; ++i)
    n[i] = num_[i] + o.num_[i];
  return RotMx(n, den_);
}

RotMx RotMx::operator-(const RotMx &o) const {
  check(den_ == o.den_);
  IntMat3 n{};
  for (int i = 0; i < 9; ++i)
    n[i] = num_[i] - o.num_[i];
  return RotMx(n, den_);
}

RotMx RotMx::operator*(const RotMx &o) const {
  return RotMx(matmul(num_, o.num_), den_ * o.den_);
}

TrVec RotMx::operator*(const TrVec &t) const {
  return TrVec(matvec(num_, t.num()), den_ * t.den());
}

IntVec3 RotMx::operator*(const IntVec3 &v) const { return matvec(num_, v); }

RotMx RotMx::operator*(int f) const {
  IntMat3 n = num_;
  for (int &x : n)
    x *= f;
  return RotMx(n, den_);
}

RotMx RotMx::operator/(int d) const {
  IntMat3 n{};
  for (int i = 0; i < 9; ++i) {
    if (num_[i] % d)
      throw RationalError("Unsuitable value for rational rotation matrix.");
    n[i] = num_[i] / d;
  }
  return RotMx(n, den_);
}

RotMx RotMx::cancel() const {
  int g = den_;
  for (int i = 0; i < 9; ++i)
    g = gcd_int(g, num_[i]);
  if (g == 0)
    return *this;
  IntMat3 n = num_;
  for (int &x : n)
    x /= g;
  return RotMx(n, den_ / g);
}

int RotMx::type() const {
  const int det = determinant(num_);
  if (det == -1 || det == 1) {
    switch (num_[0] + num_[4] + num_[8]) {
    case -3:
      return -1;
    case -2:
      return -6;
    case -1:
      return det == -1 ? -4 : 2;
    case 0:
      return det == -1 ? -3 : 3;
    case 1:
      return det == -1 ? -2 : 4;
    case 2:
      return 6;
    case 3:
      return 1;
    }
  }
  return 0;
}

int RotMx::order(int t) const {
  if (t == 0)
    t = type();
  if (t > 0)
    return t;
  if (t % 2)
    return -t * 2;
  return -t;
}

RotMx RotMx::accumulate(int t) const {
  check(den_ == 1);
  const int ord = order(t);
  if (ord == 1)
    return *this;
  check(ord != 0);
  IntMat3 result{1, 0, 0, 0, 1, 0, 0, 0, 1};
  IntMat3 a = num_;
  for (int i = 0; i < 9; ++i)
    result[i] += a[i];
  for (int i = 2; i < ord; ++i) {
    a = matmul(a, num_);
    for (int k = 0; k < 9; ++k)
      result[k] += a[k];
  }
  return RotMx(result, 1);
}

RotMxInfo RotMx::info() const { return RotMxInfo(*this); }

Mat3 RotMx::as_double() const {
  Mat3 m;
  for (int i = 0; i < 9; ++i)
    m.m[i] = static_cast<double>(num_[i]) / static_cast<double>(den_);
  return m;
}

RotMxInfo::RotMxInfo(const RotMx &r) : type(r.type()) {
  check(r.den() == 1);
  if (type == 0)
    fail("Cannot determine type of rotation matrix.");
  RotMx proper_r = r;
  int proper_order = type;
  if (proper_order < 0) {
    proper_order *= -1;
    proper_r = -proper_r;
  }
  if (proper_order > 1) {
    RotMx rmi = proper_r.minus_unit_mx();
    IntRef re{rmi.num().data(), 3, 3};
    if (row_echelon_form(re) != 2)
      fail("Cannot determine Eigenvector of rotation matrix.");
    ev = homog_rank_2(re);
    sense = sense_of_rotation(r, type, ev);
  }
}

// ---- RtMx ----

RtMx RtMx::new_denominators(int r_den, int t_den) const {
  RtMx r(*this);
  if (r_den)
    r.r_ = r.r_.new_denominator(r_den);
  if (t_den)
    r.t_ = r.t_.new_denominator(t_den);
  return r;
}

RtMx RtMx::scale(int factor_r, int factor_t) const {
  if (factor_t == 0)
    factor_t = factor_r;
  return RtMx(r_.scale(factor_r), t_.scale(factor_t));
}

RtMx RtMx::inverse() const {
  const RotMx r_inv = r_.inverse();
  const TrVec t_inv = ((-r_inv) * t_).new_denominator(t_.den());
  return RtMx(r_inv, t_inv);
}

RtMx RtMx::operator*(const RtMx &o) const {
  check(t_.den() == o.t_.den());
  return RtMx(r_ * o.r_, r_ * o.t_ + t_.scale(r_.den()));
}

TrVec RtMx::t_intrinsic_part() const {
  const int ty = r_.type();
  return (r_.accumulate(ty) * t_) / r_.order(ty);
}

bool RtMx::operator<(const RtMx &o) const { return cmp_smx(*this, o); }

std::string RtMx::as_xyz() const {
  // scitbx::matrix::rational_as_xyz, fractional constants, translation last.
  static const char letters[3] = {'x', 'y', 'z'};
  std::string out;
  for (int row = 0; row < 3; ++row) {
    if (row)
      out += ",";
    std::string term;
    for (int col = 0; col < 3; ++col) {
      const int n = r_[row * 3 + col];
      if (n == 0)
        continue;
      const int g = gcd_int(n, r_.den());
      int nn = n / g, dd = r_.den() / g;
      std::string s;
      if (nn < 0) {
        s = "-";
        nn = -nn;
      } else if (!term.empty()) {
        s = "+";
      }
      if (nn != 1)
        s += std::to_string(nn) + "*";
      s += letters[col];
      if (dd != 1)
        s += "/" + std::to_string(dd);
      term += s;
    }
    if (t_.is_valid() && t_[row] != 0) {
      const std::string f = format_rational(t_[row], t_.den());
      if (f[0] != '-' && !term.empty())
        term += "+";
      term += f;
    }
    if (term.empty())
      term = "0";
    out += term;
  }
  return out;
}

// ---- ChangeOfBasisOp ----

RotMx ChangeOfBasisOp::operator()(const RotMx &r) const {
  check(r.den() == 1);
  return (c_.r() * r * c_inv_.r()).new_denominator(1);
}

RtMx ChangeOfBasisOp::operator()(const RtMx &s) const {
  check(s.r().den() == 1);
  check(c_.t().den() % s.t().den() == 0);
  return (c_ * (s.scale(1, c_.t().den() / s.t().den()) * c_inv_))
      .new_denominators(s);
}

TrVec ChangeOfBasisOp::operator()(const TrVec &t, int sign_identity) const {
  const TrVec tf = t.new_denominator(c_inv_.t().den());
  const TrVec tw = sign_identity >= 0 ? tf + c_inv_.t() : tf - c_inv_.t();
  return (c_.r() * tw + c_.t().scale(c_.r().den())).new_denominator(t.den());
}

void ChangeOfBasisOp::update(const ChangeOfBasisOp &other) {
  c_ = (other.c() * c_).new_denominators(other.c());
  c_inv_ = (c_inv_ * other.c_inv()).new_denominators(other.c_inv());
}

void ChangeOfBasisOp::update(const TrVec &shift) {
  c_ = RtMx(c_.r(), c_.t() + shift);
  c_inv_ =
      RtMx(c_inv_.r(),
           c_inv_.t() - (c_inv_.r() * shift).new_denominator(c_inv_.t().den()));
}

ChangeOfBasisOp ChangeOfBasisOp::operator*(const ChangeOfBasisOp &rhs) const {
  return ChangeOfBasisOp((c() * rhs.c()).new_denominators(c()),
                         (rhs.c_inv() * c_inv()).new_denominators(c_inv()));
}

// ---- TrGroup ----

bool TrGroup::contains(const TrVec &t) const {
  if (!t.is_valid())
    return false;
  return std::find(elems_.begin(), elems_.end(), t.mod_positive()) !=
         elems_.end();
}

bool TrGroup::add(const TrVec &new_t) {
  if (!new_t.is_valid())
    return false;
  const TrVec t = new_t.mod_positive();
  if (std::find(elems_.begin(), elems_.end(), t) != elems_.end())
    return false;
  check(t.den() == elems_[0].den());
  elems_.push_back(t);
  return true;
}

bool TrGroup::expand(const TrVec &new_t) {
  TrVec trial = new_t;
  const std::size_t old_size = size();
  std::size_t i = old_size, j = 1;
  for (;;) {
    add(trial);
    if (j > i) {
      i++;
      j = 1;
    }
    if (i == size())
      break;
    trial = elems_[j] + elems_[i];
    j++;
  }
  return size() != old_size;
}

TrGroup TrGroup::change_basis(const ChangeOfBasisOp &cb_op) const {
  const int t_den = elems_[0].den();
  TrGroup result(t_den);
  for (int i = 0; i < 3; ++i) {
    TrVec bv(t_den);
    bv[i] = t_den;
    result.expand(cb_op(bv, 1));
  }
  for (std::size_t i = 1; i < size(); ++i)
    result.expand(cb_op(elems_[i], 1));
  return result;
}

char TrGroup::conventional_centring_type_symbol() const {
  for (const auto &e : centring_types()) {
    if (e.t.size() != size())
      continue;
    std::vector<bool> match(size(), false);
    std::size_t n_matches = 0;
    for (std::size_t i = 0; i < size(); ++i)
      for (std::size_t j = 0; j < size(); ++j)
        if (!match[j] && e.t[i].new_denominator(t_den()) == elems_[j]) {
          match[j] = true;
          n_matches++;
          break;
        }
    if (n_matches == size())
      return e.symbol;
  }
  return '\0';
}

ChangeOfBasisOp TrGroup::conventional_z2p_op(int r_den, int t_den) const {
  const RotMx m = conventional_z2p_matrix(conventional_centring_type_symbol());
  if (!m.is_valid())
    return ChangeOfBasisOp(0, 0);
  return ChangeOfBasisOp(RtMx(m.new_denominator(r_den), TrVec(t_den)));
}

TrVec TrGroup::tidy(const TrVec &t) const {
  const int lcm_den = lcm_int(elems_[0].den(), t.den());
  const int f_den = lcm_den / elems_[0].den();
  const TrVec t_lcm = t.scale(lcm_den / t.den());
  TrVec t_best = t_lcm.mod_short();
  for (std::size_t i = 1; i < size(); ++i) {
    const TrVec trial = (t_lcm + elems_[i].scale(f_den)).mod_short();
    if (cmp_i_vec(trial.num().data(), t_best.num().data(), 3))
      t_best = trial;
  }
  return t_best.new_denominator(t.den()).mod_positive();
}

void TrGroup::sort_tail() {
  std::sort(elems_.begin() + 1, elems_.end(), cmp_ltr);
}

// ---- MatrixGroup ----

const char *MatrixGroup::label() const {
  static const char *names[] = {
      "Undefined", "Unknown", "1",    "-1",    "2",   "m",    "2/m",
      "222",       "mm2",     "mmm",  "4",     "-4",  "4/m",  "422",
      "4mm",       "-4m2",    "-42m", "4/mmm", "3",   "-3",   "321",
      "312",       "32",      "3m1",  "31m",   "3m",  "-3m1", "-31m",
      "-3m",       "6",       "-6",   "6/m",   "622", "6mm",  "-6m2",
      "-62m",      "6/mmm",   "23",   "m-3",   "432", "-43m", "m-3m"};
  return names[m];
}

// ---- SpaceGroup ----

SpaceGroup::SpaceGroup(bool no_expand, int t_den) : no_expand_(no_expand) {
  reset(t_den);
}

void SpaceGroup::reset(int t_den) {
  n_lsl_ = 1;
  n_ssl_ = 1;
  f_inv_ = 1;
  ltr_ = TrGroup(t_den);
  inv_t_ = TrVec(0);
  smx_.clear();
  smx_.push_back(RtMx(1, t_den));
  is_tidy_ = false;
}

void SpaceGroup::add_inv(const TrVec &new_inv_t) {
  if (is_centric()) {
    if (ltr_.add(inv_t_ - new_inv_t))
      is_tidy_ = false;
    return;
  }
  inv_t_ = new_inv_t.mod_positive();
  f_inv_ = 2;
  is_tidy_ = false;
  if (!no_expand_)
    for (std::size_t i = 1; i < n_smx(); ++i)
      if (ltr_.add(smx_[i].r() * inv_t_ + 2 * smx_[i].t() - inv_t_))
        is_tidy_ = false;
}

void SpaceGroup::add_smx(const RtMx &new_smx) {
  const RotMx minus_r = -new_smx.r();
  for (std::size_t i = 0; i < n_smx(); ++i) {
    if (smx_[i].r() == new_smx.r()) {
      if (ltr_.add(smx_[i].t() - new_smx.t()))
        is_tidy_ = false;
      return;
    }
    if (smx_[i].r() == minus_r) {
      add_inv(smx_[i].t() + new_smx.t());
      return;
    }
  }
  const int d = determinant(new_smx.r().num());
  if (n_smx() >= n_max_repr_rot_mx || (d != -1 && d != 1))
    throw NonCrystallographic();
  check(new_smx.t().den() == smx_[0].t().den());
  smx_.push_back(new_smx.mod_positive());
  const RtMx &s = smx_.back();
  if (!no_expand_ && is_centric())
    ltr_.add(s.r() * inv_t_ + 2 * s.t() - inv_t_);
  is_tidy_ = false;
}

SpaceGroup &SpaceGroup::expand_ltr(const TrVec &new_ltr) {
  if (no_expand_) {
    if (ltr_.add(new_ltr))
      is_tidy_ = false;
    return *this;
  }
  for (std::size_t i_smx = n_ssl_; i_smx < n_smx(); ++i_smx)
    for (std::size_t i_ltr = 1; i_ltr < n_lsl_; ++i_ltr)
      if (ltr_.add(smx_[i_smx].r() * ltr_[i_ltr]))
        is_tidy_ = false;
  n_ssl_ = n_smx();
  TrVec trial = new_ltr;
  std::size_t i = n_lsl_, j = 1;
  for (;;) {
    if (ltr_.add(trial))
      is_tidy_ = false;
    for (std::size_t i_smx = 1; i_smx < n_smx(); ++i_smx)
      for (std::size_t i_ltr = n_lsl_; i_ltr < ltr_.size(); ++i_ltr)
        if (ltr_.add(smx_[i_smx].r() * ltr_[i_ltr]))
          is_tidy_ = false;
    n_lsl_ = ltr_.size();
    if (j > i) {
      i++;
      j = 1;
    }
    if (i == ltr_.size())
      break;
    trial = ltr_[j] + ltr_[i];
    j++;
  }
  return *this;
}

SpaceGroup &SpaceGroup::expand_inv(const TrVec &new_inv_t) {
  add_inv(new_inv_t);
  expand_ltr(TrVec(0));
  return *this;
}

SpaceGroup &SpaceGroup::expand_smx(const RtMx &new_smx) {
  if (new_smx.r().den() != 1)
    fail("sgtbx::space_group::expand_smx(): rotation-part denominator must "
         "be 1 (implementation limitation).");
  if (new_smx.t().den() != smx_[0].t().den())
    fail("sgtbx::space_group::expand_smx(): incompatible translation-part "
         "denominator.");
  if (no_expand_) {
    add_smx(new_smx);
    return *this;
  }
  RtMx trial = new_smx;
  std::size_t i = n_smx(), j = 1;
  for (;;) {
    add_smx(trial);
    if (j > i) {
      i++;
      j = 1;
    }
    if (i == n_smx())
      break;
    trial = smx_[j] * smx_[i];
    j++;
  }
  expand_ltr(TrVec(0));
  return *this;
}

std::size_t SpaceGroup::expand_conventional_centring_type(char symbol) {
  const CentringType *e = get_centring_type(symbol);
  if (e == nullptr)
    fail("Illegal symbol for centring type of cell.");
  for (const auto &t : e->t)
    expand_ltr(t.new_denominator(t_den()));
  return e->t.size();
}

TrVec SpaceGroup::inv_t(bool tidy) const {
  if (!tidy || !inv_t_.is_valid())
    return inv_t_;
  return ltr_.tidy(inv_t_);
}

RtMx SpaceGroup::operator()(std::size_t i_ltr, std::size_t i_inv,
                            std::size_t i_smx) const {
  check(i_ltr < ltr_.size() && i_inv < f_inv_ && i_smx < n_smx(),
        "space group index out of range");
  if (i_inv == 0)
    return smx_[i_smx] + ltr_[i_ltr];
  return (-smx_[i_smx] + inv_t_) + ltr_[i_ltr];
}

RtMx SpaceGroup::operator()(std::size_t i_op) const {
  check(i_op < order_z(), "space group index out of range");
  const std::size_t i_smx = i_op % n_smx();
  const std::size_t i_inv = (i_op / n_smx()) % f_inv_;
  const std::size_t i_ltr = i_op / (f_inv_ * n_smx());
  return (*this)(i_ltr, i_inv, i_smx);
}

SpaceGroup &SpaceGroup::make_tidy() {
  if (is_tidy_)
    return *this;
  if (is_centric()) {
    inv_t_ = inv_t(true);
    for (std::size_t i = 1; i < n_smx(); ++i)
      if (determinant(smx_[i].r().num()) < 0)
        smx_[i] = smx_[i].pre_multiply_inv_t(inv_t_);
  }
  for (std::size_t i = 1; i < n_smx(); ++i)
    smx_[i] = RtMx(smx_[i].r(), ltr_.tidy(smx_[i].t()));
  if (n_ltr() > 2)
    ltr_.sort_tail();
  if (n_smx() > 2)
    std::sort(smx_.begin() + 1, smx_.end(), cmp_smx);
  is_tidy_ = true;
  return *this;
}

bool SpaceGroup::operator==(const SpaceGroup &rhs) const {
  check(r_den() == rhs.r_den() && t_den() == rhs.t_den());
  if (n_ltr() != rhs.n_ltr() || f_inv() != rhs.f_inv() ||
      n_smx() != rhs.n_smx())
    return false;
  SpaceGroup a = *this, b = rhs;
  a.make_tidy();
  b.make_tidy();
  if (a.inv_t_ != b.inv_t_)
    return false;
  if (!(a.ltr_ == b.ltr_))
    return false;
  for (std::size_t i = 0; i < a.n_smx(); ++i)
    if (a.smx(i) != b.smx(i))
      return false;
  return true;
}

std::size_t KeyHash::operator()(const Key &key) const noexcept {
  std::size_t h = 14695981039346656037ull;
  for (int v : key) {
    h ^= static_cast<std::size_t>(static_cast<unsigned>(v));
    h *= 1099511628211ull;
  }
  return h;
}

void append_key(Key &key, const RotMx &r) {
  key.push_back(r.den());
  for (std::size_t i = 0; i < 9; ++i)
    key.push_back(r[i]);
}

void append_key(Key &key, const TrVec &t) {
  key.push_back(t.den());
  for (std::size_t i = 0; i < 3; ++i)
    key.push_back(t[i]);
}

void append_key(Key &key, const RtMx &s) {
  append_key(key, s.r());
  append_key(key, s.t());
}

void append_key(Key &key, const ChangeOfBasisOp &cb_op) {
  append_key(key, cb_op.c());
  append_key(key, cb_op.c_inv());
}

void SpaceGroup::append_key(Key &key) const {
  key.push_back(no_expand_);
  key.push_back(static_cast<int>(n_lsl_));
  key.push_back(static_cast<int>(n_ssl_));
  key.push_back(static_cast<int>(f_inv_));
  key.push_back(is_tidy_);
  key.push_back(static_cast<int>(ltr_.size()));
  for (const TrVec &t : ltr_.elems())
    sgtbx::append_key(key, t);
  sgtbx::append_key(key, inv_t_);
  key.push_back(static_cast<int>(smx_.size()));
  for (const RtMx &s : smx_)
    sgtbx::append_key(key, s);
}

SpaceGroup SpaceGroup::change_basis(const ChangeOfBasisOp &cb_op) const {
  thread_local std::unordered_map<Key, SpaceGroup, KeyHash> cache;
  Key key;
  append_key(key);
  sgtbx::append_key(key, cb_op);
  if (const auto it = cache.find(key); it != cache.end())
    return it->second;
  SpaceGroup result = compute_change_basis(cb_op);
  if (cache.size() >= kCacheLimit)
    cache.clear();
  cache.emplace(std::move(key), result);
  return result;
}

SpaceGroup
SpaceGroup::compute_change_basis(const ChangeOfBasisOp &cb_op) const {
  SpaceGroup result(no_expand_, t_den());
  result.ltr_ = ltr_.change_basis(cb_op);
  if (is_centric())
    result.expand_inv(cb_op(inv_t_, -1));
  for (std::size_t i = 1; i < n_smx(); ++i)
    result.expand_smx(cb_op(smx_[i]));
  return result;
}

SpaceGroup SpaceGroup::build_derived_group(bool discard_z, bool add_inv) const {
  SpaceGroup result(false, t_den());
  if (!discard_z)
    for (std::size_t i = 0; i < n_ltr(); ++i)
      result.expand_ltr(ltr_[i]);
  if (is_centric() || add_inv)
    result.expand_inv(TrVec(t_den()));
  for (std::size_t i = 0; i < n_smx(); ++i)
    result.expand_smx(RtMx(smx_[i].r(), t_den()));
  return result;
}

SpaceGroup SpaceGroup::build_derived_acentric_group() const {
  if (!is_centric())
    return *this;
  SpaceGroup tidy_sg = *this;
  tidy_sg.make_tidy();
  SpaceGroup result(false, tidy_sg.t_den());
  for (std::size_t i = 0; i < tidy_sg.n_ltr(); ++i)
    result.expand_ltr(tidy_sg.ltr_[i]);
  for (std::size_t i = 0; i < tidy_sg.n_smx(); ++i)
    result.expand_smx(tidy_sg.smx_[i]);
  return result;
}

namespace {

bool first_is_shorter(const IntVec3 &a, const IntVec3 &b) {
  for (int i = 0; i < 3; ++i)
    if (a[i]) {
      if (iabs(a[i]) > iabs(b[i]))
        return false;
      return true;
    }
  return true;
}

IntVec3 cross(const IntVec3 &a, const IntVec3 &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// build_list_tot_tr: the centring vectors with 1 subtracted from their
// nonzero components, one of each linearly dependent set (the shorter),
// sorted, then the unit translations.
std::vector<TrVec> build_list_tot_tr(const TrGroup &ltr, int t_den) {
  std::vector<TrVec> tlt;
  for (std::size_t i_ltr = 1; i_ltr < ltr.size(); ++i_ltr) {
    IntVec3 n_utr{1, 1, 1};
    for (int i = 0; i < 3; ++i)
      if (ltr[i_ltr][i])
        n_utr[i] = 2;
    IntVec3 u{};
    for (u[0] = 0; u[0] < n_utr[0]; u[0]++)
      for (u[1] = 0; u[1] < n_utr[1]; u[1]++)
        for (u[2] = 0; u[2] < n_utr[2]; u[2]++) {
          TrVec v = ltr[i_ltr] - TrVec(u, 1).new_denominator(ltr[0].den());
          v = v.new_denominator(t_den);
          std::size_t i = 0;
          for (; i < tlt.size(); ++i)
            if (cross(tlt[i].num(), v.num()) == IntVec3{0, 0, 0}) {
              if (!first_is_shorter(tlt[i].num(), v.num()))
                tlt[i] = v;
              break;
            }
          if (i == tlt.size())
            tlt.push_back(v);
        }
  }
  std::sort(tlt.begin(), tlt.end(), [](const TrVec &a, const TrVec &b) {
    return cmp_i_vec(a.num().data(), b.num().data(), 3);
  });
  for (int i = 0; i < 3; ++i) {
    TrVec v(t_den);
    v[i] = t_den;
    tlt.push_back(v);
  }
  return tlt;
}

} // namespace

ChangeOfBasisOp SpaceGroup::construct_z2p_op(int r_den, int t_den) const {
  const int r_den_3 = r_den * r_den * r_den;
  const std::vector<TrVec> tlt = build_list_tot_tr(ltr_, r_den);
  IntMat3 basis{};
  for (std::size_t i0 = 0; i0 + 2 < tlt.size(); ++i0) {
    set_column(basis, 0, tlt[i0].num());
    for (std::size_t i1 = i0 + 1; i1 + 1 < tlt.size(); ++i1) {
      set_column(basis, 1, tlt[i1].num());
      for (std::size_t i2 = i1 + 1; i2 < tlt.size(); ++i2) {
        set_column(basis, 2, tlt[i2].num());
        const int f = determinant(basis) * static_cast<int>(n_ltr());
        if (f == r_den_3 || -f == r_den_3) {
          if (f < 0)
            for (int i = 0; i < 3; ++i)
              basis[i * 3] *= -1;
          ChangeOfBasisOp result;
          SpaceGroup primitive_sg(false, this->t_den());
          try {
            result =
                ChangeOfBasisOp(RtMx(RotMx(basis, r_den), t_den)).inverse();
            primitive_sg = change_basis(result);
          } catch (const Error &) {
            continue;
          }
          if (primitive_sg.n_ltr() == 1) {
            check(determinant(result.c().r().num()) ==
                  static_cast<int>(n_ltr()) * r_den_3);
            return result;
          }
        }
      }
    }
  }
  fail("construct_z2p_op: no primitive basis");
}

ChangeOfBasisOp SpaceGroup::z2p_op(int r_den, int t_den) const {
  const ChangeOfBasisOp cb_op = ltr_.conventional_z2p_op(r_den, t_den);
  if (cb_op.is_valid())
    return cb_op;
  return construct_z2p_op(r_den, t_den);
}

MatrixGroup SpaceGroup::point_group_type() const {
  using namespace matrix_group;
  std::map<int, int> counter;
  for (std::size_t i = 0; i < n_smx(); ++i)
    counter[smx_[i].r().type()]++;
  const bool centric = is_centric();
  const std::size_t n = n_smx();
  if (counter[-3] + counter[3] == 8) {
    if (n == 12)
      return centric ? code_m3b : code_23;
    if (n == 24) {
      if (!centric) {
        if (counter[4] == 6)
          return code_432;
        if (counter[-4] == 6)
          return code_4b3m;
      } else {
        return code_m3bm;
      }
    }
  } else if (counter[-6] + counter[6] == 2) {
    if (n == 6) {
      if (!centric) {
        if (counter[6] == 2)
          return code_6;
        if (counter[-6] == 2)
          return code_6b;
      } else {
        return code_6_m;
      }
    } else if (n == 12) {
      if (!centric) {
        if (counter[6] == 2) {
          if (counter[2] == 7)
            return code_622;
          if (counter[-2] == 6)
            return code_6mm;
        } else if (counter[-6] == 2) {
          return code_6bm2;
        }
      } else {
        return code_6_mmm;
      }
    }
  } else if (counter[-3] + counter[3] == 2) {
    if (n == 3)
      return centric ? code_3b : code_3;
    if (n == 6) {
      if (!centric) {
        if (counter[2] == 3)
          return code_32;
        if (counter[-2] == 3)
          return code_3m;
      } else {
        return code_3bm;
      }
    }
  } else if (counter[-4] + counter[4] == 2) {
    if (n == 4) {
      if (!centric) {
        if (counter[4] == 2)
          return code_4;
        if (counter[-4] == 2)
          return code_4b;
      } else {
        return code_4_m;
      }
    } else if (n == 8) {
      if (!centric) {
        if (counter[4] == 2) {
          if (counter[2] == 5)
            return code_422;
          if (counter[-2] == 4)
            return code_4mm;
        } else if (counter[-4] == 2) {
          return code_4bm2;
        }
      } else {
        return code_4_mmm;
      }
    }
  } else if (counter[-2] + counter[2] == 3) {
    if (!centric) {
      if (counter[2] == 3)
        return code_222;
      if (counter[-2] == 2)
        return code_mm2;
    } else {
      return code_mmm;
    }
  } else if (counter[-2] + counter[2] == 1) {
    if (!centric) {
      if (counter[2] == 1)
        return code_2;
      if (counter[-2] == 1)
        return code_m;
    } else {
      return code_2_m;
    }
  } else if (n == 1) {
    return centric ? code_1b : code_1;
  }
  fail("point_group_type: internal error");
}

Sym6 SpaceGroup::average_metrical_matrix(const Sym6 &t) const {
  // average_tensor(smx, g, reciprocal_space=false):
  // sum(sym_mat3::tensor_transpose_transform(R)) / n.
  Sym6 result{0, 0, 0, 0, 0, 0};
  for (const RtMx &s : smx_) {
    const Mat3 cm = s.r().as_double();
    const auto &c = cm.m;
    const double ctt[9] = {c[0] * t[0] + c[3] * t[3] + c[6] * t[4],
                           c[3] * t[1] + c[0] * t[3] + c[6] * t[5],
                           c[6] * t[2] + c[0] * t[4] + c[3] * t[5],
                           c[1] * t[0] + c[4] * t[3] + c[7] * t[4],
                           c[4] * t[1] + c[1] * t[3] + c[7] * t[5],
                           c[7] * t[2] + c[1] * t[4] + c[4] * t[5],
                           c[2] * t[0] + c[5] * t[3] + c[8] * t[4],
                           c[5] * t[1] + c[2] * t[3] + c[8] * t[5],
                           c[8] * t[2] + c[2] * t[4] + c[5] * t[5]};
    const Sym6 v = {ctt[0] * c[0] + ctt[1] * c[3] + ctt[2] * c[6],
                    ctt[3] * c[1] + ctt[4] * c[4] + ctt[5] * c[7],
                    ctt[6] * c[2] + ctt[7] * c[5] + ctt[8] * c[8],
                    ctt[0] * c[1] + ctt[1] * c[4] + ctt[2] * c[7],
                    ctt[0] * c[2] + ctt[1] * c[5] + ctt[2] * c[8],
                    ctt[3] * c[2] + ctt[4] * c[5] + ctt[5] * c[8]};
    for (int i = 0; i < 6; ++i)
      result[i] += v[i];
  }
  const double n = static_cast<double>(smx_.size());
  for (double &x : result)
    x /= n;
  return result;
}

UnitCell SpaceGroup::average_unit_cell(const UnitCell &cell) const {
  return UnitCell::from_metric(average_metrical_matrix(cell.metric));
}

// ---- Hall symbols ----

namespace {

struct ParseString {
  std::string s;
  std::size_t i = 0, mark = 0;
  char operator()() const { return i < s.size() ? s[i] : '\0'; }
  void skip(std::size_t n = 1) { i += n; }
  const char *peek() const { return s.c_str() + i; }
};

bool hall_is_space(char c) {
  if (c == '\0')
    return false;
  if (c == '_')
    return true;
  return std::isspace(static_cast<unsigned char>(c)) != 0;
}

bool hall_is_char(char c) { return c != '\0' && !hall_is_space(c); }

int get_abs_order(char c) {
  switch (c) {
  case '1':
    return 1;
  case '2':
    return 2;
  case '3':
    return 3;
  case '4':
    return 4;
  case '6':
    return 6;
  default:
    return 0;
  }
}

int get_screw(char c) {
  if (c >= '1' && c <= '5')
    return c - '0';
  return 0;
}

char get_reference_axis(char c) {
  c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (c == 'x' || c == 'y' || c == 'z')
    return c;
  return '\0';
}

char get_direction_code(char c) {
  if (c == '\'' || c == '"' || c == '*')
    return c;
  if (c == ',' || c == '.')
    return '\'';
  if (c == ';' || c == ':')
    return '"';
  return '\0';
}

TrVec get_translation(char symbol) {
  switch (std::tolower(static_cast<unsigned char>(symbol))) {
  case 'a':
    return tr_vec_12(6, 0, 0);
  case 'b':
    return tr_vec_12(0, 6, 0);
  case 'c':
    return tr_vec_12(0, 0, 6);
  case 'n':
    return tr_vec_12(6, 6, 6);
  case 'u':
    return tr_vec_12(3, 0, 0);
  case 'v':
    return tr_vec_12(0, 3, 0);
  case 'w':
    return tr_vec_12(0, 0, 3);
  case 'd':
    return tr_vec_12(3, 3, 3);
  default:
    return TrVec(0);
  }
}

const RotMx r_1_000({1, 0, 0, 0, 1, 0, 0, 0, 1}, 1);
const RotMx r_2_001({-1, 0, 0, 0, -1, 0, 0, 0, 1}, 1);
const RotMx r_2_1b0({0, -1, 0, -1, 0, 0, 0, 0, -1}, 1);
const RotMx r_2_110({0, 1, 0, 1, 0, 0, 0, 0, -1}, 1);
const RotMx r_3_001({0, -1, 0, 1, -1, 0, 0, 0, 1}, 1);
const RotMx r_3_111({0, 0, 1, 1, 0, 0, 0, 1, 0}, 1);
const RotMx r_3i111({0, 1, 0, 0, 0, 1, 1, 0, 0}, 1);
const RotMx r_4_001({0, -1, 0, 1, 0, 0, 0, 0, 1}, 1);
const RotMx r_6_001({1, -1, 0, 1, 0, 0, 0, 0, 1}, 1);

RotMx get_rot_mx(bool improper, int abs_order, char reference_axis,
                 char direction_code) {
  struct Entry {
    int order;
    char direction_code;
    const RotMx *r;
  };
  const Entry tab[] = {
      {1, '\0', &r_1_000}, {2, '\0', &r_2_001}, {2, '\'', &r_2_1b0},
      {2, '"', &r_2_110},  {3, '\0', &r_3_001}, {3, '*', &r_3_111},
      {4, '\0', &r_4_001}, {6, '\0', &r_6_001},
  };
  for (const auto &e : tab)
    if (e.order == abs_order && e.direction_code == direction_code) {
      const RotMx r = improper ? -(*e.r) : *e.r;
      if (reference_axis == 'x')
        return r_3_111 * r * r_3i111;
      if (reference_axis == 'y')
        return r_3i111 * r * r_3_111;
      return r;
    }
  fail("get_rot_mx: internal error");
}

// The short change-of-basis form "(0 0 4)": an origin shift in twelfths.
TrVec parse_short_cb_op(ParseString &hs, int t_den) {
  check(t_den % 12 == 0);
  TrVec result(t_den);
  for (int row = 0; row < 3; ++row) {
    while (hall_is_space(hs()))
      hs.skip();
    if (row && hs() == ',') {
      hs.skip();
      while (hall_is_space(hs()))
        hs.skip();
    }
    if (hs() == ')' || hs() == '\0')
      return TrVec(0);
    char *end = nullptr;
    const long v = std::strtol(hs.peek(), &end, 10);
    if (end == hs.peek())
      return TrVec(0);
    hs.skip(static_cast<std::size_t>(end - hs.peek()));
    result[row] = static_cast<int>(v) * (t_den / 12);
  }
  return result;
}

// rt_mx.cpp's rationalize: v * den rounded half away from zero, refused
// when it is more than 0.0005 from that.
bool rationalize(double v, int &out, int den) {
  if (den == 0)
    return false;
  v *= den;
  out = v < 0. ? static_cast<int>(v - .5) : static_cast<int>(v + .5);
  v -= out;
  v /= den;
  if (v < 0.)
    v = -v;
  return !(v > .0005);
}

[[noreturn]] void parse_error(const char *what) {
  throw Error(std::string("Parse error: ") + what);
}

// cctbx's rt_mx_from_string: "x+1/2,y,z", "h,k,l" or "a+c,b,c" up to one of
// stop_chars (or the end).
struct ParsedRtMx {
  RtMx s;
  bool have_xyz = false, have_hkl = false, have_abc = false;
};

ParsedRtMx parse_rt_mx(ParseString &input, const char *stop_chars, int r_den,
                       int t_den) {
  ParsedRtMx out{RtMx(r_den, t_den)};
  RtMx &s = out.s;
  int row = 0, column = -1, sign = 1, mult = 0;
  double val_r[3] = {0., 0., 0.};
  double val_t = 0., value = 0.;
  bool have_value = false;
  constexpr unsigned P_Add = 0x01u, P_Mult = 0x02u, P_Value = 0x04u,
                     P_XYZ = 0x08u, P_Comma = 0x10u;
  unsigned p_mode = P_Add | P_Value | P_XYZ;
  auto is_stop = [&](char c) { return std::strchr(stop_chars, c) != nullptr; };
  for (;; input.skip()) {
    const char c = input();
    if (is_stop(c) || !std::isspace(static_cast<unsigned char>(c))) {
      const char k = is_stop(c) ? '\0' : c;
      switch (k) {
      case '_':
        break;
      case '+':
      case '-':
        sign = k == '+' ? 1 : -1;
        if ((p_mode & P_Add) == 0)
          parse_error("unexpected character");
        if (column >= 0)
          val_r[column] += value;
        else
          val_t += value;
        value = 0.;
        have_value = false;
        column = -1;
        mult = 0;
        p_mode = P_Value | P_XYZ;
        break;
      case '*':
        if ((p_mode & P_Mult) == 0)
          parse_error("unexpected character");
        mult = 1;
        p_mode = P_Value | P_XYZ;
        break;
      case '/':
      case ':':
        if ((p_mode & P_Mult) == 0)
          parse_error("unexpected character");
        mult = -1;
        p_mode = P_Value;
        break;
      case '.':
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9': {
        if ((p_mode & P_Value) == 0)
          parse_error("unexpected character");
        const char *begin = input.peek();
        char *end = nullptr;
        double v = std::strtod(begin, &end);
        if (end == begin)
          parse_error("unexpected character");
        input.skip(static_cast<std::size_t>(end - begin) - 1);
        if (sign == -1) {
          v = -v;
          sign = 1;
        }
        if (mult == 1)
          value *= v;
        else if (mult == -1) {
          if (v != 0.)
            value /= v;
          else if (value != 0.)
            parse_error("unexpected character");
        } else
          value = v;
        have_value = true;
        p_mode = P_Comma | P_Add | P_Mult | P_XYZ;
        break;
      }
      case 'X':
      case 'x':
      case 'Y':
      case 'y':
      case 'Z':
      case 'z':
      case 'H':
      case 'h':
      case 'K':
      case 'k':
      case 'L':
      case 'l':
      case 'A':
      case 'a':
      case 'B':
      case 'b':
      case 'C':
      case 'c': {
        const char lower =
            static_cast<char>(std::tolower(static_cast<unsigned char>(k)));
        if (lower == 'x' || lower == 'y' || lower == 'z') {
          column = lower - 'x';
          out.have_xyz = true;
        } else if (lower == 'h' || lower == 'k' || lower == 'l') {
          column = lower == 'h' ? 0 : lower == 'k' ? 1 : 2;
          out.have_hkl = true;
        } else {
          column = lower - 'a';
          out.have_abc = true;
        }
        if ((out.have_xyz && out.have_hkl) || (out.have_xyz && out.have_abc) ||
            (out.have_hkl && out.have_abc))
          parse_error("mix of notations");
        if ((p_mode & P_XYZ) == 0)
          parse_error("unexpected character");
        if (!have_value) {
          value = sign;
          sign = 1;
        }
        p_mode = P_Comma | P_Add | P_Mult;
        break;
      }
      case ',':
      case ';':
        if (row == 2)
          parse_error("too many row expressions");
        [[fallthrough]];
      case '\0':
        if ((p_mode & P_Comma) == 0)
          parse_error("unexpected end of input");
        if (column >= 0)
          val_r[column] += value;
        else
          val_t += value;
        for (int i = 0; i < 3; ++i)
          if (!rationalize(val_r[i],
                           s.r()[static_cast<std::size_t>(row * 3 + i)], r_den))
            throw RationalError(
                "Unsuitable value for rational rotation matrix.");
        if (!rationalize(val_t, s.t()[static_cast<std::size_t>(row)], t_den))
          throw RationalError(
              "Unsuitable value for rational translation vector.");
        row++;
        column = -1;
        sign = 1;
        mult = 0;
        val_r[0] = val_r[1] = val_r[2] = 0.;
        val_t = 0.;
        value = 0.;
        have_value = false;
        p_mode = P_Add | P_Value | P_XYZ;
        break;
      default:
        parse_error("unexpected character");
      }
    }
    if (is_stop(input()))
      break;
  }
  if (row != 3)
    parse_error("not enough row expressions");
  return out;
}

// change_of_basis_op(symbol, stop_chars, r_den, t_den).
ChangeOfBasisOp parse_change_of_basis_op(ParseString &input,
                                         const char *stop_chars, int r_den,
                                         int t_den) {
  const ParsedRtMx p = parse_rt_mx(input, stop_chars, r_den, t_den);
  if (p.have_hkl) {
    if (!p.s.t().is_zero())
      throw Error("h,k,l change of basis with a translation part");
    const RtMx c_inv(p.s.r().transpose(), p.s.t());
    return ChangeOfBasisOp(c_inv.inverse(), c_inv);
  }
  if (p.have_abc) {
    const RtMx c_inv(p.s.r().transpose(), p.s.t());
    return ChangeOfBasisOp(c_inv.inverse(), c_inv);
  }
  return ChangeOfBasisOp(p.s);
}

} // namespace

ChangeOfBasisOp ChangeOfBasisOp::from_string(const std::string &symbol,
                                             int r_den, int t_den) {
  ParseString ps{symbol};
  return parse_change_of_basis_op(ps, "", r_den, t_den);
}

std::size_t SpaceGroup::parse_hall_symbol(const std::string &symbol,
                                          bool pedantic,
                                          bool no_centring_type_symbol) {
  ParseString hs{symbol};
  std::size_t n_added_mx = 0;
  if (!no_centring_type_symbol) {
    while (hall_is_space(hs()))
      hs.skip();
    if (hs() == '-') {
      expand_inv(TrVec(t_den()));
      hs.skip();
      n_added_mx++;
    }
    if (hs() == '\0')
      fail("Lattice type not specified.");
    n_added_mx += expand_conventional_centring_type(hs());
    hs.skip();
  }
  const char char_after_lattice_type_symbol = hs();
  while (hall_is_space(hs()))
    hs.skip();
  if (hs() == '\0' || hs() == '(') {
    if (pedantic)
      fail("Matrix symbol expected.");
    if (hs() == '\0')
      return n_added_mx;
  }
  if (!no_centring_type_symbol && pedantic &&
      !hall_is_space(char_after_lattice_type_symbol))
    fail("Space expected after lattice type symbol.");

  int i_mx_symbol = 0;
  int first_abs_order = 0;
  char first_reference_axis = '\0';
  while (hs() != '\0' && hs() != '(') {
    bool improper = false;
    int abs_order = 0, screw = 0;
    char reference_axis = '\0', direction_code = '\0';
    TrVec smx_t(t_den());
    if (hs() == '-') {
      improper = true;
      hs.skip();
      if (!hall_is_char(hs()))
        fail("Incomplete matrix symbol.");
    }
    abs_order = get_abs_order(hs());
    if (!abs_order)
      fail("Improper symbol for rotational order.");
    hs.skip();
    screw = get_screw(hs());
    if (screw) {
      if (screw >= abs_order)
        fail("Improper screw translation.");
      hs.skip();
    }
    while (hall_is_char(hs())) {
      if (reference_axis == '\0') {
        reference_axis = get_reference_axis(hs());
        if (reference_axis != '\0') {
          if (abs_order == 1 || (abs_order == 3 && direction_code == '*'))
            fail("Inconsistent matrix symbol.");
          hs.skip();
          continue;
        }
      } else if (get_reference_axis(hs()) != '\0') {
        fail("Multiple axis symbols.");
      }
      if (direction_code == '\0') {
        direction_code = get_direction_code(hs());
        if (direction_code != '\0') {
          if (!(abs_order == 2 &&
                (direction_code == '"' || direction_code == '\'')) &&
              !(abs_order == 3 && direction_code == '*'))
            fail("Inconsistent matrix symbol.");
          if (screw)
            fail("Screw translation for non-principal direction.");
          hs.skip();
          continue;
        }
      } else if (get_direction_code(hs()) != '\0') {
        fail("Multiple axis symbols.");
      }
      const TrVec ht = get_translation(hs());
      if (ht.is_valid()) {
        smx_t = (smx_t + ht.new_denominator(t_den())).mod_positive();
        hs.skip();
        continue;
      }
      if (hs() == '(') {
        if (pedantic)
          fail("Space expected before change-of-basis operator.");
        break;
      }
      fail("Malformed matrix symbol.");
    }
    if (reference_axis == '\0') {
      if (i_mx_symbol == 0) {
        if (abs_order != 1 && !(abs_order == 3 && direction_code == '*'))
          reference_axis = 'z';
      } else if (i_mx_symbol == 1) {
        if (abs_order == 2) {
          if (first_abs_order == 2 || first_abs_order == 4) {
            if (direction_code == '\0')
              reference_axis = 'x';
          } else if (first_abs_order == 3 || first_abs_order == 6) {
            if (direction_code == '\0')
              direction_code = '\'';
            reference_axis = first_reference_axis;
          }
        } else if (abs_order == 3 &&
                   (first_abs_order == 2 || first_abs_order == 4) &&
                   direction_code == '\0') {
          direction_code = '*';
        }
      } else if (i_mx_symbol == 2) {
        if (abs_order == 3 && direction_code == '\0')
          direction_code = '*';
      }
    }
    if (reference_axis == '\0' &&
        (direction_code == '"' || direction_code == '\''))
      reference_axis = 'z';
    if (reference_axis == '\0' && abs_order != 1 && direction_code != '*')
      fail("Need explicit axis symbol.");
    const RotMx smx_r =
        get_rot_mx(improper, abs_order, reference_axis, direction_code);
    if (screw) {
      const int i = reference_axis == 'x' ? 0 : reference_axis == 'y' ? 1 : 2;
      check(smx_t.den() * screw % abs_order == 0);
      smx_t[i] += smx_t.den() * screw / abs_order;
    }
    expand_smx(RtMx(smx_r, smx_t));
    if (i_mx_symbol == 0) {
      first_abs_order = abs_order;
      first_reference_axis = reference_axis;
    }
    i_mx_symbol++;
    if (improper || abs_order != 1)
      n_added_mx++;
    while (hall_is_space(hs()))
      hs.skip();
  }
  ChangeOfBasisOp cb_op(0, 0);
  if (hs() == '(') {
    hs.skip();
    hs.mark = hs.i;
    const TrVec t = parse_short_cb_op(hs, cb_t_den);
    if (t.is_valid()) {
      cb_op = ChangeOfBasisOp(RtMx(t, cb_r_den));
    } else {
      hs.i = hs.mark;
      cb_op = parse_change_of_basis_op(hs, ")", cb_r_den, cb_t_den);
    }
    while (hall_is_space(hs()))
      hs.skip();
    if (hs() != ')')
      fail("Closing parenthesis expected after change-of-basis operator");
    hs.skip();
  }
  while (hall_is_space(hs()))
    hs.skip();
  if (hs() != '\0')
    fail("Unexpected extra character.");
  if (cb_op.is_valid())
    *this = change_basis(cb_op);
  return n_added_mx;
}

SpaceGroup SpaceGroup::from_hall(const std::string &hall, bool pedantic,
                                 bool no_centring_type_symbol, bool no_expand,
                                 int t_den) {
  SpaceGroup g(no_expand, t_den);
  g.parse_hall_symbol(hall, pedantic, no_centring_type_symbol);
  return g;
}

std::vector<RtMx> normalizer_addl_generators(int sg_number, bool affine,
                                             bool k2l, bool l2n) {
  check(0 < sg_number && sg_number <= 230);
  std::vector<RtMx> result;
  for (int i_type = 0; i_type < 2; ++i_type) {
    const char *symbol = nullptr;
    if (i_type == 0 && k2l)
      symbol = normalizer_table(sg_number).k2l;
    else if (i_type == 1 && l2n && (sg_number >= 75 || affine))
      symbol = normalizer_table(sg_number).l2n;
    if (symbol) {
      SpaceGroup g(true);
      const std::size_t n_added = g.parse_hall_symbol(symbol, true, true);
      check(n_added > 0);
      check(g.n_ltr() == 1);
      if (g.is_centric())
        result.push_back(g(0, 1, 0));
      for (std::size_t i = 1; i < g.n_smx(); ++i)
        result.push_back(g.smx(i));
    }
  }
  return result;
}

std::string bravais_symbol(const SpaceGroup &g) {
  char family = '?';
  switch (g.crystal_system()) {
  case CrystalSystem::triclinic:
    family = 'a';
    break;
  case CrystalSystem::monoclinic:
    family = 'm';
    break;
  case CrystalSystem::orthorhombic:
    family = 'o';
    break;
  case CrystalSystem::tetragonal:
    family = 't';
    break;
  case CrystalSystem::trigonal:
  case CrystalSystem::hexagonal:
    family = 'h';
    break;
  case CrystalSystem::cubic:
    family = 'c';
    break;
  default:
    fail("bravais_symbol: undefined crystal system");
  }
  std::string out(1, family);
  out += g.conventional_centring_type_symbol();
  return out;
}

std::vector<RtMx> find_affine(const SpaceGroup &group, int range) {
  check(group.n_smx() > 1, "find_affine: the P1 algorithm is not ported");
  SpaceGroup tidy_group = group;
  tidy_group.make_tidy();
  const std::size_t n_rows = group.order_z() * 9;
  std::vector<int> m;
  m.reserve(n_rows * 9);
  for (std::size_t i_op = 0; i_op < group.order_z(); ++i_op) {
    const RotMx r = group(i_op).r();
    check(r.den() == 1);
    const IntMat3 &rn = r.num();
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) {
        const std::size_t prev = m.size();
        m.resize(prev + 9, 0);
        int *row = m.data() + prev;
        for (int j = 0; j < 3; ++j) {
          if (i == j && j == k)
            continue;
          row[i * 3 + j] += rn[j * 3 + k];
          row[j * 3 + k] -= rn[i * 3 + j];
        }
      }
  }
  IntRef m_ref{m.data(), n_rows, 9};
  row_echelon_form(m_ref);
  const std::vector<std::size_t> indep = independent_indices(m_ref);
  std::vector<RtMx> out;
  const std::size_t n = indep.size();
  if (n == 0)
    return out;
  std::vector<int> loop(n, -range);
  for (;;) {
    RotMx c(1, 0);
    for (std::size_t i = 0; i < n; ++i)
      c[indep[i]] = loop[i];
    const int den =
        back_substitution_int(m_ref, nullptr, c.num().data(), nullptr);
    check(den > 0);
    if (determinant(c.num()) == 1) {
      int i = 0;
      for (; i < 9; ++i)
        if (c[i] < -range || c[i] > range)
          break;
      if (i == 9) {
        const ChangeOfBasisOp cb_op((RtMx(c, sg_t_den)));
        if (tidy_group.n_smx() == 1 ||
            tidy_group.change_basis(cb_op) == tidy_group)
          out.push_back(cb_op.c());
      }
    }
    std::size_t k = n;
    for (;;) {
      if (k == 0)
        return out;
      k--;
      loop[k]++;
      if (loop[k] < range + 1)
        break;
      loop[k] = -range;
    }
  }
}

} // namespace mxi::stills::sgtbx

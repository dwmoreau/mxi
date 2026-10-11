#include "stills/reduction.hh"

#include <cmath>
#include <string>
#include <utility>

namespace mxi::stills {

namespace {

IntMat3 multiply(const IntMat3 &a, const IntMat3 &b) {
  IntMat3 r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r[3 * i + j] =
          a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
  return r;
}

class Niggli {
public:
  Niggli(const UnitCell &cell, double relative_epsilon, std::size_t limit)
      : limit_(limit) {
    const Sym6 &g = cell.metric;
    a_ = g[0];
    b_ = g[1];
    c_ = g[2];
    d_ = 2 * g[5];
    e_ = 2 * g[4];
    f_ = 2 * g[3];
    epsilon_ = std::pow(cell.volume, 1. / 3.) * relative_epsilon;
    while (step()) {
    }
  }

  NiggliResult result() const { return {r_inv_, n_}; }

private:
  bool lt(double x, double y) const { return x < y - epsilon_; }
  bool gt(double x, double y) const { return y < x - epsilon_; }
  bool eq(double x, double y) const { return !(lt(x, y) || lt(y, x)); }

  // (number of epsilon-zero, number of epsilon-positive) among d, e, f.
  std::pair<int, int> def_test() const {
    int n_zero = 0, n_positive = 0;
    for (double x : {d_, e_, f_}) {
      if (lt(0., x))
        ++n_positive;
      else if (!lt(x, 0.))
        ++n_zero;
    }
    return {n_zero, n_positive};
  }

  bool def_gt_0() const {
    const auto [n_zero, n_positive] = def_test();
    return n_positive == 3 || (n_zero == 0 && n_positive == 1);
  }

  void update(const IntMat3 &m) {
    if (n_ == limit_)
      throw NiggliLimitExceeded(
          "Krivy-Gruber iteration limit exceeded (limit=" +
          std::to_string(limit_) + ").");
    r_inv_ = multiply(r_inv_, m);
    ++n_;
  }

  static void require(bool positive) {
    if (!positive)
      throw CellError(
          "Niggli reduction reached a non-positive squared length.");
  }

  void n1() {
    update({0, -1, 0, -1, 0, 0, 0, 0, -1});
    std::swap(a_, b_);
    std::swap(d_, e_);
  }

  void n2() {
    update({-1, 0, 0, 0, 0, -1, 0, -1, 0});
    std::swap(b_, c_);
    std::swap(e_, f_);
  }

  void n3_true() {
    IntMat3 m{1, 0, 0, 0, 1, 0, 0, 0, 1};
    if (lt(d_, 0.))
      m[0] = -1;
    if (lt(e_, 0.))
      m[4] = -1;
    if (lt(f_, 0.))
      m[8] = -1;
    update(m);
    d_ = std::abs(d_);
    e_ = std::abs(e_);
    f_ = std::abs(f_);
  }

  // The last epsilon-zero element takes the sign that makes the product of
  // the three signs positive.
  void n3_false() {
    int s[3] = {1, 1, 1};
    int z = -1;
    const double v[3] = {d_, e_, f_};
    for (int i = 0; i < 3; ++i) {
      if (gt(v[i], 0.))
        s[i] = -1;
      else if (!lt(v[i], 0.))
        z = i;
    }
    if (s[0] * s[1] * s[2] < 0) {
      if (z == -1)
        throw CellError("Niggli reduction found no epsilon-zero element.");
      s[z] = -1;
    }
    update({s[0], 0, 0, 0, s[1], 0, 0, 0, s[2]});
    d_ = -std::abs(d_);
    e_ = -std::abs(e_);
    f_ = -std::abs(f_);
  }

  void a5() {
    if (d_ > 0) {
      update({1, 0, 0, 0, 1, -1, 0, 0, 1});
      c_ += b_ - d_;
      d_ -= b_ + b_;
      e_ -= f_;
    } else {
      update({1, 0, 0, 0, 1, 1, 0, 0, 1});
      c_ += b_ + d_;
      d_ += b_ + b_;
      e_ += f_;
    }
    require(c_ > 0);
  }

  void a6() {
    if (e_ > 0) {
      update({1, 0, -1, 0, 1, 0, 0, 0, 1});
      c_ += a_ - e_;
      d_ -= f_;
      e_ -= a_ + a_;
    } else {
      update({1, 0, 1, 0, 1, 0, 0, 0, 1});
      c_ += a_ + e_;
      d_ += f_;
      e_ += a_ + a_;
    }
    require(c_ > 0);
  }

  void a7() {
    if (f_ > 0) {
      update({1, -1, 0, 0, 1, 0, 0, 0, 1});
      b_ += a_ - f_;
      d_ -= e_;
      f_ -= a_ + a_;
    } else {
      update({1, 1, 0, 0, 1, 0, 0, 0, 1});
      b_ += a_ + f_;
      d_ += e_;
      f_ += a_ + a_;
    }
    require(b_ > 0);
  }

  void a8() {
    update({1, 0, 1, 0, 1, 1, 0, 0, 1});
    c_ += a_ + b_ + d_ + e_ + f_;
    d_ += b_ + b_ + f_;
    e_ += a_ + a_ + f_;
    require(c_ > 0);
  }

  bool step() {
    if (gt(a_, b_) || (eq(a_, b_) && gt(std::abs(d_), std::abs(e_))))
      n1();
    if (gt(b_, c_) || (eq(b_, c_) && gt(std::abs(e_), std::abs(f_)))) {
      n2();
      return true;
    }
    if (def_gt_0())
      n3_true();
    else
      n3_false();
    if (gt(std::abs(d_), b_) || (eq(d_, b_) && lt(e_ + e_, f_)) ||
        (eq(d_, -b_) && lt(f_, 0.))) {
      a5();
      return true;
    }
    if (gt(std::abs(e_), a_) || (eq(e_, a_) && lt(d_ + d_, f_)) ||
        (eq(e_, -a_) && lt(f_, 0.))) {
      a6();
      return true;
    }
    if (gt(std::abs(f_), a_) || (eq(f_, a_) && lt(d_ + d_, e_)) ||
        (eq(f_, -a_) && lt(e_, 0.))) {
      a7();
      return true;
    }
    if (lt(d_ + e_ + f_ + a_ + b_, 0.) ||
        (eq(d_ + e_ + f_ + a_ + b_, 0.) && gt(a_ + a_ + e_ + e_ + f_, 0.))) {
      a8();
      return true;
    }
    return false;
  }

  std::size_t limit_;
  IntMat3 r_inv_{1, 0, 0, 0, 1, 0, 0, 0, 1};
  std::size_t n_ = 0;
  double epsilon_ = 0.0;
  double a_ = 0.0, b_ = 0.0, c_ = 0.0, d_ = 0.0, e_ = 0.0, f_ = 0.0;
};

} // namespace

NiggliResult niggli_reduction(const UnitCell &cell, double relative_epsilon,
                              std::size_t iteration_limit) {
  return Niggli(cell, relative_epsilon, iteration_limit).result();
}

Mat3 niggli_change_of_basis(const NiggliResult &r) {
  const auto &m = r.r_inv;
  const int det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
                  m[1] * (m[3] * m[8] - m[5] * m[6]) +
                  m[2] * (m[3] * m[7] - m[4] * m[6]);
  if (det != 1 && det != -1)
    throw CellError("Niggli change of basis is not unimodular.");
  // A unimodular matrix's inverse is its transposed cofactors times det.
  IntMat3 inv{m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8],
              m[1] * m[5] - m[2] * m[4], m[5] * m[6] - m[3] * m[8],
              m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
              m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7],
              m[0] * m[4] - m[1] * m[3]};
  for (int &x : inv)
    x *= det;
  return {double(inv[0]), double(inv[3]), double(inv[6]),
          double(inv[1]), double(inv[4]), double(inv[7]),
          double(inv[2]), double(inv[5]), double(inv[8])};
}

} // namespace mxi::stills

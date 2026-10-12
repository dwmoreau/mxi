// The space-group type and the change of basis to the reference setting:
// Grosse-Kunstleve (1999) section 5, then cctbx's search over the affine
// normalizer for the operator its comparator ranks first. Follows
// cctbx_project cctbx/sgtbx/space_group_type.cpp and select_generators.cpp.
// cctbx Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory. BSD-3-Clause-LBNL;
// see third_party/cctbx/LICENSE.txt.

#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <utility>

#include "stills/sgtbx.hh"

namespace mxi::stills::sgtbx {

namespace {

[[noreturn]] void fail(const char *what) { throw Error(what); }

void check(bool ok, const char *what = "sgtbx internal error") {
  if (!ok)
    fail(what);
}

int iabs(int x) { return x < 0 ? -x : x; }

int determinant(const IntMat3 &m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) -
         m[1] * (m[3] * m[8] - m[5] * m[6]) +
         m[2] * (m[3] * m[7] - m[4] * m[6]);
}

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

IntVec3 get_column(const IntMat3 &m, int c) {
  return {m[c], m[3 + c], m[6 + c]};
}

void set_column(IntMat3 &m, int c, const IntVec3 &v) {
  m[c] = v[0];
  m[3 + c] = v[1];
  m[6 + c] = v[2];
}

void swap_columns(IntMat3 &m, int a, int b) {
  for (int r = 0; r < 3; ++r)
    std::swap(m[r * 3 + a], m[r * 3 + b]);
}

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

// The row-major integer matrix view the Smith normal form works on; its shape
// changes in place over the same buffer, as scitbx's af::ref does.
struct IntRef {
  int *p = nullptr;
  std::size_t nr = 0, nc = 0;
  int &operator()(std::size_t i, std::size_t j) const { return p[i * nc + j]; }
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
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = 0; j < nc; ++j)
        (*this)(i, j) = (i == j) ? 1 : 0;
  }
  void transpose_square_in_place() const {
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = i + 1; j < nc; ++j)
        std::swap((*this)(i, j), (*this)(j, i));
  }
  void transpose_in_place() {
    if (is_square()) {
      transpose_square_in_place();
      return;
    }
    std::vector<int> mt(nr * nc);
    for (std::size_t i = 0; i < nr; ++i)
      for (std::size_t j = 0; j < nc; ++j)
        mt[j * nr + i] = (*this)(i, j);
    std::copy(mt.begin(), mt.end(), p);
    std::swap(nr, nc);
  }
};

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

std::array<IntVec3, 4> homog_rank_1(const IntRef &re) {
  check(re.nr == 1);
  bool flag_indep[3];
  check(back_substitution_int(re, nullptr, nullptr, flag_indep) >= 1);
  std::vector<std::size_t> indep;
  for (std::size_t ic = 0; ic < re.nc; ++ic)
    if (flag_indep[ic])
      indep.push_back(ic);
  check(indep.size() == 2);
  static const int trial_v[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
  std::array<IntVec3, 4> sol{};
  for (int t = 0; t < 4; ++t) {
    sol[t] = {0, 0, 0};
    for (int i = 0; i < 2; ++i)
      sol[t][indep[i]] = trial_v[t][i];
    check(back_substitution_int(re, nullptr, sol[t].data(), nullptr) > 0);
  }
  return sol;
}

// ---- construct_cb_op_r ----

MatrixGroup get_matrix_group_type(const SpaceGroup &std_sg,
                                  const MatrixGroup &point_group) {
  using namespace matrix_group;
  bool two_fold = false, mirror = false;
  if (point_group == code_4bm2 || point_group == code_6bm2) {
    two_fold = true;
  } else if (std_sg.n_ltr() == 1) {
    if (point_group == code_32)
      two_fold = true;
    else if (point_group == code_3m)
      mirror = true;
    else if (point_group == code_3bm)
      two_fold = mirror = true;
  }
  if (!(two_fold || mirror))
    return point_group;
  for (std::size_t i = 1; i < std_sg.n_smx(); ++i) {
    const int r_type = std_sg.smx(i).r().type();
    check(r_type != 0);
    if ((r_type == 2 && two_fold) || (r_type == -2 && mirror)) {
      const RotMxInfo ri(std_sg.smx(i).r());
      if (ri.ev == IntVec3{1, 0, 0}) {
        if (point_group == code_4bm2)
          return code_4b2m;
        if (point_group == code_32)
          return code_321;
        if (point_group == code_3m)
          return code_3m1;
        if (point_group == code_3bm)
          return code_3bm1;
        if (point_group == code_6bm2)
          return code_6b2m;
      }
    }
  }
  if (point_group == code_4bm2)
    return code_4bm2;
  if (point_group == code_32)
    return code_312;
  if (point_group == code_3m)
    return code_31m;
  if (point_group == code_3bm)
    return code_3b1m;
  if (point_group == code_6bm2)
    return code_6bm2;
  fail("get_matrix_group_type: internal error");
}

struct BasisBuilder {
  CrystalSystem xtal_system;
  std::size_t n_wtd = 0;
  int ord[3] = {0, 0, 0};
  RotMx r[3] = {RotMx(0), RotMx(0), RotMx(0)};
  RotMxInfo ri[3];

  explicit BasisBuilder(const MatrixGroup &laue)
      : xtal_system(laue.crystal_system()) {
    using namespace matrix_group;
    if (laue == code_1b) {
      n_wtd = 0;
    } else if (laue == code_2_m) {
      n_wtd = 1;
      ord[0] = 2;
    } else if (laue == code_mmm) {
      n_wtd = 3;
      ord[0] = ord[1] = ord[2] = 2;
    } else if (laue == code_4_m) {
      n_wtd = 1;
      ord[0] = 4;
    } else if (laue == code_4_mmm) {
      n_wtd = 2;
      ord[0] = 4;
      ord[1] = 2;
    } else if (laue == code_3b) {
      n_wtd = 1;
      ord[0] = 3;
    } else if (laue == code_3bm) {
      n_wtd = 2;
      ord[0] = 3;
      ord[1] = 2;
    } else if (laue == code_6_m) {
      n_wtd = 1;
      ord[0] = 3;
    } else if (laue == code_6_mmm) {
      n_wtd = 2;
      ord[0] = 3;
      ord[1] = 2;
    } else if (laue == code_m3b) {
      n_wtd = 2;
      ord[0] = 3;
      ord[1] = 2;
    } else if (laue == code_m3bm) {
      n_wtd = 2;
      ord[0] = 3;
      ord[1] = 4;
    } else {
      fail("basis_builder: internal error");
    }
  }

  void pick_eigen_vectors(const SpaceGroup &raw_sg) {
    for (std::size_t i = 0; i < n_wtd; ++i)
      r[i] = RotMx(0);
    std::size_t n_found = 0;
    for (;;) {
      bool restart = false;
      for (std::size_t i_smx = 0; i_smx < raw_sg.n_smx(); ++i_smx) {
        const int r_type = raw_sg.smx(i_smx).r().type();
        check(r_type != 0);
        for (std::size_t i_wtd = 0; i_wtd < n_wtd; ++i_wtd) {
          if (!r[i_wtd].is_valid() &&
              (r_type == ord[i_wtd] || -r_type == ord[i_wtd])) {
            ri[i_wtd] = raw_sg.smx(i_smx).r().info();
            if (ri[i_wtd].sense >= 0) {
              bool use_this_smx = true;
              for (std::size_t j_wtd = 0; j_wtd < n_wtd; ++j_wtd) {
                if (r[j_wtd].is_valid() && ri[i_wtd].ev == ri[j_wtd].ev) {
                  if (iabs(ord[j_wtd]) >= iabs(r_type)) {
                    use_this_smx = false;
                  } else {
                    r[j_wtd] = RotMx(0);
                    n_found--;
                    restart = true;
                  }
                  break;
                }
              }
              if (use_this_smx) {
                ord[i_wtd] = r_type;
                r[i_wtd] = raw_sg.smx(i_smx).r();
                n_found++;
                if (n_found == n_wtd)
                  return;
              }
            }
            break;
          }
        }
      }
      if (!restart)
        break;
    }
    fail("pick_eigen_vectors: internal error");
  }

  RotMx uni_axial_basis() {
    if (ord[0] < 0) {
      ord[0] *= -1;
      r[0] = r[0] * -1;
      ri[0] = r[0].info();
    }
    RotMx cum_mx = r[0].accumulate(ord[0]);
    IntRef re{cum_mx.num().data(), 3, 3};
    check(row_echelon_form(re) == 1);
    const std::array<IntVec3, 4> sol = homog_rank_1(re);
    const std::size_t n_ix = ord[0] == 2 ? 2 : 1;
    RotMx trial_basis;
    set_column(trial_basis.num(), 2, ri[0].ev);
    int min_det = 0;
    RotMx result;
    // math::loop_n_from_m<2>(4, n_ix): the n_ix-subsets of the four solutions
    // in lexicographic order.
    std::vector<std::vector<std::size_t>> subsets;
    if (n_ix == 1) {
      for (std::size_t a = 0; a < 4; ++a)
        subsets.push_back({a});
    } else {
      for (std::size_t a = 0; a < 4; ++a)
        for (std::size_t b = a + 1; b < 4; ++b)
          subsets.push_back({a, b});
    }
    for (const auto &ix : subsets) {
      for (std::size_t i = 0; i < ix.size(); ++i)
        set_column(trial_basis.num(), static_cast<int>(i), sol[ix[i]]);
      if (ix.size() == 1)
        set_column(trial_basis.num(), 1, r[0] * sol[ix[0]]);
      const int det = determinant(trial_basis.num());
      if (det != 0 && (min_det == 0 || iabs(min_det) > iabs(det))) {
        min_det = det;
        result = trial_basis;
      }
    }
    check(min_det != 0);
    if (min_det < 0)
      swap_columns(result.num(), 0, 1);
    return result;
  }

  RotMx two_axial_basis() {
    RotMx result;
    set_column(result.num(), 0, ri[1].ev);
    if (xtal_system == CrystalSystem::cubic) {
      for (int i = 0; i < 2; ++i)
        set_column(result.num(), i + 1, r[0] * get_column(result.num(), i));
      if (determinant(result.num()) < 0)
        swap_columns(result.num(), 0, 1);
    } else {
      set_column(result.num(), 2, ri[0].ev);
      if (n_wtd == 3) {
        set_column(result.num(), 1, ri[2].ev);
        if (determinant(result.num()) < 0)
          swap_columns(result.num(), 0, 1);
      } else {
        if (ord[0] > 0)
          set_column(result.num(), 1, r[0] * get_column(result.num(), 0));
        else
          set_column(result.num(), 1, (-r[0]) * get_column(result.num(), 0));
        check(determinant(result.num()) > 0);
      }
    }
    return result;
  }
};

RotMx std_basis(const SpaceGroup &raw_sg, const MatrixGroup &mgc) {
  BasisBuilder bb(mgc.laue_group_type());
  if (bb.n_wtd == 0)
    return RotMx();
  bb.pick_eigen_vectors(raw_sg);
  if (bb.n_wtd == 1)
    return bb.uni_axial_basis();
  return bb.two_axial_basis();
}

bool is_m3b_100_glide(const SpaceGroup &work_sg) {
  check(work_sg.is_centric());
  for (std::size_t i = 1; i < work_sg.n_smx(); ++i) {
    const RotMx &r = work_sg.smx(i).r();
    const int r_type = r.type();
    check(r_type != 0);
    if (iabs(r_type) == 2) {
      const RotMxInfo ri(r);
      if (ri.ev == IntVec3{1, 0, 0}) {
        const std::size_t i_inv = r_type == 2 ? 1 : 0;
        const RtMx smx = work_sg(0, i_inv, i);
        const TrVec wi = smx.t_intrinsic_part();
        return wi[2] % wi.den() != 0;
      }
    }
  }
  fail("is_m3b_100_glide: internal error");
}

RotMx get_adj_rmx(const MatrixGroup &laue, char z) {
  using namespace matrix_group;
  if (laue == code_2_m)
    return RotMx({0, 1, 0, 0, 0, 1, 1, 0, 0}, 1);
  if ((laue == code_4_m || laue == code_4_mmm) && z == 'C')
    return RotMx({1, 1, 0, 1, -1, 0, 0, 0, -1}, 1);
  if ((laue == code_4_m || laue == code_4_mmm) && z == 'F')
    return RotMx({1, 1, 0, -1, 1, 0, 0, 0, 1}, 1);
  if ((laue == code_3b || laue == code_3bm) && z == 'Q')
    return RotMx({-1, 0, 0, 0, -1, 0, 0, 0, 1}, 1);
  if ((laue == code_3bm || laue == code_6_mmm) && z == 'H')
    return RotMx({1, 1, 0, -1, 2, 0, 0, 0, 1}, 1);
  return RotMx(0);
}

// ---- select_generators::standard ----

struct Generators {
  ChangeOfBasisOp z2p_op;
  TrVec z_inv_t, p_inv_t;
  std::size_t n_gen = 0;
  RtMx z_gen[2] = {RtMx(0, 0), RtMx(0, 0)};
  RtMx p_gen[2] = {RtMx(0, 0), RtMx(0, 0)};

  Generators(const SpaceGroup &work_sg, int z2p_r_den, int z2p_t_den,
             const MatrixGroup &point_group) {
    const IntVec3 ev_001{0, 0, 1}, ev_100{1, 0, 0}, ev_110{1, 1, 0},
        ev_m10{-1, 1, 0}, ev_111{1, 1, 1};
    z2p_op = work_sg.z2p_op(z2p_r_den, z2p_t_den);
    z_inv_t = work_sg.inv_t(true);
    p_inv_t = TrVec(0);
    int principal_proper_order = 0;
    switch (point_group.crystal_system()) {
    case CrystalSystem::triclinic:
      break;
    case CrystalSystem::monoclinic:
      z_gen[0] = work_sg.smx(1);
      n_gen = 1;
      break;
    case CrystalSystem::orthorhombic:
      for (std::size_t i = 1; i < work_sg.n_smx(); ++i) {
        const RotMxInfo ri(work_sg.smx(i).r());
        if (ri.ev == ev_001) {
          z_gen[0] = work_sg.smx(i);
          n_gen++;
        } else if (ri.ev == ev_100) {
          z_gen[1] = work_sg.smx(i);
          n_gen++;
        }
      }
      check(n_gen == 2);
      break;
    case CrystalSystem::tetragonal:
      principal_proper_order = 4;
      [[fallthrough]];
    case CrystalSystem::trigonal:
      if (!principal_proper_order)
        principal_proper_order = 3;
      [[fallthrough]];
    case CrystalSystem::hexagonal:
      if (!principal_proper_order)
        principal_proper_order = 6;
      for (std::size_t i = 1; i < work_sg.n_smx(); ++i) {
        const RotMxInfo ri(work_sg.smx(i).r());
        if (iabs(ri.type) == principal_proper_order) {
          if (ri.sense > 0) {
            z_gen[0] = work_sg.smx(i);
            n_gen++;
          }
        } else if (principal_proper_order == 4) {
          if (ri.ev == ev_100) {
            z_gen[1] = work_sg.smx(i);
            n_gen++;
          }
        } else if (principal_proper_order == 3) {
          if (ri.ev == ev_m10) {
            z_gen[1] = work_sg.smx(i);
            n_gen++;
          } else if (ri.ev == ev_110) {
            z_gen[1] = work_sg.smx(i);
            n_gen++;
          }
        } else {
          if (ri.ev == ev_m10) {
            z_gen[1] = work_sg.smx(i);
            n_gen++;
          }
        }
      }
      check(n_gen == 1 || n_gen == 2);
      for (std::size_t i = 0; i < n_gen; ++i)
        check(z_gen[i].is_valid());
      break;
    case CrystalSystem::cubic:
      for (std::size_t i = 1; i < work_sg.n_smx(); ++i) {
        const RotMxInfo ri(work_sg.smx(i).r());
        if (iabs(ri.type) == 4) {
          if (ri.sense > 0 && ri.ev == ev_001) {
            if (!z_gen[0].is_valid())
              n_gen++;
            z_gen[0] = work_sg.smx(i);
          }
        } else if (iabs(ri.type) == 2) {
          if (!z_gen[0].is_valid() && ri.ev == ev_001) {
            z_gen[0] = work_sg.smx(i);
            n_gen++;
          }
        } else if (iabs(ri.type) == 3) {
          if (ri.sense > 0 && ri.ev == ev_111) {
            check(!z_gen[1].is_valid());
            z_gen[1] = work_sg.smx(i);
            n_gen++;
          }
        }
      }
      check(n_gen == 1 || n_gen == 2);
      for (std::size_t i = 0; i < n_gen; ++i)
        check(z_gen[i].is_valid());
      break;
    default:
      fail("select_generators: internal error");
    }
    if (z_inv_t.is_valid())
      for (std::size_t i = 0; i < n_gen; ++i)
        if (determinant(z_gen[i].r().num()) < 0)
          z_gen[i] = z_gen[i].pre_multiply_inv_t(z_inv_t);
    for (std::size_t i = 0; i < n_gen; ++i)
      z_gen[i] = z_gen[i].mod_positive();
  }

  void set_primitive() {
    for (std::size_t i = 0; i < n_gen; ++i)
      p_gen[i] = z2p_op(z_gen[i]).mod_positive();
    if (z_inv_t.is_valid())
      p_inv_t = z2p_op(z_inv_t, -1).mod_positive();
  }
};

// ---- construct_cb_op_t ----

// solve_inhomog_mod_z over flat buffers, reading the Smith normal form's
// diagonal by the same stride cctbx does.
bool solve_inhomog_mod_z(int *m, std::size_t nr, std::size_t nc, int *b,
                         int den, int *x) {
  constexpr std::size_t maxr = 9, maxc = 6;
  check(nr <= maxr && nc <= maxc);
  int p[maxr * maxr];
  int q[maxc * maxc];
  IntRef m_ref{m, nr, nc};
  IntRef p_ref{p, nr, nr};
  IntRef q_ref{q, nc, nc};
  p_ref.set_identity();
  q_ref.set_identity();
  for (;;) {
    row_echelon_form_t(m_ref, p_ref);
    if (m_ref.is_diagonal())
      break;
    m_ref.transpose_in_place();
    row_echelon_form_t(m_ref, q_ref);
    if (m_ref.is_diagonal())
      break;
    m_ref.transpose_in_place();
  }
  q_ref.transpose_square_in_place();
  check(m_ref.is_square());
  const std::size_t nd = m_ref.nr;
  check(nd <= nc);
  int pb[maxr];
  for (std::size_t i = 0; i < nr; ++i) {
    int s = 0;
    for (std::size_t k = 0; k < nr; ++k)
      s += p[i * nr + k] * b[k];
    pb[i] = s;
  }
  for (std::size_t i = nd; i < nr; ++i)
    if (pb[i] % den != 0)
      return false;
  if (x) {
    int xp[maxc];
    for (std::size_t i = 0; i < nc; ++i) {
      xp[i] = 0;
      const int d = m[i * nd + i];
      if (d) {
        check(pb[i] % d == 0);
        xp[i] = pb[i] / d;
      }
    }
    for (std::size_t i = 0; i < nc; ++i) {
      int s = 0;
      for (std::size_t k = 0; k < nc; ++k)
        s += q[i * nc + k] * xp[k];
      x[i] = s;
    }
  }
  return true;
}

TrVec find_origin_shift(const Generators &tab, const Generators &tst,
                        int t_den) {
  RotMx rmi[3];
  TrVec delta_t[3];
  for (std::size_t i = 0; i < tab.n_gen; ++i) {
    check(tst.p_gen[i].r() == tab.p_gen[i].r());
    rmi[i] = tst.p_gen[i].r() - RotMx();
    delta_t[i] = (tst.p_gen[i].t() - tab.p_gen[i].t()).new_denominator(t_den);
  }
  check(tst.p_inv_t.is_valid() == tab.p_inv_t.is_valid());
  std::size_t n_gen = tab.n_gen;
  if (tst.p_inv_t.is_valid()) {
    rmi[n_gen] = RotMx(1, -2);
    delta_t[n_gen] = (tst.p_inv_t - tab.p_inv_t).new_denominator(t_den);
    n_gen++;
  }
  const std::size_t nr_snf = n_gen * 3;
  int snf[9 * 3] = {};
  int v[3 * 3] = {};
  for (std::size_t g = 0; g < n_gen; ++g) {
    for (int i = 0; i < 9; ++i)
      snf[g * 9 + i] = rmi[g][i];
    for (int i = 0; i < 3; ++i)
      v[g * 3 + i] = delta_t[g][i];
  }
  IntVec3 x{};
  if (solve_inhomog_mod_z(snf, nr_snf, 3, v, t_den, x.data())) {
    const TrVec cb_t = tst.z2p_op.c_inv().r() * TrVec(x, t_den);
    return cb_t.new_denominator(t_den);
  }
  return TrVec(0);
}

ChangeOfBasisOp match_generators(int r_den, int t_den,
                                 const SpaceGroup &work_sg,
                                 const MatrixGroup &point_group, char tab_z,
                                 const Generators &tab) {
  if (tab.n_gen == 0 && !tab.z_inv_t.is_valid())
    return ChangeOfBasisOp(r_den, t_den);
  RotMx r_2fold(0), r_3fold(0);
  if (point_group.crystal_system() == CrystalSystem::monoclinic) {
    r_2fold = RotMx({0, 0, 1, 0, -1, 0, 1, 0, 0}, 1);
    r_3fold = RotMx({-1, 0, 1, 0, 1, 0, -1, 0, 0}, 1);
  } else if (point_group.crystal_system() == CrystalSystem::orthorhombic) {
    r_2fold = RotMx({0, 1, 0, 1, 0, 0, 0, 0, -1}, 1);
    r_3fold = RotMx({0, 0, 1, 1, 0, 0, 0, 1, 0}, 1);
  }
  if (r_2fold.is_valid()) {
    const ChangeOfBasisOp cb_op_2fold(
        RtMx(r_2fold.new_denominator(r_den), t_den));
    const ChangeOfBasisOp cb_op_3fold(
        RtMx(r_3fold.new_denominator(r_den), t_den));
    ChangeOfBasisOp cb_op(r_den, t_den);
    for (int i2fold = 0; i2fold < 2; ++i2fold) {
      if (i2fold)
        cb_op = cb_op_2fold;
      for (int i3fold = 0; i3fold < 3; ++i3fold) {
        if (i3fold)
          cb_op.update(cb_op_3fold);
        const SpaceGroup tst_sg = work_sg.change_basis(cb_op);
        const char tst_z = tst_sg.ltr().conventional_centring_type_symbol();
        check(tst_z != '\0' && tst_z != 'Q');
        if (tst_z != tab_z)
          continue;
        Generators tst(tst_sg, r_den, t_den, point_group);
        check(tst.n_gen == tab.n_gen);
        if (tab.n_gen != 2 || (tab.z_gen[0].r()[8] == tst.z_gen[0].r()[8] &&
                               tab.z_gen[1].r()[0] == tst.z_gen[1].r()[0])) {
          tst.set_primitive();
          const TrVec cb_t = find_origin_shift(tab, tst, t_den);
          if (cb_t.is_valid()) {
            cb_op.update(cb_t);
            return cb_op;
          }
        }
      }
    }
  } else {
    Generators tst(work_sg, r_den, t_den, point_group);
    check(tst.n_gen == tab.n_gen);
    tst.set_primitive();
    const TrVec cb_t = find_origin_shift(tab, tst, t_den);
    if (cb_t.is_valid())
      return ChangeOfBasisOp(RtMx(cb_t, r_den));
  }
  return ChangeOfBasisOp(0, 0);
}

void tidy_cb_op_t(const Generators &target, const SpaceGroup &given_sg,
                  const MatrixGroup &point_group, ChangeOfBasisOp &trial) {
  trial = ChangeOfBasisOp(RtMx(trial.c().r(), trial.c().t().den()),
                          RtMx(trial.c_inv().r(), trial.c_inv().t().den()));
  if (given_sg.n_smx() == 1 && !given_sg.is_centric())
    return;
  const SpaceGroup transformed_sg = given_sg.change_basis(trial);
  Generators transformed(transformed_sg, trial.c().r().den(),
                         trial.c().t().den(), point_group);
  transformed.set_primitive();
  const TrVec cb_t =
      find_origin_shift(target, transformed, trial.c().t().den());
  check(cb_t.is_valid());
  trial.update(cb_t.mod_short());
}

// cmp_change_of_basis_mx: true if a ranks before b.
bool cmp_change_of_basis_mx(const RtMx &a, const RtMx &b) {
  const RotMx &ar = a.r(), &br = b.r();
  const TrVec &at = a.t(), &bt = b.t();
  bool ba = ar.is_unit_mx(), bb = br.is_unit_mx();
  if (ba && !bb)
    return true;
  if (!ba && bb)
    return false;
  auto is_diag = [](const IntMat3 &m) {
    return m[1] == 0 && m[2] == 0 && m[3] == 0 && m[5] == 0 && m[6] == 0 &&
           m[7] == 0;
  };
  ba = is_diag(ar.num());
  bb = is_diag(br.num());
  if (ba && !bb)
    return true;
  if (!ba && bb)
    return false;
  ba = at.is_zero();
  bb = bt.is_zero();
  if (ba && !bb)
    return true;
  if (!ba && bb)
    return false;
  std::size_t na = 0, nb = 0;
  for (int i = 0; i < 9; ++i)
    if (ar[i] == 0)
      na++;
  for (int i = 0; i < 9; ++i)
    if (br[i] == 0)
      nb++;
  if (na > nb)
    return true;
  if (na < nb)
    return false;
  na = nb = 0;
  for (int i = 0; i < 9; ++i)
    if (iabs(ar[i]) == ar.den())
      na++;
  for (int i = 0; i < 9; ++i)
    if (iabs(br[i]) == br.den())
      nb++;
  if (na > nb)
    return true;
  if (na < nb)
    return false;
  na = nb = 0;
  for (int i = 0; i < 9; ++i)
    if (ar[i] > 0)
      na++;
  for (int i = 0; i < 9; ++i)
    if (br[i] > 0)
      nb++;
  if (na > nb)
    return true;
  if (na < nb)
    return false;
  if (cmp_i_vec(at.num().data(), bt.num().data(), 3))
    return true;
  if (cmp_i_vec(bt.num().data(), at.num().data(), 3))
    return false;
  return cmp_i_vec(br.num().data(), ar.num().data(), 9);
}

void monoclinic_affine_trial_ranges(const RotMx &cb_mx_r, int &r00, int &r22) {
  r00 = 1;
  r22 = 1;
  for (int i = 0; i < 3; ++i) {
    const int l = lcm_int(cb_mx_r[i], cb_mx_r[6 + i]);
    if (cb_mx_r[i]) {
      const int n = iabs(l / cb_mx_r[i]);
      if (r00 < n)
        r00 = n;
    }
    if (cb_mx_r[i + 6]) {
      const int n = iabs(l / cb_mx_r[6 + i]);
      if (r22 < n)
        r22 = n;
    }
  }
  r00++;
  r22++;
}

bool monoclinic_affine_restrictions(int sg_number, const RotMx &r) {
  switch (sg_number) {
  case 3:
  case 4:
  case 6:
  case 10:
  case 11:
    break;
  case 5:
  case 8:
  case 12:
  case 9:
  case 15:
    if (r[0] % (2 * r.den()) == 0)
      return false;
    if (r[6] % (2 * r.den()) != 0)
      return false;
    if (r[8] % (2 * r.den()) == 0)
      return false;
    break;
  case 7:
  case 13:
  case 14:
    if (r[0] % (2 * r.den()) == 0)
      return false;
    if (r[2] % (2 * r.den()) != 0)
      return false;
    if (r[8] % (2 * r.den()) == 0)
      return false;
    break;
  default:
    fail("check_monoclinic_affine_restrictions: internal error");
  }
  return true;
}

ChangeOfBasisOp find_best_cb_op(const SpaceGroup &given_sg,
                                const MatrixGroup &point_group, int sg_number,
                                const ChangeOfBasisOp &reference_cb_op,
                                const SpaceGroup &target_sg,
                                const Generators &target,
                                const ChangeOfBasisOp &raw_cb_op) {
  const std::vector<RtMx> addl_g =
      normalizer_addl_generators(sg_number, true, true, true);
  const ChangeOfBasisOp cb_op = reference_cb_op.inverse();
  SpaceGroup norm_sg = target_sg;
  std::size_t old_order_p = norm_sg.order_p();
  for (const RtMx &g : addl_g) {
    const RtMx s = cb_op(g.new_denominators(norm_sg.smx()[0]));
    norm_sg.expand_smx(s);
    const std::size_t new_order_p = norm_sg.order_p();
    check(old_order_p * static_cast<std::size_t>(s.r().order()) == new_order_p);
    old_order_p = new_order_p;
  }
  ChangeOfBasisOp best_cb_op = raw_cb_op;
  const bool cmp_inv = determinant(best_cb_op.c().r().num()) <
                       determinant(best_cb_op.c_inv().r().num());
  for (std::size_t i_inv = 0; i_inv < norm_sg.f_inv(); ++i_inv)
    for (std::size_t i_smx = 0; i_smx < norm_sg.n_smx(); ++i_smx) {
      const RtMx s = norm_sg(0, i_inv, i_smx);
      if (determinant(s.r().num()) < 0)
        continue;
      const ChangeOfBasisOp norm_cb_op(s.new_denominators(raw_cb_op.c()));
      ChangeOfBasisOp trial_cb_op = norm_cb_op * raw_cb_op;
      if (sg_number < 3 || sg_number > 15) {
        tidy_cb_op_t(target, given_sg, point_group, trial_cb_op);
        if (!cmp_change_of_basis_mx(best_cb_op.select(cmp_inv),
                                    trial_cb_op.select(cmp_inv)))
          best_cb_op = trial_cb_op;
      } else {
        const int r_den = raw_cb_op.c().r().den();
        const int t_den = raw_cb_op.c().t().den();
        int r00 = 0, r22 = 0;
        monoclinic_affine_trial_ranges(trial_cb_op.c().r(), r00, r22);
        IntMat3 r{}, r_inv{};
        for (r[0] = -r00 * r_den; r[0] <= r00 * r_den; r[0] += r_den)
          for (r[2] = -r22 * r_den; r[2] <= r22 * r_den; r[2] += r_den)
            for (r[6] = -r00 * r_den; r[6] <= r00 * r_den; r[6] += r_den)
              for (r[8] = -r22 * r_den; r[8] <= r22 * r_den; r[8] += r_den) {
                if (!monoclinic_affine_restrictions(sg_number, RotMx(r, r_den)))
                  continue;
                r[4] = r[0] * r[8] - r[2] * r[6];
                if (r[4] != -r_den * r_den && r[4] != r_den * r_den)
                  continue;
                r[4] /= r_den;
                const int f = r[4] / r_den;
                r_inv[0] = f * r[8];
                r_inv[2] = -f * r[2];
                r_inv[4] = r[4];
                r_inv[6] = -f * r[6];
                r_inv[8] = f * r[0];
                ChangeOfBasisOp r_trial =
                    ChangeOfBasisOp(RtMx(RotMx(r, r_den), t_den),
                                    RtMx(RotMx(r_inv, r_den), t_den)) *
                    trial_cb_op;
                tidy_cb_op_t(target, given_sg, point_group, r_trial);
                if (!cmp_change_of_basis_mx(best_cb_op.select(cmp_inv),
                                            r_trial.select(cmp_inv)))
                  best_cb_op = r_trial;
              }
      }
    }
  return best_cb_op;
}

} // namespace

SpaceGroupType::SpaceGroupType(const SpaceGroup &group, bool tidy_cb_op,
                               int r_den, int t_den)
    : group_(group), number_(0), cb_op_(r_den, t_den) {
  thread_local std::unordered_map<Key, std::pair<int, ChangeOfBasisOp>, KeyHash>
      cache;
  Key key{tidy_cb_op, r_den, t_den};
  group_.append_key(key);
  if (const auto it = cache.find(key); it != cache.end()) {
    number_ = it->second.first;
    cb_op_ = it->second.second;
    return;
  }
  determine(tidy_cb_op, r_den, t_den);
  if (cache.size() >= kCacheLimit)
    cache.clear();
  cache.emplace(std::move(key), std::make_pair(number_, cb_op_));
}

void SpaceGroupType::determine(bool tidy_cb_op, int r_den, int t_den) {
  const MatrixGroup point_group = group_.point_group_type();
  const MatrixGroup laue_group = point_group.laue_group_type();
  const CrystalSystem xtal_system = point_group.crystal_system();

  const SpaceGroup z_point_group_sg = group_.build_derived_group(false, false);
  SpaceGroup work_sg = z_point_group_sg;
  char centring_type = '\0';
  int run_away_counter = 0;
  do {
    check(run_away_counter++ < 10);
    ChangeOfBasisOp addl_cb_op = work_sg.z2p_op(r_den, t_den);
    cb_op_.update(addl_cb_op);
    work_sg = z_point_group_sg.change_basis(cb_op_);
    check(work_sg.n_ltr() == 1);

    RotMx basis = std_basis(work_sg, laue_group).new_denominator(r_den);
    addl_cb_op = ChangeOfBasisOp(RtMx(basis, t_den)).inverse();
    cb_op_.update(addl_cb_op);
    work_sg = z_point_group_sg.change_basis(cb_op_);
    centring_type = work_sg.ltr().conventional_centring_type_symbol();

    RotMx adj = get_adj_rmx(laue_group, centring_type);
    if (adj.is_valid()) {
      adj = adj.new_denominator(r_den);
      addl_cb_op = ChangeOfBasisOp(RtMx(adj, t_den));
      cb_op_.update(addl_cb_op);
      work_sg = z_point_group_sg.change_basis(cb_op_);
      centring_type = work_sg.ltr().conventional_centring_type_symbol();
    }
  } while (centring_type == '\0');
  check(centring_type != 'Q');

  work_sg = group_.change_basis(cb_op_);

  if (point_group == matrix_group::code_m3b && centring_type == 'P' &&
      is_m3b_100_glide(work_sg)) {
    const RotMx r_4_001({0, -1, 0, 1, 0, 0, 0, 0, 1}, 1);
    const ChangeOfBasisOp addl_cb_op(
        RtMx(r_4_001.new_denominator(r_den), t_den));
    cb_op_.update(addl_cb_op);
    work_sg = group_.change_basis(cb_op_);
  }

  const MatrixGroup mx_group = get_matrix_group_type(work_sg, point_group);
  check(point_group == mx_group.point_group_type());

  const bool match_centring_type_symbol =
      xtal_system != CrystalSystem::monoclinic &&
      (xtal_system != CrystalSystem::orthorhombic || centring_type == 'I' ||
       centring_type == 'F');

  for (number_ = 1; number_ <= 230; ++number_) {
    const char *hall = hall_symbol_table(number_);
    if (match_centring_type_symbol && centring_type != hall[1])
      continue;
    if ((centring_type == 'P') != (hall[1] == 'P'))
      continue;
    if (matrix_group_code_table(number_) != mx_group)
      continue;
    SpaceGroup tab_sg;
    try {
      tab_sg = SpaceGroup::from_hall(hall, true, false, false, work_sg.t_den());
    } catch (const Error &) {
      fail("space_group_type: internal error parsing a reference symbol");
    }
    if (tab_sg.n_ltr() != work_sg.n_ltr())
      continue;
    Generators tab(tab_sg, r_den, t_den, point_group);
    tab.set_primitive();
    const ChangeOfBasisOp addl_cb_op =
        match_generators(r_den, t_den, work_sg, point_group, hall[1], tab);
    if (addl_cb_op.is_valid()) {
      cb_op_.update(addl_cb_op);
      if (tidy_cb_op)
        cb_op_ =
            find_best_cb_op(group_, point_group, number_,
                            ChangeOfBasisOp(r_den, t_den), tab_sg, tab, cb_op_);
      return;
    }
  }
  fail("space_group_type: no matching reference setting");
}

std::vector<RtMx>
SpaceGroupType::addl_generators_of_euclidean_normalizer(bool k2l,
                                                        bool l2n) const {
  std::vector<RtMx> result =
      normalizer_addl_generators(number_, false, k2l, l2n);
  const ChangeOfBasisOp cb_op_inv = cb_op_.inverse();
  for (RtMx &s : result)
    s = cb_op_inv(s);
  return result;
}

SpaceGroup
SpaceGroupType::expand_addl_generators_of_euclidean_normalizer(bool k2l,
                                                               bool l2n) const {
  SpaceGroup result = group_;
  for (const RtMx &s : addl_generators_of_euclidean_normalizer(k2l, l2n))
    result.expand_smx(s);
  return result;
}

} // namespace mxi::stills::sgtbx

// dials' SymmetryHandler and apply_symmetry against dials' own calls on the
// corpus: the target the handler derives, and the constrained crystal, its
// change of basis and the calls that find nothing.

#include <stdexcept>
#include <string>

#include "check.hh"
#include "real_symmetry.hh"
#include "stills/apply_symmetry.hh"

using namespace mxi::stills;
namespace sg = mxi::stills::sgtbx;

namespace {

SymmetryHandler handler_for(const std::string &dataset) {
  for (const auto &h : real_symmetry::handler_cases)
    if (h.dataset == dataset)
      return SymmetryHandler(UnitCell::from_parameters(h.cell),
                             sg::SpaceGroup::from_hall(h.hall), h.max_delta);
  throw std::invalid_argument("no handler for " + dataset);
}

} // namespace

TEST(stills_symmetry_handler_targets_match_dials) {
  // The handler's reference-setting target is what filter_known_symmetry
  // gets; its derived target must be dials' (target_cases).
  for (const auto &t : real_symmetry::target_cases) {
    const SymmetryHandler h = handler_for(t.dataset);
    check::is_true(h.reference_cell().has_value(), t.dataset + ": cell");
    for (std::size_t i = 0; i < 6; ++i)
      check::is_true(h.reference_cell()->params[i] == t.cell[i],
                     t.dataset + ": reference cell[" + std::to_string(i) + "]");
    const KnownSymmetry k = h.known_symmetry(t.rel, t.ang, t.max_delta);
    check::is_true(k.target_bravais == t.bravais, t.dataset + ": Bravais");
    for (std::size_t i = 0; i < 6; ++i)
      check::is_true(k.target_unit_cell->params[i] == t.target_unit_cell[i],
                     t.dataset + ": target cell[" + std::to_string(i) + "]");
  }
}

TEST(stills_apply_symmetry_matches_dials) {
  int n = 0;
  for (const auto &c : real_symmetry::apply_cases) {
    const std::string name = c.dataset + " call " + std::to_string(n++);
    const SymmetryHandler h = handler_for(c.dataset);
    mxi::Mat3 A;
    A.m = c.A;
    const auto r = h.apply_symmetry(A, UnitCell::from_parameters(c.cell));
    check::is_true(r.has_value() == c.found, name + ": outcome");
    if (!r)
      continue;
    const mxi::Mat3 out = r->model.A();
    for (std::size_t i = 0; i < 9; ++i)
      check::is_true(out.m[i] == c.out_A[i],
                     name + ": A[" + std::to_string(i) + "]");
    for (std::size_t i = 0; i < 6; ++i)
      check::is_true(r->model.cell.params[i] == c.out_cell[i],
                     name + ": cell[" + std::to_string(i) + "]");
    const sg::RtMx &cb = r->cb_op.c();
    for (std::size_t i = 0; i < 9; ++i)
      check::equal(cb.r()[i], c.cb_op_c[i], name + ": cb_op r");
    check::equal(cb.r().den(), c.cb_op_c[9], name + ": cb_op r_den");
    for (std::size_t i = 0; i < 3; ++i)
      check::equal(cb.t()[i], c.cb_op_c[10 + i], name + ": cb_op t");
    check::equal(cb.t().den(), c.cb_op_c[13], name + ": cb_op t_den");
  }
}

TEST(stills_constrain_orientation_reaches_the_system) {
  // A slightly sheared tetragonal direct matrix comes out with a = b and the
  // angles at 90 degrees to well within the fit's tolerance.
  mxi::Mat3 direct;
  direct.m = {78.0, 0.01, 0.02, -0.03, 78.1, 0.0, 0.05, 0.0, 37.0};
  const mxi::Mat3 d = constrain_orientation(direct, "tetragonal");
  const mxi::Vec3 a = d.row(0), b = d.row(1), c = d.row(2);
  check::close(a.dot(a), b.dot(b), 1e-6 * a.dot(a), "a = b");
  check::close(a.dot(b), 0, 1e-6 * a.dot(a), "gamma 90");
  check::close(a.dot(c), 0, 1e-6 * a.dot(a), "beta 90");
  check::close(b.dot(c), 0, 1e-6 * a.dot(a), "alpha 90");
}

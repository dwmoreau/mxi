// Python bindings for the space-group machinery and dials' known-symmetry
// filter. Groups cross as tuples in the order the group holds them, so a
// caller can compare that order with cctbx's:
//   group  = (ltr, f_inv, inv_t, smx)
//   ltr    = ((n0, n1, n2, den), ...)
//   inv_t  = (n0, n1, n2, den), or None when f_inv is 1
//   smx    = ((r0..r8, r_den, t0, t1, t2, t_den), ...)
//   cb_op  = (c, c_inv), each in the smx form.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <optional>
#include <string>

#include "stills/apply_symmetry.hh"
#include "stills/known_symmetry.hh"
#include "stills/sgtbx.hh"

namespace py = pybind11;
using namespace mxi::stills;
namespace sg = mxi::stills::sgtbx;

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

py::tuple rt_tuple(const sg::RtMx &s) {
  py::list l;
  for (int i = 0; i < 9; ++i)
    l.append(s.r()[static_cast<std::size_t>(i)]);
  l.append(s.r().den());
  for (int i = 0; i < 3; ++i)
    l.append(s.t()[static_cast<std::size_t>(i)]);
  l.append(s.t().den());
  return py::tuple(l);
}

sg::RtMx rt_from(const py::sequence &v) {
  if (py::len(v) != 14)
    throw py::value_error("expected (r0..r8, r_den, t0, t1, t2, t_den)");
  sg::IntMat3 r{};
  for (int i = 0; i < 9; ++i)
    r[static_cast<std::size_t>(i)] = v[static_cast<std::size_t>(i)].cast<int>();
  const sg::TrVec t(v[10].cast<int>(), v[11].cast<int>(), v[12].cast<int>(),
                    v[13].cast<int>());
  return sg::RtMx(sg::RotMx(r, v[9].cast<int>()), t);
}

py::tuple tr_tuple(const sg::TrVec &t) {
  return py::make_tuple(t[0], t[1], t[2], t.den());
}

py::tuple group_tuple(const sg::SpaceGroup &g) {
  py::list ltr, smx;
  for (const auto &t : g.ltr().elems())
    ltr.append(tr_tuple(t));
  for (const auto &s : g.smx())
    smx.append(rt_tuple(s));
  py::object inv = py::none();
  if (g.is_centric())
    inv = tr_tuple(g.inv_t());
  return py::make_tuple(py::tuple(ltr), g.f_inv(), inv, py::tuple(smx));
}

// A group rebuilt from its lattice translations, inversion and operations by
// expansion. The set is the caller's; the order need not be.
sg::SpaceGroup group_from(const py::tuple &g) {
  sg::SpaceGroup out;
  for (const auto &t : g[0].cast<py::sequence>()) {
    const auto v = t.cast<py::sequence>();
    out.expand_ltr(sg::TrVec(v[0].cast<int>(), v[1].cast<int>(),
                             v[2].cast<int>(), v[3].cast<int>()));
  }
  if (!g[2].is_none()) {
    const auto v = g[2].cast<py::sequence>();
    out.expand_inv(sg::TrVec(v[0].cast<int>(), v[1].cast<int>(),
                             v[2].cast<int>(), v[3].cast<int>()));
  }
  for (const auto &s : g[3].cast<py::sequence>())
    out.expand_smx(rt_from(s.cast<py::sequence>()));
  return out;
}

py::tuple cb_tuple(const sg::ChangeOfBasisOp &cb) {
  return py::make_tuple(rt_tuple(cb.c()), rt_tuple(cb.c_inv()));
}

sg::ChangeOfBasisOp cb_from(const py::tuple &cb) {
  return sg::ChangeOfBasisOp(rt_from(cb[0].cast<py::sequence>()),
                             rt_from(cb[1].cast<py::sequence>()));
}

UnitCell cell_from(const py::sequence &p) {
  if (py::len(p) != 6)
    throw py::value_error("expected six cell parameters");
  std::array<double, 6> params{};
  for (std::size_t i = 0; i < 6; ++i)
    params[i] = p[i].cast<double>();
  return UnitCell::from_parameters(params);
}

py::tuple cell_tuple(const UnitCell &c) {
  return py::make_tuple(c.params[0], c.params[1], c.params[2], c.params[3],
                        c.params[4], c.params[5]);
}

sg::SpaceGroup hall_group(const std::string &hall,
                          const std::optional<py::tuple> &cb) {
  sg::SpaceGroup g = sg::SpaceGroup::from_hall(hall);
  if (cb)
    g = g.change_basis(cb_from(*cb));
  return g;
}

mxi::Mat3 mat3_from(const Array &a, const char *name) {
  if (a.ndim() != 2 || a.shape(0) != 3 || a.shape(1) != 3)
    throw py::value_error(std::string(name) + ": expected shape (3, 3)");
  mxi::Mat3 m;
  for (py::ssize_t r = 0; r < 3; ++r)
    for (py::ssize_t c = 0; c < 3; ++c)
      m(static_cast<std::size_t>(r), static_cast<std::size_t>(c)) = a.at(r, c);
  return m;
}

Array mat3_array(const mxi::Mat3 &m) {
  Array out({py::ssize_t{3}, py::ssize_t{3}});
  auto o = out.mutable_unchecked<2>();
  for (py::ssize_t r = 0; r < 3; ++r)
    for (py::ssize_t c = 0; c < 3; ++c)
      o(r, c) = m(static_cast<std::size_t>(r), static_cast<std::size_t>(c));
  return out;
}

py::dict matching_dict(const MatchingSymmetry &m) {
  py::dict d;
  d["max_angular_difference"] = m.max_angular_difference;
  d["subsym_cell"] = cell_tuple(m.subsym.cell);
  d["ref_subsym_cell"] = cell_tuple(m.ref_subsym.cell);
  d["best_subsym_cell"] = cell_tuple(m.best_subsym.cell);
  d["subsym_group"] = group_tuple(m.subsym.group);
  d["best_subsym_group"] = group_tuple(m.best_subsym.group);
  return d;
}

} // namespace

void bind_known_symmetry(py::module_ &m) {
  py::register_exception<sg::Error>(m, "SgtbxError", PyExc_RuntimeError);
  py::register_exception<IncompatibleUnitCell>(m, "IncompatibleUnitCell",
                                               PyExc_AssertionError);

  m.def(
      "sgtbx_group",
      [](const std::string &hall, std::optional<py::tuple> cb, bool tidy) {
        sg::SpaceGroup g = hall_group(hall, cb);
        if (tidy)
          g.make_tidy();
        return group_tuple(g);
      },
      py::arg("hall"), py::arg("cb_op") = py::none(), py::arg("tidy") = false,
      "cctbx's space_group(hall).change_basis(cb_op), optionally tidied, as "
      "a group tuple in its stored order.");

  m.def(
      "sgtbx_type",
      [](const std::string &hall, std::optional<py::tuple> cb,
         bool tidy_cb_op) {
        const sg::SpaceGroupType t(hall_group(hall, cb), tidy_cb_op);
        return py::make_tuple(t.number(), cb_tuple(t.cb_op()));
      },
      py::arg("hall"), py::arg("cb_op") = py::none(),
      py::arg("tidy_cb_op") = true,
      "cctbx's space_group_type of space_group(hall).change_basis(cb_op): "
      "the number and the change of basis to the reference setting.");

  m.def(
      "sgtbx_type_of_group",
      [](const py::tuple &g, bool tidy_cb_op) {
        const sg::SpaceGroupType t(group_from(g), tidy_cb_op);
        return py::make_tuple(t.number(), cb_tuple(t.cb_op()));
      },
      py::arg("group"), py::arg("tidy_cb_op") = true);

  m.def(
      "metric_supergroup",
      [](const std::string &hall, std::optional<py::tuple> cb) {
        return group_tuple(metric_supergroup(hall_group(hall, cb)));
      },
      py::arg("hall"), py::arg("cb_op") = py::none(),
      "dials' metric_supergroup of space_group(hall).change_basis(cb_op).");

  m.def(
      "bravais_symbol",
      [](const std::string &hall, std::optional<py::tuple> cb) {
        return sg::bravais_symbol(hall_group(hall, cb));
      },
      py::arg("hall"), py::arg("cb_op") = py::none());

  m.def(
      "minimum_cell",
      [](const py::sequence &cell) {
        const CrystalSymmetry cs =
            CrystalSymmetry::make(cell_from(cell), sg::SpaceGroup());
        const sg::ChangeOfBasisOp cb = cs.change_of_basis_op_to_minimum_cell();
        return py::make_tuple(cb_tuple(cb),
                              cell_tuple(cs.change_basis(cb).cell));
      },
      py::arg("cell"),
      "dials' minimum cell of a P1 symmetry: the change of basis and the "
      "cell.");

  m.def(
      "lattice_symmetry_group",
      [](const py::sequence &cell, double max_delta, bool enforce) {
        return group_tuple(
            lattice_symmetry_group(cell_from(cell), max_delta, enforce));
      },
      py::arg("reduced_cell"), py::arg("max_delta") = 3.0,
      py::arg("enforce_max_delta_for_generated_two_folds") = true);

  m.def(
      "acentric_subgroups",
      [](const py::sequence &reduced_cell, double max_delta,
         const std::string &target_bravais) {
        const sg::SpaceGroup lg =
            lattice_symmetry_group(cell_from(reduced_cell), max_delta, true);
        const AcentricSubgroups a = calc_acentric_subgroups(lg, target_bravais);
        py::list out;
        for (std::size_t i = 0; i < a.subgroups.size(); ++i)
          out.append(py::make_tuple(group_tuple(a.subgroups[i]),
                                    cb_tuple(a.cb_ops[i]),
                                    group_tuple(a.supergroups[i])));
        return out;
      },
      py::arg("reduced_cell"), py::arg("max_delta"), py::arg("target_bravais"),
      "dials' calc_acentric_subgroups of the lattice group of a reduced cell: "
      "(subgroup, cb_op to reference, metric supergroup) in dials' order.");

  m.def(
      "find_matching_symmetry",
      [](const py::sequence &cell, double max_delta,
         const std::string &target_bravais) -> py::object {
        const auto r =
            find_matching_symmetry(cell_from(cell), max_delta, target_bravais);
        if (!r)
          return py::none();
        return matching_dict(*r);
      },
      py::arg("unit_cell"), py::arg("max_delta"), py::arg("target_bravais"),
      "dials' find_matching_symmetry(unit_cell, None, max_delta, "
      "target_bravais_str=...): None, or the best subgroup's "
      "max_angular_difference and cells.");

  m.def(
      "iotbx_converter",
      [](const py::sequence &cell, double max_delta,
         const std::string &target_bravais) {
        py::list out;
        for (const auto &item :
             iotbx_converter(cell_from(cell), max_delta, target_bravais)) {
          py::dict d;
          d["max_angular_difference"] = item.max_angular_difference;
          d["bravais"] = item.bravais;
          d["system"] = item.system;
          d["number"] = item.number;
          d["cb_op_inp_best"] = cb_tuple(item.cb_op_inp_best);
          d["subsym_cell"] = cell_tuple(item.subsym.cell);
          d["ref_subsym_cell"] = cell_tuple(item.ref_subsym.cell);
          d["best_subsym_cell"] = cell_tuple(item.best_subsym.cell);
          out.append(d);
        }
        return out;
      },
      py::arg("unit_cell"), py::arg("max_delta"), py::arg("target_bravais"),
      "rstbx's iotbx_converter(unit_cell, max_delta, target_bravais=...): the "
      "metric subgroups in its sorted order.");

  py::class_<SymmetryHandler>(m, "SymmetryHandler")
      .def(py::init([](std::optional<py::sequence> cell,
                       std::optional<std::string> hall, double max_delta) {
             std::optional<UnitCell> c;
             if (cell)
               c = cell_from(*cell);
             std::optional<sg::SpaceGroup> g;
             if (hall)
               g = sg::SpaceGroup::from_hall(*hall);
             return SymmetryHandler(c, g, max_delta);
           }),
           py::arg("unit_cell"), py::arg("hall"), py::arg("max_delta") = 5.0,
           "dials' SymmetryHandler(unit_cell, space_group, max_delta), the "
           "space group given as a Hall symbol.")
      .def_property_readonly("reference_cell",
                             [](const SymmetryHandler &h) -> py::object {
                               if (!h.reference_cell())
                                 return py::none();
                               return cell_tuple(*h.reference_cell());
                             })
      .def_property_readonly("reference_group",
                             [](const SymmetryHandler &h) -> py::object {
                               if (!h.reference_group())
                                 return py::none();
                               return group_tuple(*h.reference_group());
                             })
      .def_property_readonly("target_bravais", &SymmetryHandler::target_bravais)
      .def("known_symmetry", &SymmetryHandler::known_symmetry,
           py::arg("relative_length_tolerance") = 0.1,
           py::arg("absolute_angle_tolerance") = 5.0,
           py::arg("max_delta") = 5.0)
      .def(
          "apply_symmetry",
          [](const SymmetryHandler &h, const Array &A,
             const py::sequence &cell) -> py::object {
            const auto r = h.apply_symmetry(mat3_from(A, "A"), cell_from(cell));
            if (!r)
              return py::none();
            py::dict d;
            d["A"] = mat3_array(r->model.A());
            d["cell"] = cell_tuple(r->model.cell);
            d["space_group"] = group_tuple(r->space_group);
            d["cb_op"] = cb_tuple(r->cb_op);
            return d;
          },
          py::arg("A"), py::arg("unit_cell"),
          "dials' apply_symmetry on a crystal given by its A and cell: None, "
          "or the constrained crystal (A, cell) in the input setting, its "
          "space group, and the change of basis to the user's target.");

  py::class_<KnownSymmetry>(m, "KnownSymmetry")
      .def(
          py::init([](std::optional<py::sequence> cell, const std::string &hall,
                      double rel, double ang, double max_delta) {
            std::optional<UnitCell> c;
            if (cell)
              c = cell_from(*cell);
            return KnownSymmetry::make(c, sg::SpaceGroup::from_hall(hall), rel,
                                       ang, max_delta);
          }),
          py::arg("unit_cell"), py::arg("hall"),
          py::arg("relative_length_tolerance") = 0.1,
          py::arg("absolute_angle_tolerance") = 5.0, py::arg("max_delta") = 5.0,
          "The target dials' filter_known_symmetry compares with: the target "
          "symmetry's cell as dials holds it (or None) and its group as a "
          "Hall symbol.")
      .def_property_readonly("target_unit_cell",
                             [](const KnownSymmetry &k) -> py::object {
                               if (!k.target_unit_cell)
                                 return py::none();
                               return cell_tuple(*k.target_unit_cell);
                             })
      .def_property_readonly("volume_window",
                             [](const KnownSymmetry &k) -> py::object {
                               if (!k.volume_window)
                                 return py::none();
                               return py::make_tuple(k.volume_window->first,
                                                     k.volume_window->second);
                             })
      .def_readonly("target_bravais", &KnownSymmetry::target_bravais)
      .def(
          "filter",
          [](const KnownSymmetry &k, const Array &cells) {
            const py::ssize_t n = cells.ndim() == 2 ? cells.shape(0) : -1;
            if (n < 0 || cells.shape(1) != 6)
              throw py::value_error("cells: expected shape (n, 6)");
            py::array_t<int> outcome(n);
            Array max_ang(n);
            Array best({n, py::ssize_t{6}});
            auto o = outcome.mutable_unchecked<1>();
            auto a = max_ang.mutable_unchecked<1>();
            auto b = best.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < n; ++i) {
              std::array<double, 6> p{};
              for (py::ssize_t j = 0; j < 6; ++j)
                p[static_cast<std::size_t>(j)] = cells.at(i, j);
              const auto r =
                  filter_known_symmetry(UnitCell::from_parameters(p), k);
              o(i) = static_cast<int>(r.outcome);
              a(i) = r.best ? r.best->max_angular_difference : std::nan("");
              for (py::ssize_t j = 0; j < 6; ++j)
                b(i, j) = r.best ? r.best->best_subsym.cell
                                       .params[static_cast<std::size_t>(j)]
                                 : std::nan("");
            }
            return py::make_tuple(outcome, max_ang, best);
          },
          py::arg("cells"),
          "Each cell through filter_known_symmetry: the outcome (0 accepted, "
          "1 outside the volume window, 2 no subgroup of the target Bravais "
          "type, 3 best subgroup not similar to the target), and the best "
          "subgroup's max_angular_difference and best cell (NaN without "
          "one).");
}

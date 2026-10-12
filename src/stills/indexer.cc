// Derived from dials.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.

#include "stills/indexer.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <utility>

#include "stills/candidates.hh"
#include "stills/crystal.hh"
#include "stills/method_list.hh"
#include "stills/outliers.hh"
#include "stills/predict.hh"

namespace mxi::stills {

namespace {

constexpr std::uint64_t kPredicted = 1u << 0;
constexpr std::uint64_t kIndexed = 1u << 2;
constexpr std::uint64_t kUsedInRefinement = 1u << 3;

// An error that ends one indexing method, with the message dials' exception
// carries.
struct MethodFailure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// An error that ends one candidate's chain.
struct ChainFailure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

CellConstraint constraint_of(const sgtbx::SpaceGroup &group) {
  std::vector<std::array<int, 9>> rotations;
  for (const sgtbx::RtMx &s : group.smx()) {
    const sgtbx::RotMx &r = s.r();
    std::array<int, 9> m{};
    for (std::size_t i = 0; i < 9; ++i)
      m[i] = r[i] / r.den();
    rotations.push_back(m);
  }
  return CellConstraint::from_rotations(rotations);
}

// cctbx's change_of_basis_op.apply(miller_index): h times the inverse's
// rotation part, which must divide exactly by its denominator.
std::array<int, 3> apply_cb_op(const sgtbx::ChangeOfBasisOp &cb,
                               const std::array<int, 3> &h) {
  const sgtbx::RotMx &r = cb.c_inv().r();
  std::array<int, 3> out{};
  for (std::size_t j = 0; j < 3; ++j) {
    const int v = h[0] * r[j] + h[1] * r[3 + j] + h[2] * r[6 + j];
    if (v % r.den() != 0)
      throw MethodFailure("Change of basis yields non-integral Miller index.");
    out[j] = v / r.den();
  }
  return out;
}

Mat3 cb_op_matrix(const sgtbx::ChangeOfBasisOp &cb) {
  return cb.c().r().transpose().as_double();
}

std::vector<std::string> operators_of(const sgtbx::SpaceGroup &group) {
  std::vector<std::string> out;
  for (std::size_t i = 0; i < group.order_z(); ++i)
    out.push_back(group(i).as_xyz());
  return out;
}

std::string format_double(const char *format, double v) {
  char buffer[128];
  std::snprintf(buffer, sizeof buffer, format, v);
  return buffer;
}

} // namespace

struct StillsIndexer::Frame {
  Vec3 s0;
  double wavelength = 0;
  const StrongSpots *strong = nullptr;
  MappedSpots mapped;
  std::vector<double> d; // 1 / |rlp|
  std::array<double, 2> pixel_size_0{};
  double d_min = 0;
};

struct StillsIndexer::ChainResult {
  CrystalModel crystal;
  Mat3 raw_B; // the refinement's B, before dxtbx's set_B
  NaveResult nave;
  std::vector<std::size_t> rows;
  std::vector<std::array<int, 3>> miller;
  Predictions predictions;
};

StillsIndexer::StillsIndexer(StillsDetector detector, IndexerSettings settings)
    : detector_(std::move(detector)), settings_(std::move(settings)),
      user_group_(sgtbx::SpaceGroup::from_hall(settings_.space_group_hall)),
      symmetry_(UnitCell::from_parameters(settings_.unit_cell), user_group_,
                settings_.max_delta),
      known_symmetry_(symmetry_.known_symmetry(
          settings_.relative_length_tolerance,
          settings_.absolute_angle_tolerance, settings_.max_delta)) {
  if (detector_.size() == 0)
    throw std::invalid_argument("StillsIndexer: a detector with no panels");
  const auto &p = symmetry_.primitive_cell()->params;
  max_cell_ = settings_.max_cell_multiplier * std::max({p[0], p[1], p[2]});
  for (const std::string &m : settings_.methods) {
    if (m == "real_space_grid_search") {
      if (!grid_search_)
        grid_search_ = std::make_unique<RealSpaceGridSearch>(
            p, settings_.rsgs_characteristic_grid);
    } else if (m != "fft1d") {
      throw std::invalid_argument("StillsIndexer: unknown method " + m);
    }
  }
}

std::vector<Vec3> StillsIndexer::search(const std::string &method,
                                        const std::vector<Vec3> &rlp) const {
  if (method == "fft1d") {
    try {
      return fft1d_basis_vectors(rlp, max_cell_,
                                 settings_.fft1d_characteristic_grid, nullptr,
                                 settings_.fft1d_grid);
    } catch (const TooFewSpots &) {
      // rstbx asserts, with no message.
      throw MethodFailure("");
    }
  }
  return grid_search_->find_basis_vectors(rlp, settings_.rsgs_max_vectors);
}

StillsIndexer::ChainResult
StillsIndexer::refine_chain(const Frame &frame, const CrystalModel &start,
                            const sgtbx::SpaceGroup &group,
                            const std::vector<std::size_t> &rows,
                            const std::vector<std::array<int, 3>> &miller,
                            std::vector<StageTrace> *stages) const {
  const CellConstraint constraint = constraint_of(group);
  const MappedSpots &m = frame.mapped;
  const Spots &spots = frame.strong->spots;

  auto spots_of = [&](const std::vector<std::size_t> &r,
                      const std::vector<std::array<int, 3>> &h) {
    IndexedSpots s;
    s.miller = h;
    for (std::size_t i : r) {
      s.panel.push_back(spots.panel[i]);
      s.x_mm.push_back(m.x_mm[i]);
      s.y_mm.push_back(m.y_mm[i]);
      s.var_x_mm.push_back(m.var_x_mm[i]);
      s.var_y_mm.push_back(m.var_y_mm[i]);
      s.flags.push_back(frame.strong->flags.empty() ? 0
                                                    : frame.strong->flags[i]);
    }
    return s;
  };
  Mat3 last_raw_B;
  auto refine = [&](const CrystalModel &from, const IndexedSpots &s) {
    const RefineResult r = refine_crystal(detector_, frame.s0, from.U, from.B,
                                          constraint, s, settings_.refine);
    stages->push_back({"refine", s.size(),
                       static_cast<std::size_t>(r.n_iterations), r.rmsd[0],
                       r.rmsd[1]});
    if (failed(r.status))
      throw ChainFailure("refinement: " + r.message);
    last_raw_B = r.B;
    return CrystalModel::from_U_B(r.U, r.B);
  };
  auto predict = [&](const CrystalModel &c, const IndexedSpots &s) {
    Predictions p =
        predict_spots(detector_, frame.s0, c.U, c.B, s.miller, s.panel);
    for (PredictStatus st : p.status)
      if (st != PredictStatus::ok)
        throw ChainFailure("prediction: a reflection cannot be predicted");
    return p;
  };
  // identify_outliers: refine, predict, Sauter-Poon on positions in mm formed
  // as pixels times panel 0's pixel size.
  auto outliers = [&](const CrystalModel &from,
                      const std::vector<std::size_t> &r,
                      const IndexedSpots &s) {
    const CrystalModel refined = refine(from, s);
    const Predictions p = predict(refined, s);
    const std::size_t n = r.size();
    std::vector<double> xo(n), yo(n), xc(n), yc(n);
    for (std::size_t k = 0; k < n; ++k) {
      const auto [px, py] =
          detector_.panel(s.panel[k]).mm_to_px(p.x_mm[k], p.y_mm[k]);
      xo[k] = spots.x_px[r[k]] * frame.pixel_size_0[0];
      yo[k] = spots.y_px[r[k]] * frame.pixel_size_0[1];
      xc[k] = px * frame.pixel_size_0[0];
      yc[k] = py * frame.pixel_size_0[1];
    }
    const SauterPoonResult sp = sauter_poon(xo, yo, xc, yc);
    std::size_t n_kept = 0;
    for (bool k : sp.kept)
      n_kept += k;
    stages->push_back({"sauter_poon", n, n_kept, sp.sigma, 0});
    if (sp.status == SauterPoonStatus::too_few)
      throw ChainFailure("No reflections left after removing outliers!");
    if (sp.status != SauterPoonStatus::ok)
      throw ChainFailure(std::string("Sauter-Poon: ") + to_string(sp.status));
    return sp.kept;
  };
  auto nave = [&](const CrystalModel &c, const IndexedSpots &s) {
    const Predictions p = predict(c, s);
    const std::vector<double> d = cell_d_spacings(c.B, s.miller);
    NaveResult r =
        nave_parameters(p.delpsi, d, two_theta_deg(d, frame.wavelength),
                        cell_volume(c.B), frame.wavelength, settings_.nave);
    std::size_t n_kept = 0;
    for (bool k : r.accepted)
      n_kept += k;
    stages->push_back({"nave", s.size(), n_kept, r.alpha, r.eta});
    if (r.status == NaveStatus::ewald_volume)
      throw ChainFailure(
          format_double("Ewald proximity volume too high, %f", r.ewald_volume));
    if (r.status == NaveStatus::diverged)
      throw ChainFailure("Model has diverged, cannot continue");
    if (r.status == NaveStatus::outside_domain)
      throw ChainFailure("max likelihood outside of math domain");
    if (failed(r.status))
      throw ChainFailure(std::string("Nave: ") + to_string(r.status));
    return r;
  };
  auto select = [](const std::vector<std::size_t> &r,
                   const std::vector<std::array<int, 3>> &h,
                   const std::vector<bool> &keep) {
    std::pair<std::vector<std::size_t>, std::vector<std::array<int, 3>>> out;
    for (std::size_t k = 0; k < r.size(); ++k)
      if (keep[k]) {
        out.first.push_back(r[k]);
        out.second.push_back(h[k]);
      }
    return out;
  };

  const std::vector<bool> kept_1 =
      outliers(start, rows, spots_of(rows, miller));
  const auto [rows_1, miller_1] = select(rows, miller, kept_1);
  const IndexedSpots spots_1 = spots_of(rows_1, miller_1);

  const CrystalModel refined_b = refine(start, spots_1);
  const std::vector<bool> kept_2 = outliers(refined_b, rows_1, spots_1);
  const NaveResult nave_0 = nave(refined_b, spots_1);
  std::vector<bool> keep(rows_1.size());
  for (std::size_t k = 0; k < keep.size(); ++k)
    keep[k] = kept_2[k] && nave_0.accepted[k];
  const auto [rows_2, miller_2] = select(rows_1, miller_1, keep);
  const IndexedSpots spots_2 = spots_of(rows_2, miller_2);

  ChainResult out;
  out.crystal = refine(refined_b, spots_2);
  out.raw_B = last_raw_B;
  out.nave = nave(out.crystal, spots_2);
  out.rows = rows_2;
  out.miller = miller_2;
  out.predictions = predict(out.crystal, spots_2);
  return out;
}

void StillsIndexer::index_with(const std::string &method, Frame &frame,
                               IndexedStill &out, MethodTrace &trace) const {
  const auto t_search = std::chrono::steady_clock::now();
  const std::vector<Vec3> &rlp = frame.mapped.rlp;
  const std::size_t n_all = rlp.size();

  // BasisVectorSearch.find_candidate_basis_vectors.
  frame.d_min = settings_.d_min_start;
  std::vector<std::size_t> searched;
  for (std::size_t i = 0; i < n_all; ++i)
    if (frame.d_min <= 0 || frame.d[i] > frame.d_min)
      searched.push_back(i);
  std::vector<Vec3> rlp_searched;
  for (std::size_t i : searched)
    rlp_searched.push_back(rlp[i]);
  const std::vector<Vec3> basis = search(method, rlp_searched);
  trace.n_basis_vectors = basis.size();
  if (frame.d_min <= 0) {
    if (searched.empty())
      throw MethodFailure("flex.min: empty array");
    frame.d_min = frame.d[searched[0]];
    for (std::size_t i : searched)
      frame.d_min = std::min(frame.d_min, frame.d[i]);
  }
  trace.d_min = frame.d_min;
  trace.seconds_search = seconds_since(t_search);

  // choose_best_orientation_matrix over the candidates, made lazily.
  const auto t_candidates = std::chrono::steady_clock::now();
  const bool known_symmetry_path =
      settings_.refine_candidates_with_known_symmetry;
  const std::vector<int> unindexed(n_all, -1);
  struct Candidate {
    CrystalModel crystal;
    sgtbx::SpaceGroup group;
    NaveResult nave;
    double rmsd = 0;
  };
  std::vector<Candidate> candidates;
  CandidateStream stream(basis);
  for (std::size_t icm = 0;; ++icm) {
    std::optional<CrystalModel> cm;
    while ((cm = stream.next()))
      if (filter_known_symmetry(cm->cell, known_symmetry_).outcome ==
          KnownSymmetryOutcome::accepted)
        break;
    if (!cm || icm >= settings_.max_refine)
      break;

    CandidateTrace ct;
    const Assignment a = assign_indices_global(
        rlp, unindexed, {cm->A()}, frame.d_min, settings_.hkl_tolerance);
    std::vector<std::size_t> rows;
    std::vector<std::array<int, 3>> miller;
    for (std::size_t i = 0; i < n_all; ++i)
      if (a.id[i] >= 0) {
        rows.push_back(i);
        miller.push_back(a.miller_index[i]);
      }
    ct.n_indexed = rows.size();

    CrystalModel crystal = *cm;
    sgtbx::SpaceGroup group; // P1
    if (known_symmetry_path) {
      const std::optional<AppliedSymmetry> applied =
          symmetry_.apply_symmetry(crystal);
      if (!applied) {
        ct.status = "cannot convert to target symmetry";
        trace.candidates.push_back(ct);
        continue;
      }
      crystal = applied->model.change_basis(cb_op_matrix(applied->cb_op));
      group = applied->space_group.change_basis(applied->cb_op);
      if (!applied->cb_op.is_identity_op())
        for (auto &h : miller)
          h = apply_cb_op(applied->cb_op, h);
    }

    try {
      ChainResult chain =
          refine_chain(frame, crystal, group, rows, miller, &ct.stages);
      ct.n_indexed = chain.rows.size();
      ct.cell = chain.crystal.cell.params;
      ct.half_mosaicity_deg = chain.nave.half_mosaicity_deg;
      ct.domain_size_ang = chain.nave.domain_size_ang;
      if (!known_symmetry_path && !symmetry_.apply_symmetry(chain.crystal)) {
        ct.status = "P1 refinement yielded model diverged from target";
        trace.candidates.push_back(ct);
        continue;
      }
      // calc_2D_rmsd_and_displacements, in pixels.
      double sum = 0;
      const Spots &spots = frame.strong->spots;
      for (std::size_t k = 0; k < chain.rows.size(); ++k) {
        const std::size_t i = chain.rows[k];
        const auto [px, py] =
            detector_.panel(spots.panel[i])
                .mm_to_px(chain.predictions.x_mm[k], chain.predictions.y_mm[k]);
        const double dx = spots.x_px[i] - px, dy = spots.y_px[i] - py;
        sum += dx * dx + dy * dy;
      }
      ct.rmsd = std::sqrt(sum * (1. / static_cast<double>(chain.rows.size())));
      ct.status = "ok";
      candidates.push_back({chain.crystal, group, chain.nave, ct.rmsd});
    } catch (const std::exception &e) {
      ct.status = e.what();
    }
    trace.candidates.push_back(ct);
  }
  trace.seconds_candidates = seconds_since(t_candidates);

  if (candidates.empty())
    throw MethodFailure("No suitable indexing solution found");
  std::size_t best = 0;
  for (std::size_t c = 1; c < candidates.size(); ++c)
    if (candidates[c].rmsd < candidates[best].rmsd)
      best = c;
  {
    std::size_t k = 0;
    for (std::size_t c = 0; c < trace.candidates.size(); ++c)
      if (trace.candidates[c].status == "ok" && k++ == best)
        trace.best = static_cast<int>(c);
  }
  if (candidates[best].rmsd > settings_.xy_rmsd_threshold)
    throw MethodFailure(
        format_double("RMSD too high, %f", candidates[best].rmsd));

  // StillsIndexer.index: re-index every spot with the chosen crystal.
  const auto t_refine = std::chrono::steady_clock::now();
  CrystalModel crystal = candidates[best].crystal;
  sgtbx::SpaceGroup group = candidates[best].group;
  trace.chosen_A = crystal.A();
  const Assignment a = assign_indices_global(
      rlp, unindexed, {crystal.A()}, frame.d_min, settings_.hkl_tolerance);
  std::vector<std::size_t> rows;
  std::vector<std::array<int, 3>> miller;
  for (std::size_t i = 0; i < n_all; ++i)
    if (a.id[i] >= 0) {
      rows.push_back(i);
      miller.push_back(a.miller_index[i]);
    }
  if (!known_symmetry_path) {
    const std::optional<AppliedSymmetry> applied =
        symmetry_.apply_symmetry(crystal);
    if (!applied)
      throw MethodFailure("'NoneType' object has no attribute 'change_basis'");
    crystal = applied->model.change_basis(cb_op_matrix(applied->cb_op));
    group = user_group_;
    if (!applied->cb_op.is_identity_op())
      for (auto &h : miller)
        h = apply_cb_op(applied->cb_op, h);
  }

  ChainResult chain;
  try {
    chain =
        refine_chain(frame, crystal, group, rows, miller, &trace.post_index);
  } catch (const std::exception &e) {
    throw MethodFailure(e.what());
  }
  trace.seconds_refine = seconds_since(t_refine);

  // _unit_cell_volume_sanity_check.
  const double v0 = crystal.cell.volume, v1 = chain.crystal.cell.volume;
  const double change = std::abs(v0 - v1) / v0;
  if (change > settings_.volume_change_limit)
    throw MethodFailure(format_double(
        "Unrealistic unit cell volume increase during refinement of %.1f%%.",
        100 * change));

  out.crystal = chain.crystal;
  out.raw_B = chain.raw_B;
  out.space_group_ops = operators_of(group);
  out.half_mosaicity_deg = chain.nave.half_mosaicity_deg;
  out.domain_size_ang = chain.nave.domain_size_ang;
  out.ml_half_mosaicity_deg = chain.nave.eta * 180.0 / (2.0 * M_PI);
  out.ml_domain_size_ang = 2.0 / chain.nave.alpha;
  out.row = chain.rows;
  out.miller_index = chain.miller;
  const Predictions &p = chain.predictions;
  const std::size_t n = chain.rows.size();
  out.x_cal_mm = p.x_mm;
  out.y_cal_mm = p.y_mm;
  out.delpsi = p.delpsi;
  out.s1 = p.s1;
  out.x_cal_px.resize(n);
  out.y_cal_px.resize(n);
  out.flags.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = chain.rows[k];
    const auto [px, py] = detector_.panel(frame.strong->spots.panel[i])
                              .mm_to_px(p.x_mm[k], p.y_mm[k]);
    out.x_cal_px[k] = px;
    out.y_cal_px[k] = py;
    out.flags[k] = (frame.strong->flags.empty() ? 0 : frame.strong->flags[i]) |
                   kPredicted | kIndexed | kUsedInRefinement;
  }
}

IndexedStill StillsIndexer::index(const Vec3 &s0, double wavelength,
                                  const StrongSpots &strong) const {
  IndexedStill out;
  Frame frame;
  frame.s0 = s0;
  frame.wavelength = wavelength;
  frame.strong = &strong;
  frame.pixel_size_0 = {detector_.panels[0].pixel_size[0],
                        detector_.panels[0].pixel_size[1]};
  if (!strong.flags.empty() && strong.flags.size() != strong.spots.size())
    throw std::invalid_argument("StillsIndexer: one flag per spot");
  frame.mapped = map_spots(detector_, s0, wavelength, strong.spots);
  out.mapped = frame.mapped;
  frame.d.resize(frame.mapped.rlp.size());
  for (std::size_t i = 0; i < frame.d.size(); ++i) {
    const Vec3 &q = frame.mapped.rlp[i];
    frame.d[i] = 1 / std::sqrt(q.dot(q));
  }

  try {
    index_with_method_list<bool>(
        settings_.methods,
        [&](const std::string &) {
          if (strong.spots.size() == 0)
            throw MethodFailure("No reflections left to index!");
        },
        [&](const std::string &method) {
          out.methods.push_back({});
          MethodTrace &trace = out.methods.back();
          trace.method = method;
          try {
            index_with(method, frame, out, trace);
          } catch (const std::exception &e) {
            trace.error = e.what();
            throw;
          }
          out.method = method;
          return true;
        });
    out.indexed = true;
  } catch (const std::exception &e) {
    out.indexed = false;
    out.reason = e.what();
  }
  return out;
}

} // namespace mxi::stills

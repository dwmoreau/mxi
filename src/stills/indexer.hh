// One still's indexing as dials' stills indexer does it (StillsIndexer.index
// with one lattice, through stills_process's method_list):
//
//   map       the spots' pixel centroids to millimetres, s1 and the
//             reciprocal lattice point
//   search    basis vectors by fft1d or the real-space grid search, on the
//             spots beyond d_min_start; without d_min_start, d_min becomes the
//             spacing of the highest-resolution spot searched
//   candidates  crystal models from triplets of those vectors, Niggli
//             reduced, filtered against the known symmetry, taken in order as
//             they are made and at most max_refine of them
//   per candidate  index the spots; on the known-symmetry path, put the
//             crystal into the target's lattice symmetry first; Sauter-Poon
//             outlier rejection after a refinement, refinement, a second
//             Sauter-Poon pass, a Nave fit, the spots both keep, refinement, a
//             second Nave fit; on the P1 path, a check that the refined
//             crystal still has the target's lattice symmetry; the RMSD in
//             pixels. A candidate whose chain fails is skipped
//   select    the lowest RMSD, the earlier candidate on ties; none, or an
//             RMSD above the threshold, fails the method
//   re-index  every spot with the chosen crystal; on the P1 path, the crystal
//             is put into the target's symmetry and given the user's space
//             group, and the indices follow the change of basis
//   refine    the same chain once more on the indexed spots
//   check     a cell volume that moved by more than half fails the method
//
// The first method that indexes the still gives the result; a method's
// failure is kept and the next tried, and the last failure is the still's.
//
// Follows dials/algorithms/indexing/stills_indexer.py (StillsIndexer.index,
// choose_best_orientation_matrix, identify_outliers, refine), indexer.py and
// lattice_search/__init__.py (BasisVectorSearch), and
// command_line/stills_process.py (Processor._index_impl).
//
// Derived from dials.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "stills/apply_symmetry.hh"
#include "stills/cell.hh"
#include "stills/fft1d.hh"
#include "stills/geometry.hh"
#include "stills/known_symmetry.hh"
#include "stills/nave.hh"
#include "stills/refine.hh"
#include "stills/rsgs.hh"
#include "stills/sgtbx.hh"

namespace mxi::stills {

struct IndexerSettings {
  // indexing.known_symmetry: the cell (a, b, c, alpha, beta, gamma) and the
  // space group as a Hall symbol. The cell is required: max_cell and the grid
  // search come from it.
  std::array<double, 6> unit_cell{};
  std::string space_group_hall;
  double relative_length_tolerance = 0.1;
  double absolute_angle_tolerance = 5;
  double max_delta = 5;

  bool refine_candidates_with_known_symmetry = false;
  // indexing.stills.method_list, or {indexing.method}.
  std::vector<std::string> methods{"fft1d"};
  // refinement_protocol.d_min_start; zero or less is none.
  double d_min_start = 0;
  std::size_t max_refine = 5;
  double xy_rmsd_threshold = 2.0; // pixels
  double hkl_tolerance = 0.3;
  double max_cell_multiplier = 1.3;
  double volume_change_limit = 0.5;

  std::optional<double> fft1d_characteristic_grid;
  Fft1dGrid fft1d_grid;
  double rsgs_characteristic_grid = 0.02;
  std::size_t rsgs_max_vectors = 30;

  RefineOptions refine;
  NaveOptions nave;
};

// A still's strong spots with their reflection flags.
struct StrongSpots {
  Spots spots;
  std::vector<std::uint64_t> flags;
};

// One step of a refinement chain: "refine" (a = the X RMSD, b = the Y RMSD,
// in mm; n_kept = iterations), "sauter_poon" (a = sigma) or "nave" (a = alpha,
// b = eta).
struct StageTrace {
  std::string step;
  std::size_t n_in = 0, n_kept = 0;
  double a = 0, b = 0;
};

// How one refined candidate ended.
struct CandidateTrace {
  std::string status; // "ok", or the stage that failed and why
  std::size_t n_indexed = 0;
  double rmsd = 0;
  std::array<double, 6> cell{};
  double half_mosaicity_deg = 0, domain_size_ang = 0;
  std::vector<StageTrace> stages;
};

struct MethodTrace {
  std::string method;
  std::string error; // empty if the method indexed the still
  std::size_t n_basis_vectors = 0;
  double d_min = 0;
  std::vector<CandidateTrace> candidates;
  int best = -1;
  Mat3 chosen_A; // the chosen candidate's crystal, before re-indexing
  std::vector<StageTrace> post_index;
  double seconds_search = 0, seconds_candidates = 0, seconds_refine = 0;
};

struct IndexedStill {
  bool indexed = false;
  // The failure dials would raise, as its message; empty when indexed.
  std::string reason;
  std::string method;
  std::vector<MethodTrace> methods;

  // The spots mapped to the laboratory frame, every input row.
  MappedSpots mapped;

  // The crystal: U, B, its cell, its space group's operators as x,y,z
  // strings, and the Nave model (widened and as fitted).
  CrystalModel crystal;
  Mat3 raw_B; // B as refinement left it, before dxtbx's set_B
  std::vector<std::string> space_group_ops;
  double half_mosaicity_deg = 0, domain_size_ang = 0;
  double ml_half_mosaicity_deg = 0, ml_domain_size_ang = 0;

  // The refined reflections, one per row kept: its row in the input, its
  // Miller index, and the final prediction.
  std::vector<std::size_t> row;
  std::vector<std::array<int, 3>> miller_index;
  std::vector<double> x_cal_mm, y_cal_mm, x_cal_px, y_cal_px, delpsi;
  std::vector<Vec3> s1;
  std::vector<std::uint64_t> flags;
};

class StillsIndexer {
public:
  StillsIndexer(StillsDetector detector, IndexerSettings settings);

  IndexedStill index(const Vec3 &s0, double wavelength,
                     const StrongSpots &strong) const;

  const StillsDetector &detector() const { return detector_; }
  double max_cell() const { return max_cell_; }
  const IndexerSettings &settings() const { return settings_; }

private:
  struct Frame;
  struct Chosen;
  struct ChainResult;

  void index_with(const std::string &method, Frame &frame, IndexedStill &out,
                  MethodTrace &trace) const;
  std::vector<Vec3> search(const std::string &method,
                           const std::vector<Vec3> &rlp) const;
  ChainResult refine_chain(const Frame &frame, const CrystalModel &start,
                           const sgtbx::SpaceGroup &group,
                           const std::vector<std::size_t> &rows,
                           const std::vector<std::array<int, 3>> &miller,
                           std::vector<StageTrace> *stages) const;

  StillsDetector detector_;
  IndexerSettings settings_;
  sgtbx::SpaceGroup user_group_;
  SymmetryHandler symmetry_;
  KnownSymmetry known_symmetry_;
  double max_cell_ = 0;
  std::unique_ptr<RealSpaceGridSearch> grid_search_;
};

} // namespace mxi::stills

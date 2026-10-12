// mxi_stills_replay: run StillsIndexer over stills written by r10_dump.py,
// without Python, for profiling and for checking that a change leaves every
// result bit for bit the same.
//
//   mxi_stills_replay FRAMES [--out RESULTS] [--times TIMES] [--repeat N]
//                     [--only TAG]
//
// FRAMES holds the indexer settings and, per still, its detector, beam and
// strong spots, with doubles in hexadecimal. RESULTS gets each still's full
// result with doubles in hexadecimal, so two builds can be compared with cmp.
// TIMES gets one line per still and repetition: the seconds index() took and
// the seconds of each method's stages. The indexer is built once and rebuilt
// only when a still's detector differs from the last one's, as the Python
// Workspace does.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "stills/indexer.hh"

using namespace mxi::stills;

namespace {

struct Frame {
  std::string tag;
  StillsDetector detector;
  mxi::Vec3 s0;
  double wavelength = 0;
  StrongSpots strong;
};

struct Input {
  IndexerSettings settings;
  std::vector<Frame> frames;
};

double number(const std::string &token) {
  char *end = nullptr;
  const double v = std::strtod(token.c_str(), &end);
  if (end == token.c_str() || *end != '\0')
    throw std::runtime_error("not a number: " + token);
  return v;
}

std::istringstream expect(std::istream &in, const std::string &key) {
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("unexpected end of input, wanted " + key);
  if (line.compare(0, key.size() + 1, key + " ") != 0)
    throw std::runtime_error("wanted " + key + ", got: " + line.substr(0, 60));
  return std::istringstream(line.substr(key.size() + 1));
}

std::string rest(std::istream &in, const std::string &key) {
  std::istringstream s = expect(in, key);
  return s.str();
}

double read_double(std::istream &s) {
  std::string t;
  if (!(s >> t))
    throw std::runtime_error("missing number");
  return number(t);
}

std::optional<double> read_optional(std::istream &in, const std::string &key) {
  std::istringstream s = expect(in, key);
  std::string t;
  s >> t;
  if (t == "none")
    return std::nullopt;
  return number(t);
}

Input read_input(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  Input out;
  IndexerSettings &s = out.settings;
  expect(in, "mxi_stills_replay");
  expect(in, "dataset");
  expect(in, "config");
  {
    std::istringstream l = expect(in, "unit_cell");
    for (double &v : s.unit_cell)
      v = read_double(l);
  }
  s.space_group_hall = rest(in, "space_group_hall");
  {
    std::istringstream l = expect(in, "methods");
    s.methods.clear();
    for (std::string m; l >> m;)
      s.methods.push_back(m);
  }
  s.refine_candidates_with_known_symmetry =
      std::stoi(rest(in, "refine_candidates_with_known_symmetry")) != 0;
  s.d_min_start = read_optional(in, "d_min_start").value_or(0.0);
  s.max_refine = std::stoul(rest(in, "max_refine"));
  s.xy_rmsd_threshold = number(rest(in, "xy_rmsd_threshold"));
  s.fft1d_characteristic_grid = read_optional(in, "fft1d_characteristic_grid");
  s.rsgs_characteristic_grid = number(rest(in, "rsgs_characteristic_grid"));
  const std::string axis = rest(in, "axis_derivative");
  if (axis == "dials")
    s.refine.axis_derivative = AxisDerivative::dials;
  else if (axis == "analytic")
    s.refine.axis_derivative = AxisDerivative::analytic;
  else
    throw std::runtime_error("axis_derivative: " + axis);
  s.nave.safe_likelihood = std::stoi(rest(in, "safe_likelihood")) != 0;
  s.nave.ewald_proximal_volume_max =
      number(rest(in, "ewald_proximal_volume_max"));
  s.nave.ewald_proximity_resolution_cutoff =
      number(rest(in, "ewald_proximity_resolution_cutoff"));

  const std::size_t n_frames = std::stoul(rest(in, "frames"));
  out.frames.resize(n_frames);
  std::string line;
  for (Frame &f : out.frames) {
    f.tag = rest(in, "frame");
    const std::size_t n_panels = std::stoul(rest(in, "panels"));
    f.detector.panels.resize(n_panels);
    for (StillsPanel &p : f.detector.panels) {
      std::getline(in, line);
      std::istringstream l(line);
      for (mxi::Vec3 *v : {&p.origin, &p.fast, &p.slow}) {
        const double x = read_double(l), y = read_double(l), z = read_double(l);
        *v = mxi::Vec3{x, y, z};
      }
      p.pixel_size[0] = read_double(l);
      p.pixel_size[1] = read_double(l);
      p.mu = read_double(l);
      p.t0 = read_double(l);
    }
    {
      std::istringstream l = expect(in, "s0");
      const double x = read_double(l), y = read_double(l), z = read_double(l);
      f.s0 = mxi::Vec3{x, y, z};
    }
    f.wavelength = number(rest(in, "wavelength"));
    const std::size_t n = std::stoul(rest(in, "spots"));
    Spots &sp = f.strong.spots;
    sp.panel.resize(n);
    sp.x_px.resize(n);
    sp.y_px.resize(n);
    sp.var_x_px.resize(n);
    sp.var_y_px.resize(n);
    f.strong.flags.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      std::getline(in, line);
      std::istringstream l(line);
      std::string t;
      l >> t;
      sp.panel[i] = std::stoll(t);
      sp.x_px[i] = read_double(l);
      sp.y_px[i] = read_double(l);
      sp.var_x_px[i] = read_double(l);
      sp.var_y_px[i] = read_double(l);
      l >> t;
      f.strong.flags[i] = std::stoull(t);
    }
  }
  return out;
}

bool same_detector(const StillsDetector &a, const StillsDetector &b) {
  if (a.size() != b.size())
    return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const StillsPanel &p = a.panels[i], &q = b.panels[i];
    for (int k = 0; k < 3; ++k)
      if (p.origin[k] != q.origin[k] || p.fast[k] != q.fast[k] ||
          p.slow[k] != q.slow[k])
        return false;
    if (p.pixel_size[0] != q.pixel_size[0] ||
        p.pixel_size[1] != q.pixel_size[1] || p.mu != q.mu || p.t0 != q.t0)
      return false;
  }
  return true;
}

void put(std::ostream &o, double v) {
  char b[40];
  std::snprintf(b, sizeof b, " %a", v);
  o << b;
}

void put_stages(std::ostream &o, const std::vector<StageTrace> &stages) {
  for (const StageTrace &t : stages) {
    o << "  stage " << t.step << ' ' << t.n_in << ' ' << t.n_kept;
    put(o, t.a);
    put(o, t.b);
    o << '\n';
  }
}

void write_result(std::ostream &o, const std::string &tag,
                  const IndexedStill &r) {
  o << "frame " << tag << '\n'
    << "indexed " << r.indexed << '\n'
    << "reason " << r.reason << '\n'
    << "method " << r.method << '\n';
  for (const MethodTrace &m : r.methods) {
    o << " method " << m.method << " error " << m.error << " basis "
      << m.n_basis_vectors << " best " << m.best << " d_min";
    put(o, m.d_min);
    o << '\n';
    for (const CandidateTrace &c : m.candidates) {
      o << "  candidate " << c.status << " n " << c.n_indexed << " rmsd";
      put(o, c.rmsd);
      for (double v : c.cell)
        put(o, v);
      put(o, c.half_mosaicity_deg);
      put(o, c.domain_size_ang);
      o << '\n';
      put_stages(o, c.stages);
    }
    o << "  post_index\n";
    put_stages(o, m.post_index);
    o << "  chosen_A";
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        put(o, m.chosen_A(i, j));
    o << '\n';
  }
  o << "mapped";
  for (std::size_t i = 0; i < r.mapped.x_mm.size(); ++i) {
    put(o, r.mapped.x_mm[i]);
    put(o, r.mapped.y_mm[i]);
    put(o, r.mapped.var_x_mm[i]);
    put(o, r.mapped.var_y_mm[i]);
    for (int k = 0; k < 3; ++k) {
      put(o, r.mapped.s1[i][k]);
      put(o, r.mapped.rlp[i][k]);
    }
  }
  o << '\n';
  if (!r.indexed)
    return;
  o << "U B B_raw";
  for (const mxi::Mat3 *m : {&r.crystal.U, &r.crystal.B, &r.raw_B})
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        put(o, (*m)(i, j));
  o << "\ncell";
  for (double v : r.crystal.cell.params)
    put(o, v);
  o << "\nops";
  for (const std::string &op : r.space_group_ops)
    o << ' ' << op;
  o << "\nnave";
  put(o, r.half_mosaicity_deg);
  put(o, r.domain_size_ang);
  put(o, r.ml_half_mosaicity_deg);
  put(o, r.ml_domain_size_ang);
  o << '\n';
  for (std::size_t k = 0; k < r.row.size(); ++k) {
    o << "refl " << r.row[k] << ' ' << r.miller_index[k][0] << ' '
      << r.miller_index[k][1] << ' ' << r.miller_index[k][2] << ' '
      << r.flags[k];
    put(o, r.x_cal_mm[k]);
    put(o, r.y_cal_mm[k]);
    put(o, r.x_cal_px[k]);
    put(o, r.y_cal_px[k]);
    put(o, r.delpsi[k]);
    for (int c = 0; c < 3; ++c)
      put(o, r.s1[k][c]);
    o << '\n';
  }
}

} // namespace

int main(int argc, char **argv) {
  std::string frames_path, out_path, times_path, only;
  int repeat = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto value = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::runtime_error(a + " needs a value");
      return argv[++i];
    };
    if (a == "--out")
      out_path = value();
    else if (a == "--times")
      times_path = value();
    else if (a == "--repeat")
      repeat = std::stoi(value());
    else if (a == "--only")
      only = value();
    else if (frames_path.empty())
      frames_path = a;
    else {
      std::cerr << "mxi_stills_replay: unexpected argument " << a << '\n';
      return 2;
    }
  }
  if (frames_path.empty()) {
    std::cerr << "usage: mxi_stills_replay FRAMES [--out RESULTS] "
                 "[--times TIMES] [--repeat N] [--only TAG]\n";
    return 2;
  }

  const Input input = read_input(frames_path);
  std::ofstream out, times;
  if (!out_path.empty())
    out.open(out_path);
  if (!times_path.empty()) {
    times.open(times_path);
    times << "tag\trepeat\tindexed\tseconds\tmethod\tsearch\tcandidates\t"
             "refine\n";
  }

  std::optional<StillsIndexer> indexer;
  const StillsDetector *built_for = nullptr;
  std::size_t n_builds = 0, n_indexed = 0, n_run = 0;
  double total = 0;
  for (int rep = 0; rep < repeat; ++rep) {
    for (const Frame &f : input.frames) {
      if (!only.empty() && f.tag != only)
        continue;
      if (!built_for || !same_detector(*built_for, f.detector)) {
        indexer.emplace(f.detector, input.settings);
        built_for = &f.detector;
        ++n_builds;
      }
      const auto t0 = std::chrono::steady_clock::now();
      const IndexedStill r = indexer->index(f.s0, f.wavelength, f.strong);
      const double s =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
              .count();
      total += s;
      ++n_run;
      n_indexed += r.indexed;
      if (out.is_open() && rep == 0)
        write_result(out, f.tag, r);
      if (times.is_open())
        for (const MethodTrace &m : r.methods)
          times << f.tag << '\t' << rep << '\t' << r.indexed << '\t' << s
                << '\t' << m.method << '\t' << m.seconds_search << '\t'
                << m.seconds_candidates << '\t' << m.seconds_refine << '\n';
    }
  }
  std::printf("%zu stills, %zu indexed, %zu indexer builds, %.4f s/still\n",
              n_run, n_indexed, n_builds,
              n_run ? total / static_cast<double>(n_run) : 0.0);
  return 0;
}

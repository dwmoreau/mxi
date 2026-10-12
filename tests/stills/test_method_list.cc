// The method_list driver with stub searches: order, the first success, the
// last error rethrown, and set-up errors ending the attempt.

#include <stdexcept>
#include <string>
#include <vector>

#include "check.hh"
#include "stills/method_list.hh"

using namespace mxi::stills;

namespace {

struct SearchFailed : std::runtime_error {
  using std::runtime_error::runtime_error;
};

} // namespace

TEST(stills_method_list_first_success_wins) {
  std::vector<std::string> set_up, ran;
  std::vector<MethodRun> runs;
  const int result = index_with_method_list<int>(
      {"fft1d", "real_space_grid_search"},
      [&](const std::string &m) { set_up.push_back(m); },
      [&](const std::string &m) {
        ran.push_back(m);
        return m == "fft1d" ? 1 : 2;
      },
      &runs);
  check::equal(result, 1, "fft1d's result");
  check::equal(static_cast<long long>(ran.size()), 1, "the second never runs");
  check::equal(static_cast<long long>(set_up.size()), 1, "nor is it set up");
  check::is_true(runs.size() == 1 && !runs[0].error, "one clean run");
}

TEST(stills_method_list_falls_back_in_order) {
  std::vector<std::string> ran;
  std::vector<MethodRun> runs;
  const int result = index_with_method_list<int>(
      {"fft1d", "real_space_grid_search"}, [](const std::string &) {},
      [&](const std::string &m) {
        ran.push_back(m);
        if (m == "fft1d")
          throw SearchFailed("no lattice");
        return 2;
      },
      &runs);
  check::equal(result, 2, "the grid search's result");
  check::is_true(
      ran == std::vector<std::string>{"fft1d", "real_space_grid_search"},
      "tried in order");
  check::is_true(runs.size() == 2 && runs[0].error && !runs[1].error,
                 "the first failed, the second did not");
}

TEST(stills_method_list_rethrows_the_last_error) {
  std::string caught;
  try {
    index_with_method_list<int>(
        {"fft1d", "real_space_grid_search"}, [](const std::string &) {},
        [](const std::string &m) -> int { throw SearchFailed(m + " failed"); });
  } catch (const SearchFailed &e) {
    caught = e.what();
  }
  check::is_true(caught == "real_space_grid_search failed",
                 "the last method's error, unchanged");
}

TEST(stills_method_list_set_up_errors_end_the_attempt) {
  std::vector<std::string> ran;
  bool thrown = false;
  try {
    index_with_method_list<int>(
        {"fft1d", "real_space_grid_search"},
        [](const std::string &m) {
          if (m == "fft1d")
            throw std::runtime_error("set-up failed");
        },
        [&](const std::string &m) {
          ran.push_back(m);
          return 0;
        });
  } catch (const std::runtime_error &) {
    thrown = true;
  }
  check::is_true(thrown && ran.empty(),
                 "a set-up error is not caught and nothing runs");
}

TEST(stills_method_list_passes_on_non_standard_errors) {
  bool thrown = false;
  try {
    index_with_method_list<int>(
        {"fft1d", "real_space_grid_search"}, [](const std::string &) {},
        [](const std::string &) -> int { throw 42; });
  } catch (int) {
    thrown = true;
  }
  check::is_true(thrown, "not caught: the second method is not tried");
}

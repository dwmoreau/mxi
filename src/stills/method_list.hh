// dials.stills_process's indexing.stills.method_list: the indexing methods
// tried in order on one still. Each method's indexer is set up (an error there
// ends the attempt at once) and then run; an error from the run is kept and the
// next method tried. The first method to succeed gives the result; if none
// does, the last method's error is rethrown.
//
// Follows dials/command_line/stills_process.py, Processor._index_impl.
//
// Derived from dials.
// dials: Copyright (c) 2012 Diamond Light Source, Lawrence Berkeley National
// Laboratory, United Kingdom Research and Innovation, and contributors. All
// rights reserved. BSD-3-Clause; see third_party/dials/LICENSE.

#pragma once

#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace mxi::stills {

struct MethodRun {
  std::string method;
  // Empty if the method succeeded.
  std::exception_ptr error;
};

// set_up(method) builds or refreshes that method's indexer; index(method)
// runs it. Errors that dials would not catch (anything not derived from
// std::exception, as dials catches only Exception) pass straight through. runs,
// if given, receives each method tried and how it ended.
template <typename Result>
Result
index_with_method_list(const std::vector<std::string> &methods,
                       const std::function<void(const std::string &)> &set_up,
                       const std::function<Result(const std::string &)> &index,
                       std::vector<MethodRun> *runs = nullptr) {
  if (methods.empty())
    throw std::invalid_argument("index_with_method_list: no methods");
  std::exception_ptr last_error;
  for (const std::string &method : methods) {
    set_up(method);
    try {
      Result result = index(method);
      if (runs)
        runs->push_back({method, nullptr});
      return result;
    } catch (const std::exception &) {
      last_error = std::current_exception();
      if (runs)
        runs->push_back({method, last_error});
    }
  }
  std::rethrow_exception(last_error);
}

} // namespace mxi::stills

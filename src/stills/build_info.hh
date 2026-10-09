// What build of the stills library is running.

#pragma once

#include <string>

namespace mxi::stills {

struct BuildInfo {
  std::string version;
  std::string git_revision;
  std::string compiler;
};

BuildInfo build_info();

}  // namespace mxi::stills

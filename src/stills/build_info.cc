#include "build_info.hh"

namespace mxi::stills {

BuildInfo build_info() {
  return {MXI_STILLS_VERSION, MXI_STILLS_GIT_REVISION,
#if defined(__clang__)
          "clang " __clang_version__
#elif defined(__GNUC__)
          "gcc " __VERSION__
#else
          "unknown"
#endif
  };
}

} // namespace mxi::stills

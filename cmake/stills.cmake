# Stills indexing for dials.stills_process, included from the top-level
# CMakeLists.txt after the mxi library is declared.

# The revision is compiled in so that a running stills_process can say which
# build it loaded: the module is pip-installed, and nothing in the dials,
# dxtbx or cctbx checkouts records it.
execute_process(
  COMMAND git describe --always --dirty --abbrev=12
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE MXI_STILLS_GIT_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT MXI_STILLS_GIT_REVISION)
  set(MXI_STILLS_GIT_REVISION unknown)
endif()

add_library(mxi_stills STATIC
  src/stills/build_info.cc
  src/stills/crystal.cc
  src/stills/geometry.cc
  src/stills/index.cc
  src/stills/lbfgs.cc
  src/stills/nave.cc
  src/stills/outliers.cc
  src/stills/predict.cc
  src/stills/refine.cc
)
target_link_libraries(mxi_stills PUBLIC mxi)
target_compile_definitions(mxi_stills PRIVATE
  MXI_STILLS_VERSION="${PROJECT_VERSION}"
  MXI_STILLS_GIT_REVISION="${MXI_STILLS_GIT_REVISION}")
# dials' geometry is compiled without fused multiply-add, and matching it to
# the last bit needs the same.
target_compile_options(mxi_stills PRIVATE -Wall -Wextra -Wpedantic
  -ffp-contract=off)

add_executable(mxi_stills_tests
  tests/main.cc
  tests/stills/test_geometry.cc
  tests/stills/test_outliers.cc
  tests/stills/test_refine.cc
)
target_include_directories(mxi_stills_tests PRIVATE tests)
target_link_libraries(mxi_stills_tests PRIVATE mxi_stills)
target_compile_options(mxi_stills_tests PRIVATE -Wall -Wextra)
add_test(NAME mxi_stills_tests COMMAND mxi_stills_tests)

if(MXI_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_mxi_stills bindings/module.cc)
  target_link_libraries(_mxi_stills PRIVATE mxi_stills)
  target_compile_options(_mxi_stills PRIVATE -Wall -Wextra)
  install(TARGETS _mxi_stills LIBRARY DESTINATION mxi_stills COMPONENT python)
endif()

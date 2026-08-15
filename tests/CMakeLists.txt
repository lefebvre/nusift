# One test binary, with test_<module>.cpp files named 1:1 with the modules they cover.
#
# CTest labels partition the suite by what a test NEEDS, which is what makes the inner loop
# fast and a clean checkout green:
#
#   unit     -- no data files at all. Must pass on a fresh clone with no store staged.
#   data     -- needs the real staged nuclear-data store.
#   staging  -- needs ENDFtk and ENDF tapes.
#
# `ctest -L unit` is the inner loop; `ctest -LE 'data|staging'` is what CI runs by default.

find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
  message(STATUS "Fetching GoogleTest ${NUSIFT_GTEST_TAG}")
  include(FetchContent)
  # Match the CRT linkage of the rest of the build; without this MSVC fails to link
  # GoogleTest against a dynamically-linked runtime.
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_Declare(googletest
    GIT_REPOSITORY ${NUSIFT_GTEST_REPO}
    GIT_TAG        ${NUSIFT_GTEST_TAG}
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

# Produces the small hand-entered store the end-to-end test and manual CLI checks run
# against. A separate binary rather than test-only code so the same fixture can be generated
# outside the suite, which is what makes `nusift rank` runnable before any ENDF tape exists.
add_executable(nusift_make_test_store tools/make_test_store.cpp)
target_link_libraries(nusift_make_test_store PRIVATE nusift::nusift nusift_warnings)

add_executable(nusift_unit_tests
  unit/test_data_store.cpp
  unit/test_decay_engine.cpp
  unit/test_forecast.cpp
  unit/test_inventory.cpp
  unit/test_inventory_io.cpp
  unit/test_point_source.cpp
  unit/test_report.cpp
  unit/test_seed_fission.cpp
  unit/test_nuclide.cpp
  unit/test_nuclide_name.cpp
  unit/test_time_spec.cpp
  unit/test_triage.cpp
  unit/test_units.cpp)

# Test helpers are included by path relative to unit/, so the directory has to be on the
# include path in its own right. HDF5 is here only so test_data_store can forge a store with
# a future version attribute -- there is no production path that writes one, and the reader's
# forward-compatibility guard is otherwise untestable.
target_include_directories(nusift_unit_tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/unit ${HDF5_C_INCLUDE_DIRS})
target_link_libraries(nusift_unit_tests PRIVATE ${HDF5_C_LIBRARIES})

target_link_libraries(nusift_unit_tests
  PRIVATE nusift::nusift GTest::gtest GTest::gtest_main nusift_warnings)

nusift_clang_tidy(nusift_unit_tests)

# PRE_TEST discovery rather than the default POST_BUILD: it defers enumerating tests until
# ctest runs, which keeps cross-compiled and multi-config (MSVC) builds working. cram uses
# the same mode for the same reason.
gtest_discover_tests(nusift_unit_tests
  TEST_PREFIX nusift.
  PROPERTIES LABELS unit
  DISCOVERY_MODE PRE_TEST)

# ---------------------------------------------------------------------------
# Packaging smoke test
#
# Builds an external consumer against the INSTALLED nusift package. This catches export
# regressions that an in-tree build structurally cannot see: a source-tree path leaking into
# the exported interface, a PUBLIC dependency missing from the export set, or a public header
# that was added but never installed. Modeled on cram's own install-consume CI job.
#
# Only meaningful when there is something to install, so it is gated on NUSIFT_INSTALL --
# which is itself off for fetched-cram and staging builds.
# ---------------------------------------------------------------------------
if(NUSIFT_INSTALL)
  add_test(NAME nusift.package.find_package
    COMMAND ${CMAKE_COMMAND}
      -DNUSIFT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
      -DNUSIFT_BINARY_DIR=${CMAKE_BINARY_DIR}
      -DNUSIFT_CONFIG=$<CONFIG>
      -DNUSIFT_GENERATOR=${CMAKE_GENERATOR}
      -DNUSIFT_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DNUSIFT_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/package/run_package_test.cmake)
  set_tests_properties(nusift.package.find_package PROPERTIES LABELS package)
endif()

cmake_minimum_required(VERSION 3.28)
# The single source of truth for NuSIFT's version. Four things derive from it and nothing
# else declares one: the generated nusift/version.hpp, the nusiftConfigVersion.cmake that
# answers find_package(nusift <version>), `nusift --version`, and the nusift_version root
# attribute stamped into every staged data store. The release tag is checked against it by
# the version-tag-check CI job, so a tag without the matching bump here fails rather than
# shipping a binary that misreports itself. Same discipline as cram.
project(nusift VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Defines CMAKE_INSTALL_INCLUDEDIR / LIBDIR / BINDIR / DATADIR. Included this early because
# target_include_directories() below needs INCLUDEDIR for its INSTALL_INTERFACE, not just
# the install rules at the end of this file.
include(GNUInstallDirs)

# Cache FetchContent downloads outside the build tree so they survive clean builds and can
# be shared across build directories. Same rationale, and same directory name, as cram.
set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/external CACHE PATH
    "Directory for FetchContent downloads")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/dep_versions.cmake)

# ---------------------------------------------------------------------------
# Options
#
# All NUSIFT_-prefixed, because cram is pulled into this build as a sub-build and carries
# identically-named options (ENABLE_TESTS, ENABLE_COVERAGE, ENABLE_SANITIZERS, WITH_ENDFTK).
# An unprefixed option here would silently drive cram's.
# ---------------------------------------------------------------------------
option(NUSIFT_ENABLE_TESTS      "Build the NuSIFT test suite"                      ON)
option(NUSIFT_BUILD_CLI         "Build the nusift command-line tool"               ON)
option(NUSIFT_WITH_STAGING      "Build the ENDF->store staging tool (needs ENDFtk)" OFF)
option(NUSIFT_BUILD_PYTHON      "Build the nanobind Python extension"              OFF)
option(NUSIFT_ENABLE_COVERAGE   "Instrument for code coverage (gcov/gcovr)"        OFF)
option(NUSIFT_ENABLE_SANITIZERS "Build with Address + UB sanitizers"               OFF)
option(NUSIFT_CLANG_TIDY        "Run clang-tidy during the build"                  OFF)
option(NUSIFT_BUILD_DOCS        "Generate the Doxygen API reference"               OFF)
# NUSIFT_INSTALL is declared further down: its default depends on whether cram and Eigen end
# up found or fetched, which is not known until after the dependency block below.

# ---------------------------------------------------------------------------
# Whole-build compile settings
#
# The two settings here are deliberately GLOBAL, unlike the per-target instrumentation further
# down, and they are set here -- above the dependency block -- because a fetched dependency is
# configured at the point it is made available. Set any later and cram, and through it Eigen,
# would already have been configured without them.
# ---------------------------------------------------------------------------
if(NUSIFT_ENABLE_SANITIZERS)
  if(MSVC)
    message(FATAL_ERROR "NUSIFT_ENABLE_SANITIZERS requires a GCC/Clang toolchain.")
  endif()
  # A sanitizer build is only sound when every object linked together is instrumented, and with
  # Eigen in the mix an uninstrumented dependency does not merely miss diagnostics -- it invents
  # one. Eigen's Memory.h keys EIGEN_MALLOC_ALREADY_ALIGNED off __SANITIZE_ADDRESS__, which GCC
  # defines per translation unit, so instrumenting NuSIFT alone left cram allocating through
  # plain malloc while NuSIFT freed through the hand-made aligned allocator, which looks eight
  # bytes in front of the block for a pointer that was never written there. Every test that ran
  # a solve died on a heap-buffer-overflow attributed to Eigen rather than to the mismatch.
  add_compile_options(-fsanitize=address,undefined -g)
  add_link_options(-fsanitize=address,undefined)
endif()

# The Python extension is a shared module, so every static library linked into it -- nusift's
# and cram's alike -- has to be position-independent. Without this an ELF link fails outright
# ("relocation R_X86_64_PC32 against symbol ... can not be used when making a shared object"),
# which is a link-time error in the extension rather than anything visible in the libraries
# themselves. Windows has no such distinction, so this is invisible there in both directions.
if(NUSIFT_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# ---------------------------------------------------------------------------
# cram-depletion: find first, fetch as fallback
#
# The two modes must be interchangeable, so everything downstream links cram::cram and
# nothing branches on which one happened. The only asymmetry is that cram namespaces its
# targets on INSTALL, not in the build tree, so the fetch branch adds the alias itself.
# ---------------------------------------------------------------------------
find_package(cram ${NUSIFT_CRAM_MIN_VERSION} CONFIG QUIET)

if(cram_FOUND)
  set(NUSIFT_CRAM_FETCHED FALSE)
  message(STATUS "Using installed cram ${cram_VERSION}")
  if(NUSIFT_WITH_STAGING)
    # An installed cram cannot have been built with ENDFtk: cram makes CRAM_ENABLE_INSTALL
    # and CRAM_WITH_ENDFTK mutually exclusive, because ENDFtk is always fetched and so can
    # never appear in an installed cram's export set. Its loadDecayData/loadFissionYields
    # are therefore no-op stubs returning 0 -- staging against it would silently write an
    # EMPTY store rather than fail, which is the worst possible outcome. Diagnosed here.
    message(FATAL_ERROR
      "NUSIFT_WITH_STAGING is ON but an installed cram was found. An installed cram is "
      "built without ENDFtk (cram forbids CRAM_ENABLE_INSTALL with CRAM_WITH_ENDFTK), so "
      "its ENDF readers are no-op stubs and staging would write an empty store. Configure "
      "a separate staging build tree with -DCMAKE_DISABLE_FIND_PACKAGE_cram=ON so cram is "
      "fetched and built with CRAM_WITH_ENDFTK=ON.")
  endif()
else()
  set(NUSIFT_CRAM_FETCHED TRUE)
  message(STATUS "Fetching cram ${NUSIFT_CRAM_VERSION} from ${NUSIFT_CRAM_REPO}")

  # cram's options are driven here, BEFORE MakeAvailable, because they are read at its
  # configure time. FORCE because a stale cache entry from a previous configure would
  # otherwise win over the value this build needs.
  set(CRAM_ENABLE_TESTS      OFF CACHE BOOL "cram: tests off in sub-build"      FORCE)
  set(CRAM_ENABLE_BENCHMARKS OFF CACHE BOOL "cram: benchmarks off in sub-build" FORCE)
  # ENDFtk is reached only transitively, and only for staging. It is never a NuSIFT
  # dependency and never present in a runtime build.
  if(NUSIFT_WITH_STAGING)
    set(CRAM_WITH_ENDFTK ON  CACHE BOOL "cram: ENDF readers needed for staging" FORCE)
  else()
    set(CRAM_WITH_ENDFTK OFF CACHE BOOL "cram: ENDF readers not needed"         FORCE)
  endif()

  include(FetchContent)
  FetchContent_Declare(cram
    GIT_REPOSITORY ${NUSIFT_CRAM_REPO}
    GIT_TAG        ${NUSIFT_CRAM_VERSION}
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(cram)

  # cram applies its cram:: namespace at install(EXPORT) time only, so the alias does not
  # exist in a sub-build. Defining it here is what lets every consumer below -- and the
  # install branch -- write cram::cram unconditionally. Guarded because a future cram may
  # define it itself, at which point this becomes dead code rather than a redefinition error.
  if(NOT TARGET cram::cram)
    add_library(cram::cram ALIAS cram)
  endif()
endif()

find_package(HDF5 REQUIRED COMPONENTS C)

# ---------------------------------------------------------------------------
# Install policy
#
# Mirrors cram's, for the identical structural reason: a FETCHED dependency is a target in
# THIS build tree and is never installed, so exporting a target whose interface names it is
# not merely non-relocatable but a hard install(EXPORT) error --
#
#   install(EXPORT "nusiftTargets" ...) includes target "nusift" which requires target
#   "cram" that is not in any export set.
#
# Three conditions have to hold:
#   top-level      -- an add_subdirectory'd nusift must not inject its headers and package
#                     config into the parent project's install tree;
#   cram found     -- a fetched cram is a build-tree target (and it in turn requires a found
#                     Eigen, which cram enforces on its own side);
#   not staging    -- a staging build fetches cram with ENDFtk, so neither can be exported.
#
# Defaulting off rather than erroring keeps `cmake -S . -B build` working for the common
# developer setup with no system Eigen. Forcing it ON anyway is diagnosed below with the
# remedy, rather than left to surface as the error above.
# ---------------------------------------------------------------------------
if(PROJECT_IS_TOP_LEVEL AND NOT NUSIFT_CRAM_FETCHED AND NOT NUSIFT_WITH_STAGING)
  set(_nusift_install_default ON)
else()
  set(_nusift_install_default OFF)
endif()
option(NUSIFT_INSTALL "Generate install and export rules" ${_nusift_install_default})

if(NUSIFT_INSTALL)
  if(NUSIFT_CRAM_FETCHED)
    message(FATAL_ERROR
      "NUSIFT_INSTALL is ON but cram was fetched rather than found. A fetched cram is a "
      "target in this build tree and is never installed, so exporting nusift against it "
      "cannot produce a usable package. Install cram (which itself needs a system Eigen) "
      "and reconfigure, or build with -DNUSIFT_INSTALL=OFF.")
  endif()
  if(NUSIFT_WITH_STAGING)
    message(FATAL_ERROR
      "NUSIFT_INSTALL and NUSIFT_WITH_STAGING are mutually exclusive. Staging requires a "
      "cram built with CRAM_WITH_ENDFTK=ON, and ENDFtk is always fetched and never "
      "installed, so such a build can never be exported. Use two build trees: a runtime "
      "one that installs, and a staging one that produces the .h5 store.")
  endif()
endif()

# ---------------------------------------------------------------------------
# Shared instrumentation flags (coverage)
# INTERFACE target rather than global flags so instrumentation is opt-in per target and the
# Python extension and staging tool can stay out of it. cram's pattern.
#
# Coverage stays per-target where sanitizers could not: a coverage report over vendored code
# is noise rather than a wrong answer, and gcov has no equivalent of Eigen's per-TU allocator
# switch. See the whole-build settings above for why the sanitizers moved.
# ---------------------------------------------------------------------------
add_library(nusift_quality_flags INTERFACE)
if(NUSIFT_ENABLE_COVERAGE)
  if(MSVC)
    message(FATAL_ERROR "NUSIFT_ENABLE_COVERAGE requires a GCC/Clang toolchain (gcov).")
  endif()
  target_compile_options(nusift_quality_flags INTERFACE --coverage -O0 -g)
  target_link_options(nusift_quality_flags INTERFACE --coverage)
endif()

# --- warnings ---------------------------------------------------------------
# On for every first-party target and gated in CI. Kept separate from quality_flags so
# warnings apply even to builds that opt out of instrumentation.
add_library(nusift_warnings INTERFACE)
if(MSVC)
  target_compile_options(nusift_warnings INTERFACE /W4 /permissive-)
else()
  target_compile_options(nusift_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()

# --- clang-tidy -------------------------------------------------------------
# The absolute --header-filter override is deliberate: the repo root is named "nusift", so
# the .clang-tidy file's relative 'nusift/.*' regex would also match vendored dependencies
# under external/ and _deps/. An absolute path rooted at the source tree cannot collide with
# those. This form works with clang-tidy 18, unlike ExcludeHeaderFilterRegex.
function(nusift_clang_tidy target)
  if(NUSIFT_CLANG_TIDY)
    set_target_properties(${target} PROPERTIES
      CXX_CLANG_TIDY "${NUSIFT_CLANG_TIDY_COMMAND}")
  endif()
endfunction()
if(NUSIFT_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy REQUIRED)
  set(NUSIFT_CLANG_TIDY_COMMAND
      "${CLANG_TIDY_EXE};--header-filter=${CMAKE_SOURCE_DIR}/nusift/.*\\.hpp$")
endif()

# ---------------------------------------------------------------------------
# Generated version header
#
# Written into the build tree rather than the source tree so the version can never be edited
# in two places and drift. configure_file also makes the generated header depend on both the
# template and the CMake run, so bumping project(VERSION) rebuilds everything including it.
# ---------------------------------------------------------------------------
set(NUSIFT_GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in
  ${NUSIFT_GENERATED_INCLUDE_DIR}/nusift/version.hpp
  @ONLY)

add_subdirectory(nusift)

if(NUSIFT_BUILD_CLI OR NUSIFT_WITH_STAGING)
  add_subdirectory(nusift_apps)
endif()

if(NUSIFT_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(NUSIFT_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Install / export
# ---------------------------------------------------------------------------
if(NUSIFT_INSTALL)
  include(CMakePackageConfigHelpers)

  install(TARGETS nusift nusift_quality_flags nusift_warnings
    EXPORT nusiftTargets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  if(NUSIFT_BUILD_CLI)
    # Installed but deliberately kept out of the export set: it is an application, and
    # exporting it would put a nusift::nusift_cli target into every consumer's namespace
    # for no reason. cram does the same with its `deplete` driver.
    install(TARGETS nusift_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()

  # Headers are installed by nusift/CMakeLists.txt, which is where the public/internal
  # distinction is made and therefore the only place that can be kept honest.

  install(EXPORT nusiftTargets
    FILE nusiftTargets.cmake
    NAMESPACE nusift::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nusift)

  # Where cram was found, for nusiftConfig.cmake to record as a last-resort search hint (see
  # the rationale there). cram_DIR points at <prefix>/lib/cmake/cram, so the prefix is three
  # levels up. Empty when cram came from somewhere without a config dir, in which case the
  # generated config simply omits the hint.
  set(NUSIFT_CRAM_PREFIX_HINT "")
  if(cram_DIR)
    get_filename_component(NUSIFT_CRAM_PREFIX_HINT "${cram_DIR}/../../.." ABSOLUTE)
  endif()

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/nusiftConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nusiftConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nusift)

  # SameMajorVersion: past 1.0 the major is the breaking-change axis, so a consumer asking
  # for 1.0 is satisfied by any 1.x at or above it. Below 1.0 CMake treats each MINOR as
  # breaking, which is what we want while the API is still moving.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/nusiftConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nusiftConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nusiftConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nusift)

  # The staged nuclear-data store, when one has been built or committed. The CLI's store
  # locator searches <prefix>/share/nusift, so installing it here is what makes an
  # installed nusift work with no --store argument and no environment variable.
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/data/nusift_b8.1.h5)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/data/nusift_b8.1.h5
            DESTINATION ${CMAKE_INSTALL_DATADIR}/nusift)
  endif()
endif()

message(STATUS "")
message(STATUS "NuSIFT ${PROJECT_VERSION}")
if(NUSIFT_CRAM_FETCHED)
  message(STATUS "  cram          : fetched ${NUSIFT_CRAM_VERSION} (build tree; not installable)")
else()
  message(STATUS "  cram          : found ${cram_VERSION} at ${cram_DIR}")
endif()
message(STATUS "  configuration : ${CMAKE_BUILD_TYPE}")
message(STATUS "  tests         : ${NUSIFT_ENABLE_TESTS}")
message(STATUS "  cli           : ${NUSIFT_BUILD_CLI}")
message(STATUS "  staging       : ${NUSIFT_WITH_STAGING}")
message(STATUS "  python        : ${NUSIFT_BUILD_PYTHON}")
message(STATUS "  install/export: ${NUSIFT_INSTALL}")
message(STATUS "")

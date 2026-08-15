# Application binaries. Two, deliberately separate:
#
#   nusift             -- the runtime driver, always built
#   nusift_stage_data  -- the offline ENDF->HDF5 staging tool, gated behind
#                         NUSIFT_WITH_STAGING because it is the only thing in the tree that
#                         needs ENDFtk, and a build carrying ENDFtk can never be installed
#                         (see the install policy in the top-level CMakeLists).

# CLI11 is needed by the applications, never by the library, so it is fetched here rather
# than at the top level. Both apps parse arguments, so this sits outside either gate --
# a staging-only build (NUSIFT_BUILD_CLI=OFF, NUSIFT_WITH_STAGING=ON) is a legitimate and
# in fact expected configuration, and it needs CLI11 just the same.
# find_package first so a system or vcpkg CLI11 wins over a fetched copy.
find_package(CLI11 CONFIG QUIET)
if(NOT CLI11_FOUND)
  message(STATUS "Fetching CLI11 ${NUSIFT_CLI11_TAG}")
  include(FetchContent)
  FetchContent_Declare(cli11
    GIT_REPOSITORY ${NUSIFT_CLI11_REPO}
    GIT_TAG        ${NUSIFT_CLI11_TAG}
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(cli11)
endif()

if(NUSIFT_BUILD_CLI)
  # OUTPUT_NAME rather than naming the target `nusift`: that name is taken by the library.
  add_executable(nusift_cli nusift.cpp)
  set_target_properties(nusift_cli PROPERTIES OUTPUT_NAME nusift)
  target_link_libraries(nusift_cli PRIVATE nusift::nusift CLI11::CLI11 nusift_warnings)
  nusift_clang_tidy(nusift_cli)
endif()

if(NUSIFT_WITH_STAGING)
  # The only target in the tree that touches ENDFtk. It is reached transitively through cram
  # (configured with CRAM_WITH_ENDFTK=ON above), so ENDFtk never appears in NuSIFT's own
  # dependency list and never reaches a runtime build.
  #
  # cram::cram is linked directly here rather than through nusift::nusift because this tool
  # needs cram's ENDF readers and its chain, which the library deliberately PIMPLs away.
  add_executable(nusift_stage_data nusift_stage_data.cpp)
  target_link_libraries(nusift_stage_data
    PRIVATE nusift::nusift cram::cram CLI11::CLI11 nusift_warnings)
  nusift_clang_tidy(nusift_stage_data)
endif()

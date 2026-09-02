# The NuSIFT library.
#
# Headers and sources are co-located per module (no include/ vs src/ split) and included
# fully qualified from the repo root -- #include "nusift/core/nuclide_name.hpp" -- which is
# why the include root below is the directory CONTAINING nusift/, not nusift/ itself.

add_library(nusift
  core/nuclide_name.cpp
  nucdata/store_arrays.cpp
  nucdata/coefficient_pack.cpp
  nucdata/data_store.cpp
  nucdata/nuclear_data.cpp
  engine/inventory.cpp
  engine/decay_engine.cpp
  engine/adjoint_engine.cpp
  engine/reconcile.cpp
  engine/sensitivity.cpp
  nucdata/store_locator.cpp
  nucdata/fission_yield.cpp
  seed/seed_fission.cpp
  exposure/air_coefficients.cpp
  exposure/dose_coefficients.cpp
  exposure/point_source.cpp
  triage/response.cpp
  triage/ranking.cpp
  triage/attribution.cpp
  triage/forecast.cpp
  triage/events.cpp
  triage/allowable.cpp
  triage/intervention.cpp
  triage/spectrum.cpp
  triage/triage_set.cpp
  triage/task_plan.cpp
  triage/uncertainty.cpp
  io/time_spec.cpp
  io/inventory_io.cpp
  io/report.cpp
  io/source_report.cpp)

add_library(nusift::nusift ALIAS nusift)

# A bare ${CMAKE_SOURCE_DIR} here would not merely be wrong for consumers but fatal:
# install(EXPORT) rejects an absolute path into the source tree outright. The generated
# include root carries nusift/version.hpp only; both roots collapse to the same
# ${CMAKE_INSTALL_INCLUDEDIR} once installed, where the generated header sits beside the
# checked-in ones. Same arrangement as cram.
target_include_directories(nusift PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<BUILD_INTERFACE:${NUSIFT_GENERATED_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(nusift PUBLIC cxx_std_20)

# cram and HDF5 are PRIVATE, and that is the whole point of the PIMPL in nusift/nucdata and
# nusift/engine: no installed NuSIFT header includes a cram, Eigen, or HDF5 header, so a
# consumer never compiles one. They remain link requirements (nusiftConfig.cmake re-finds
# them), but not compile requirements.
target_link_libraries(nusift
  PUBLIC  nusift_quality_flags
  PRIVATE nusift_warnings cram::cram ${HDF5_C_LIBRARIES})
target_include_directories(nusift PRIVATE ${HDF5_C_INCLUDE_DIRS})

nusift_clang_tidy(nusift)

# ---------------------------------------------------------------------------
# Header installation
#
# Listed per module rather than globbed, so adding a header is a deliberate decision about
# whether it is public. Anything named *_internal.hpp is by construction not installed: those
# expose cram/Eigen types across the PIMPL seam and are included by a handful of .cpp files
# only. The generated version.hpp is installed alongside, matching its include path.
# ---------------------------------------------------------------------------
if(NUSIFT_INSTALL)
  install(FILES
    units.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift)

  install(FILES
    core/element_symbols.hpp
    core/error.hpp
    core/nuclide.hpp
    core/nuclide_name.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/core)

  # nuclear_data_internal.hpp and decay_engine_internal.hpp are deliberately absent: they
  # expose cram and Eigen types across the PIMPL seam and are included only by the three
  # translation units that assemble or solve a matrix. Installing either would put Eigen back
  # on every consumer's include path.
  install(FILES
    nucdata/coefficient_pack.hpp
    nucdata/data_store.hpp
    nucdata/nuclear_data.hpp
    nucdata/fission_yield.hpp
    nucdata/photon_lines.hpp
    nucdata/store_arrays.hpp
    nucdata/store_locator.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/nucdata)

  install(FILES
    engine/adjoint_engine.hpp
    engine/decay_engine.hpp
    engine/decay_result.hpp
    engine/inventory.hpp
    engine/reconcile.hpp
    engine/sensitivity.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/engine)

  install(FILES
    exposure/air_coefficients.hpp
    exposure/dose_coefficients.hpp
    exposure/point_source.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/exposure)

  install(FILES
    triage/allowable.hpp
    triage/attribution.hpp
    triage/events.hpp
    triage/forecast.hpp
    triage/intervention.hpp
    triage/ranking.hpp
    triage/response.hpp
    triage/spectrum.hpp
    triage/task_plan.hpp
    triage/triage_set.hpp
    triage/uncertainty.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/triage)

  install(FILES
    seed/fission_energy.hpp
    seed/seed_fission.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/seed)

  install(FILES
    io/inventory_io.hpp
    io/number_format.hpp
    io/report.hpp
    io/source_report.hpp
    io/time_spec.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift/io)

  install(FILES
    ${NUSIFT_GENERATED_INCLUDE_DIR}/nusift/version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nusift)
endif()

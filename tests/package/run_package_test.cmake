# Drives the packaging smoke test as a CTest case.
#
# Four steps, each of which must succeed:
#   1. install nusift from the already-built tree into a scratch prefix
#   2. configure tests/package/ against that prefix, and nothing else
#   3. build it
#   4. run it
#
# Step 2 is the one that does the real work. It runs in a fresh binary directory with only
# CMAKE_PREFIX_PATH pointing at the staged prefix, so if nusiftConfig.cmake forgets a
# find_dependency, or nusiftTargets.cmake carries an absolute path into the source tree, or
# a public header was never installed, this is where it fails -- immediately, and with a
# message that names the cause.
#
# Invoked with -P, so the variables below arrive as -D definitions from tests/CMakeLists.txt.

if(NOT DEFINED NUSIFT_SOURCE_DIR OR NOT DEFINED NUSIFT_BINARY_DIR)
  message(FATAL_ERROR "run_package_test.cmake requires NUSIFT_SOURCE_DIR and NUSIFT_BINARY_DIR")
endif()

set(_prefix "${NUSIFT_BINARY_DIR}/package_test/prefix")
set(_consumer_build "${NUSIFT_BINARY_DIR}/package_test/consumer")

# Start from a clean slate so a previous run's staged headers cannot mask one that the
# install rules have since stopped installing.
file(REMOVE_RECURSE "${_prefix}" "${_consumer_build}")

function(run_step description)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _rc OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "package test: ${description} failed (exit ${_rc})\n"
                        "--- stdout ---\n${_out}\n--- stderr ---\n${_err}")
  endif()
endfunction()

# 1. Install into the scratch prefix.
run_step("install"
  ${CMAKE_COMMAND} --install "${NUSIFT_BINARY_DIR}"
                   --prefix "${_prefix}"
                   --config "${NUSIFT_CONFIG}")

# 2. Configure the consumer. CMAKE_PREFIX_PATH points at the staged nusift prefix and NOTHING
#    ELSE -- deliberately. cram's prefix is withheld so this asserts that nusiftConfig.cmake is
#    self-sufficient: it must locate its own dependencies from what it recorded at configure
#    time, because a downstream consumer has no reason to know nusift links cram at all. The
#    toolchain file is forwarded because on Windows it is what locates HDF5 and Eigen.
set(_configure_args
  -S "${NUSIFT_SOURCE_DIR}/tests/package"
  -B "${_consumer_build}"
  -G "${NUSIFT_GENERATOR}"
  "-DCMAKE_PREFIX_PATH=${_prefix}"
  "-DCMAKE_BUILD_TYPE=${NUSIFT_CONFIG}")
if(NUSIFT_CXX_COMPILER)
  list(APPEND _configure_args "-DCMAKE_CXX_COMPILER=${NUSIFT_CXX_COMPILER}")
endif()
if(NUSIFT_TOOLCHAIN_FILE)
  list(APPEND _configure_args "-DCMAKE_TOOLCHAIN_FILE=${NUSIFT_TOOLCHAIN_FILE}")
endif()
run_step("configure consumer" ${CMAKE_COMMAND} ${_configure_args})

# 3. Build it.
run_step("build consumer"
  ${CMAKE_COMMAND} --build "${_consumer_build}" --config "${NUSIFT_CONFIG}")

# 4. Run it. Multi-config generators put the binary in a per-config subdirectory; single-
#    config ones do not, so accept either rather than guessing from the generator name.
find_program(_use_nusift
  NAMES use_nusift
  PATHS "${_consumer_build}" "${_consumer_build}/${NUSIFT_CONFIG}"
  NO_DEFAULT_PATH REQUIRED)
run_step("run consumer" "${_use_nusift}")

message(STATUS "package test: find_package(nusift) consumer built and ran against ${_prefix}")

# Run clang-format over the repository, in check or fix mode.
#
# Invoked in script mode (cmake -P) by the `format-check` and `format` targets rather than being
# expanded at configure time, and that is the whole point: the file list comes from `git ls-files`
# at the moment the target runs, so a source added since the last `cmake` configure is still
# checked. A configure-time glob would go stale exactly when it matters -- on the commit that adds
# a file -- and would send a green local check into a red CI.
#
# This IS the CI gate: .github/workflows/ci.yml invokes this same script in script mode, needing
# no configure step to do it. So a local `format-check` and the CI job cannot disagree about which
# files are checked or how -- which is the failure this target exists to prevent, not merely a
# faster version of it.
#
# One consequence of `git ls-files` worth knowing: a file that has never been `git add`ed is not
# checked, here or in CI, because neither can see it. Stage a new source before trusting a pass.
#
#   cmake -DNUSIFT_SOURCE_DIR=<dir> [-DNUSIFT_FORMAT_FIX=ON] [-DCLANG_FORMAT_EXE=<path>]
#         -P cmake/format.cmake
#
# clang-format is located here rather than at configure time so that installing it does not
# require re-running cmake -- the targets exist either way, and say what to install when it is
# missing.

if(NOT DEFINED NUSIFT_SOURCE_DIR)
  message(FATAL_ERROR "format.cmake: NUSIFT_SOURCE_DIR is required")
endif()

# The clang-format major version the CI gate runs, from ubuntu-24.04's package. Formatting is not
# stable across major versions, so a developer on a different one can be green locally and red in
# CI over rules neither of them wrote. Warned about rather than enforced: a mismatch is usually
# harmless, and refusing to run would be worse than saying so.
set(NUSIFT_CI_CLANG_FORMAT_MAJOR 18)

if(NOT CLANG_FORMAT_EXE OR NOT EXISTS "${CLANG_FORMAT_EXE}")
  find_program(CLANG_FORMAT_EXE NAMES clang-format clang-format-18 clang-format-17)
endif()
if(NOT CLANG_FORMAT_EXE)
  message(FATAL_ERROR
    "clang-format was not found on PATH.\n"
    "  Debian/Ubuntu : sudo apt-get install clang-format\n"
    "  Fedora/RHEL   : sudo dnf install clang-tools-extra\n"
    "  macOS         : brew install clang-format\n"
    "  Python        : pip install clang-format\n"
    "CI runs version ${NUSIFT_CI_CLANG_FORMAT_MAJOR}.")
endif()

execute_process(
  COMMAND "${CLANG_FORMAT_EXE}" --version
  OUTPUT_VARIABLE version_text
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
string(REGEX MATCH "([0-9]+)\\.[0-9]+\\.[0-9]+" _match "${version_text}")
set(found_major "${CMAKE_MATCH_1}")

# `git ls-files` is the source of truth, not a glob: it is what CI uses, and it excludes the
# fetched dependencies under external/ without needing a single exclusion pattern to be kept in
# sync with them.
find_package(Git QUIET)
if(NOT Git_FOUND)
  message(FATAL_ERROR
    "format.cmake: git was not found, and the file list comes from `git ls-files` so that this "
    "target and the CI gate check exactly the same files. Run clang-format by hand in a source "
    "tree without git.")
endif()

execute_process(
  COMMAND "${GIT_EXECUTABLE}" ls-files "*.cpp" "*.hpp"
  WORKING_DIRECTORY "${NUSIFT_SOURCE_DIR}"
  OUTPUT_VARIABLE tracked
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE git_result
  ERROR_VARIABLE git_error)
if(NOT git_result EQUAL 0)
  message(FATAL_ERROR "format.cmake: `git ls-files` failed: ${git_error}")
endif()
string(REPLACE "\n" ";" files "${tracked}")
list(REMOVE_ITEM files "")
if(NOT files)
  message(STATUS "no tracked C++ sources to format")
  return()
endif()
list(LENGTH files count)

if(NUSIFT_FORMAT_FIX)
  execute_process(
    COMMAND "${CLANG_FORMAT_EXE}" -i ${files}
    WORKING_DIRECTORY "${NUSIFT_SOURCE_DIR}"
    RESULT_VARIABLE result)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "clang-format failed while rewriting sources")
  endif()
  # Deliberately silent about WHICH files changed: clang-format -i does not say, and `git diff`
  # answers it better than any summary invented here would.
  message(STATUS "formatted ${count} files with clang-format ${found_major}")
  return()
endif()

execute_process(
  COMMAND "${CLANG_FORMAT_EXE}" --dry-run --Werror ${files}
  WORKING_DIRECTORY "${NUSIFT_SOURCE_DIR}"
  RESULT_VARIABLE result
  ERROR_VARIABLE diagnostics)
if(result EQUAL 0)
  message(STATUS "clang-format ${found_major}: ${count} files, all formatted")
  if(found_major AND NOT found_major STREQUAL "${NUSIFT_CI_CLANG_FORMAT_MAJOR}")
    message(WARNING
      "clang-format ${found_major} here against ${NUSIFT_CI_CLANG_FORMAT_MAJOR} in CI. Formatting "
      "is not stable across major versions, so this pass does not guarantee the CI gate passes.")
  endif()
  return()
endif()

# The diagnostics carry the file, line and a diff, which is the useful part; the FATAL_ERROR that
# follows carries only what to do about it.
message("${diagnostics}")
message(FATAL_ERROR
  "clang-format ${found_major} found unformatted files (the diffs are above).\n"
  "Fix them in place with:  cmake --build <build dir> --target format")

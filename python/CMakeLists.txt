# The nanobind extension.
#
# nanobind rather than pybind11 for one reason that matters here: nb::ndarray hands out
# genuine zero-copy NumPy views with an owner capsule. NuSIFT's outputs are (times x
# contributors) matrices, and copying a sixty-by-a-thousand table into Python would cost more
# than the solve that produced it. The smaller binaries and faster compiles are a bonus.
#
# nanobind is found from the installed Python package rather than fetched, which is how
# scikit-build-core expects it to be located and what keeps the wheel build reproducible.

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR
  RESULT_VARIABLE NB_RESULT)
if(NOT NB_RESULT EQUAL 0)
  message(FATAL_ERROR
    "NUSIFT_BUILD_PYTHON is ON but nanobind is not importable from ${Python_EXECUTABLE}. "
    "Install it with: ${Python_EXECUTABLE} -m pip install nanobind")
endif()
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(nusift_ext STABLE_ABI NB_STATIC src/nusift_ext.cpp)
set_target_properties(nusift_ext PROPERTIES OUTPUT_NAME _core)
target_link_libraries(nusift_ext PRIVATE nusift::nusift)

# Warnings are deliberately NOT applied here: nanobind's generated glue trips several of them
# and none of the findings would be in NuSIFT's code.

# Land the extension next to the Python sources so an in-tree build is importable without
# installing, which is what the test suite and any notebook in the source tree use.
set_target_properties(nusift_ext PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/nusift"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/nusift")
foreach(config IN ITEMS DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
  set_target_properties(nusift_ext PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY_${config} "${CMAKE_CURRENT_SOURCE_DIR}/nusift"
    RUNTIME_OUTPUT_DIRECTORY_${config} "${CMAKE_CURRENT_SOURCE_DIR}/nusift")
endforeach()

# scikit-build-core installs from here; an in-tree build just leaves it beside the sources.
if(DEFINED SKBUILD)
  install(TARGETS nusift_ext LIBRARY DESTINATION nusift)

  # The staged store, into the directory _data.default_store_path() looks in. It has to be
  # installed here rather than by the root install block: a wheel build fetches cram and so
  # runs with NUSIFT_INSTALL=OFF, which takes that block -- and the store with it -- out of
  # the build entirely. Without this, `NuclearData.open()` with no argument works only from a
  # directory that happens to have a ./data beside it.
  #
  # One named file rather than a glob over data/: that directory also holds small fixtures,
  # and default_store_path() takes the first store it finds, so shipping them would make
  # which evaluation a wheel opens depend on alphabetical order.
  if(EXISTS ${PROJECT_SOURCE_DIR}/data/nusift_b8.1.h5)
    install(FILES ${PROJECT_SOURCE_DIR}/data/nusift_b8.1.h5 DESTINATION nusift/data)
  else()
    message(WARNING
      "No data/nusift_b8.1.h5 in the source tree, so this wheel will ship without a nuclear-"
      "data store. nusift.NuclearData.open() will then need $NUSIFT_DATA_STORE or an explicit "
      "path.")
  endif()
endif()

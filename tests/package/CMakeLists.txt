# An external consumer of the INSTALLED nusift package.
#
# Deliberately minimal and deliberately NOT part of the main build: it is configured from
# scratch against a staged install prefix by run_package_test.cmake, so it sees exactly what
# a downstream project sees and nothing of NuSIFT's own build tree. That is the only way to
# catch a source-tree path leaking into the exported interface, a dependency missing from the
# export set, or a public header that was added but never installed.

cmake_minimum_required(VERSION 3.28)
project(nusift_package_test LANGUAGES CXX)

find_package(nusift REQUIRED)

add_executable(use_nusift main.cpp)
target_link_libraries(use_nusift PRIVATE nusift::nusift)

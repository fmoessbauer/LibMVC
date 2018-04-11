cmake_minimum_required(VERSION 3.2)

project(LibMVC)

set(LIBMVC_ENABLE_BENCH   0 CACHE STRING "Build Benchmarks")
set(LIBMVC_ENABLE_TESTING 0 CACHE STRING "Build Tests")
set(LIBMVC_ENABLE_DOCS    1 CACHE STRING "Build Documentation")

if(${LIBMVC_ENABLE_BENCH})
  message("Build Benchmarks")

  # build google benchmark
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
  add_subdirectory(vendor/google/benchmark)
  target_include_directories(
    benchmark PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/google/benchmark/include>)

  # build LibMVC benchmarks
  add_subdirectory(bench)
endif()

#build tests
if(${LIBMVC_ENABLE_TESTING})
  message("Build Testsuite")
  add_subdirectory(vendor/google/googletest/googletest)

  add_subdirectory(test)
endif()

if(${LIBMVC_ENABLE_DOCS})
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    message("Build Documentation - build using make doc")
    set(DOXYGEN_OUTPUT_DIRECTORY docs)
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
    doxygen_add_docs(
      "doc"
      README.md LibNuMVC LibFastVC ParallelSolverAdapter)
  endif()
endif()

message("Build Standalone Solvers")
add_subdirectory(LibNuMVC)
add_subdirectory(LibFastVC)
add_subdirectory(ParallelSolverAdapter)


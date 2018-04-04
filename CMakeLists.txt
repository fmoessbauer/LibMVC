cmake_minimum_required(VERSION 3.2)

project(LibMVC)

set(LIBMVC_ENABLE_BENCH 0 CACHE STRING "Build Benchmarks")
set(LIBMVC_ENABLE_TESTING 0 CACHE STRING "Build Tests")

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

message("Build Standalone Solvers")
add_subdirectory(LibNuMVC)
add_subdirectory(LibFastVC)


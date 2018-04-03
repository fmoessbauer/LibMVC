cmake_minimum_required(VERSION 3.2)

project(LibMVC)

# build google benchmark
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
add_subdirectory(vendor/google/benchmark)
target_include_directories(
  benchmark PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/google/benchmark/include>)

# build LibMVC benchmarks
add_subdirectory(bench)


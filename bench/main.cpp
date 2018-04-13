#include <benchmark/benchmark.h>

#include <limits>
#include <fstream>
#include <chrono>

#include "../LibNuMVC/numvc.hpp"
#include "../LibFastVC/fastvc.hpp"
#include "../ParallelSolverAdapter/ParallelSolverAdapter.hpp"

using namespace libmvc;

using duration_s = std::chrono::duration<std::chrono::seconds>;

template<
  typename Solver>
static void BM_Solver(
    benchmark::State& state,
    const std::string filename,
    const int coversize,
    const std::chrono::seconds time,
    unsigned int seed)
{
  int successes = 0;
  int min_cover = std::numeric_limits<int>::max();
  for (auto _ : state){
    state.PauseTiming();
    std::ifstream file(filename, std::ios::in);
    if(!file){
      auto error = std::string("No such file: ") + filename;
      state.SkipWithError(error.c_str());
    }
    // If no seed given, use a fresh seed per iteration
    if(seed == 0){
      seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    state.ResumeTiming();

    Solver solver(file, coversize, time, false, seed);
    solver.cover_LS();

    state.PauseTiming();
    auto current_cover = solver.get_best_cover_size();
    if(current_cover == coversize) ++successes;
    min_cover = std::min(min_cover, current_cover);
    file.close();
    state.ResumeTiming();
  }
  state.counters["best-cover"]   = min_cover;
  state.counters["success-rate"] = static_cast<double>(successes) / state.iterations();
}

int main(int argc, char** argv){
  benchmark::Initialize(&argc, argv);

  if(argc < 3 || argc > 5){
    std::cout << "usage: <filename> <coversize> <time(s)> [<seed>]" << std::endl;
    return 1;
  }

  unsigned int seed = 0;
  if(argc == 5) {
    seed = std::stoi(argv[4]);
  }

  std::string filepath(argv[1]);

  int coversize = std::stoi(argv[2]);
  auto time_limit(std::chrono::seconds(std::stoi(argv[3])));

  std::string filename(
      filepath.begin()+filepath.find_last_of("/")+1,
      filepath.end());

  std::string config(filename + "/"
      + std::to_string(coversize) + "/" 
      + std::to_string(time_limit.count()));

  // Register benchmarks
  benchmark::RegisterBenchmark(
      (std::string("NuMVC/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<NuMVC>(st, filepath, coversize, time_limit, seed);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("FastVC/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<FastVC>(st, filepath, coversize, time_limit, seed);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("ParallelSolverAdapter<NuMVC>/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<ParallelSolverAdapter<NuMVC>>(st, filepath, coversize, time_limit, seed);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("ParallelSolverAdapter<FastVC>/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<ParallelSolverAdapter<FastVC>>(st, filepath, coversize, time_limit, seed);
        })->Unit(benchmark::kMillisecond);

  benchmark::RunSpecifiedBenchmarks();
}

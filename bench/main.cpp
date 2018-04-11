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
    const std::chrono::seconds time)
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
    state.ResumeTiming();

    Solver solver(file, coversize, time);
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
  if(argc < 3){
    std::cout << "usage: <filename> <coversize> <time(s)>" << std::endl;
    return 1;
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

  benchmark::Initialize(&argc, argv);

  // Register benchmarks
  benchmark::RegisterBenchmark(
      (std::string("NuMVC/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<NuMVC>(st, filepath, coversize, time_limit);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("FastVC/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<FastVC>(st, filepath, coversize, time_limit);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("ParallelSolverAdapter<NuMVC>/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<ParallelSolverAdapter<NuMVC>>(st, filepath, coversize, time_limit);
        })->Unit(benchmark::kMillisecond);

  benchmark::RegisterBenchmark(
      (std::string("ParallelSolverAdapter<FastVC>/") + config).c_str(), [&](benchmark::State & st){
          return BM_Solver<ParallelSolverAdapter<FastVC>>(st, filepath, coversize, time_limit);
        })->Unit(benchmark::kMillisecond);

  benchmark::RunSpecifiedBenchmarks();
}

/**
 * This script is a prototypic implementation of a parallelized MVC solver.
 *
 * It works with all solvers implementing the LibMVC interface by running
 * n instances in parallel. The program stops if one of the solvers
 * has found a cover of the desired size.
 *
 * Due to its embarrassing parallel nature it should almost perfectly scale.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <functional>
#include <atomic>
#include <cstdlib>

#include "../LibNuMVC/numvc.hpp"

// Use this solver
using SOLVER = NuMVC;

// Global state for collective solving
class SolverState {
  public:
    std::mutex mon_mx;
    std::atomic<int> best_cover;
    int best_solver   = 0;
    int optimal_cover = 0;
  SolverState(int num_vertices, int opt_cover)
    : best_cover(num_vertices),
      optimal_cover(opt_cover){ }
};

// acts as monitor for all solver instances
bool monitor(
    const SOLVER & solver,
    bool better_cover_found,
    unsigned int tid, SolverState * state)
{
  auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        solver.get_best_duration());
  auto cover_size = solver.get_best_cover_size();

  // only lock if actual change present
  if(better_cover_found && cover_size < state->best_cover){
    std::lock_guard<std::mutex> lock(state->mon_mx);
    state->best_cover  = cover_size;
    state->best_solver = tid;
    std::cout << std::setw(2) <<tid
              << ": Better MVC found.\tSize: " << solver.get_best_cover_size()
              << "\tTime: " << std::fixed << std::setw(4)
              << std::setprecision(4) << time_ms.count() << "ms" << std::endl;
  }

  // false -> continue if best cover not found yet
  // true  -> stop calculation
  return (state->best_cover == state->optimal_cover);
}

// invoke threads using this function
void start_solver(
    SOLVER * solver,
    unsigned int seed,
    unsigned int tid,
    SolverState * state)
{
  solver->set_random_seed(seed + tid);
  solver->cover_LS(
      std::bind(monitor, std::placeholders::_1, std::placeholders::_2, tid, state));
  std::cout << "-- " << tid << " terminated" << std::endl;
}

// print solution and stats of best solver
void print_solution(const SOLVER & solver){
  // check solution
  if (solver.check_solution()) {
    auto solution(solver.get_independent_set());
    std::cout << "c Best found vertex cover size = " << solver.get_best_cover_size()
         << std::endl;

    std::cout << "c independent set:" << std::endl;
    for (auto v : solution) std::cout << v << ' ';
    std::cout << std::endl;

    auto solver_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        solver.get_best_duration());
    auto total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        solver.get_total_duration());
    double performance = static_cast<double>(solver.get_best_step()) /
                         solver_time_ms.count() / 1000.0;
    std::cout << "c searchSteps = " << solver.get_best_step() << std::endl;
    std::cout << "c solveTime = " << solver_time_ms.count() << "ms" << std::endl;
    std::cout << "c totalTime = " << total_time_ms.count() << "ms" << std::endl;
    std::cout << "c performance = " << performance << "MT/s" << std::endl;
  } else {
    std::cout << "the solution is wrong." << std::endl;
  }
}

int main(int argc, char * argv[]){

  const auto         num_solvers = std::thread::hardware_concurrency();
  const std::string  filename(argv[1]);
  const int          cover_size  = std::atoi(argv[2]);
  const int          cutoff_time = std::atoi(argv[3]);
  const unsigned int base_seed   = (argc == 5 ?
      std::atoi(argv[4]) :
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  // read in file
  std::ifstream file(filename, std::ios::in);

  // create single solver, copy it
  SOLVER master(file, cover_size, std::chrono::seconds(cutoff_time), false);

  std::vector<std::thread> workers;
  std::vector<SOLVER> solvers(num_solvers, master);
  SolverState global_state(master.get_vertex_count(), cover_size);

  std::cout << "Using " << num_solvers << " parallel instances" << std::endl;
  // start threads
  for(unsigned int i=0; i<num_solvers; ++i){
    workers.emplace_back(start_solver, &solvers[i], base_seed, i, &global_state);
  }

  for(auto & w : workers){
    w.join();
  }

  //print solution
  print_solution(solvers[global_state.best_solver]);
}


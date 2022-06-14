#include "gtest/gtest.h"
#include "numvc.hpp"
#include "fastvc.hpp"
#include "parallelsolveradapter.hpp"

#include <algorithm>
#include <fstream>

using namespace libmvc;

template<typename SOLVER>
class InterfaceTests: public ::testing::Test { };

using SolverTypes = ::testing::Types<NuMVC, FastVC, ParallelSolverAdapter<NuMVC>, ParallelSolverAdapter<FastVC>>;
TYPED_TEST_SUITE(InterfaceTests, SolverTypes, );

extern std::string filename;
extern int cover_size;
extern int instance_size;
extern int cutoff_time;

using Edge     = std::pair<int,int>;
using EdgeList = std::pair<int, std::vector<Edge>>;

TYPED_TEST(InterfaceTests, Constructors){
  using SOLVER=TypeParam;

  std::ifstream file(filename, std::ios::in);
  SOLVER solver_one(
      file,
      cover_size,
      std::chrono::seconds(cutoff_time),
      false,
      42);

  EdgeList edgelist = solver_one.get_instance_as_edgelist();
  EXPECT_EQ(edgelist.first, instance_size);

  SOLVER solver_two(
      edgelist.second,
      edgelist.first,
      cover_size,
      std::chrono::milliseconds(cutoff_time * 1000), // test duration_cast
      false,
      42);

  // test solvers
  solver_one.cover_LS();
  solver_two.cover_LS();

  EXPECT_EQ(solver_one.check_solution(), true);
  EXPECT_EQ(solver_two.check_solution(), true);
}

TYPED_TEST(InterfaceTests, GetterSetter){
  using SOLVER=TypeParam;

  std::ifstream file(filename, std::ios::in);
  SOLVER solver(
      file,
      cover_size,
      std::chrono::seconds(cutoff_time));

  solver.set_cutoff_time(std::chrono::milliseconds(cutoff_time*2));
  solver.set_optimal_size(cover_size+1);
  solver.set_random_seed(42);

  // solve instance
  solver.cover_LS();

  std::vector<char> flaglist(std::move(solver.get_cover_as_flaglist()));
  std::vector<int> cover(std::move(solver.get_cover()));

  // check equality
  unsigned int num_true_flags = std::count(flaglist.begin(), flaglist.end(), true);
  EXPECT_EQ(num_true_flags, cover.size());
  EXPECT_EQ(solver.get_best_cover_size(), static_cast<int>(cover.size()));

  // Is are all vertices which are not in cover
  std::vector<int> ind_set(std::move(solver.get_independent_set()));
  EXPECT_EQ(instance_size - num_true_flags, ind_set.size());

  // Check remaining getters
  EdgeList edgelist(std::move(solver.get_instance_as_edgelist()));
  int vertex_count  = solver.get_vertex_count();
  int edge_count    = solver.get_edge_count();

  EXPECT_EQ(vertex_count, instance_size);
  EXPECT_EQ(vertex_count, static_cast<int>(flaglist.size()));
  EXPECT_EQ(edge_count, static_cast<int>(edgelist.second.size()));

  long best_step = solver.get_best_step();
  EXPECT_GT(best_step, 1);

  std::chrono::milliseconds best_dur  = solver.get_best_duration();
  std::chrono::milliseconds total_dur = solver.get_total_duration();
  EXPECT_GE(best_dur.count(), 0);
  EXPECT_GE(total_dur.count(), best_dur.count());
}

// Execute verbose-only code lines, useful for coverage
TYPED_TEST(InterfaceTests, VerboseOutput){
  using SOLVER=TypeParam;

  std::ifstream file(filename, std::ios::in);
  SOLVER solver(
      file,
      cover_size,
      std::chrono::seconds(cutoff_time),
      true);

  solver.cover_LS(SOLVER::default_stats_printer);
}


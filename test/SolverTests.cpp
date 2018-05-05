#include "gtest/gtest.h"
#include "../LibNuMVC/numvc.hpp"
#include "../LibFastVC/fastvc.hpp"
#include "../ParallelSolverAdapter/ParallelSolverAdapter.hpp"

#include <fstream>
#include <algorithm>

using namespace libmvc;

template<typename SOLVER>
class SolverTests : public ::testing::Test { };

typedef ::testing::Types<NuMVC, FastVC, ParallelSolverAdapter<NuMVC>, ParallelSolverAdapter<FastVC>> SolverTypes;
TYPED_TEST_CASE(SolverTests, SolverTypes);

extern std::string filename;
extern int cover_size;
extern int instance_size;
extern int cutoff_time;

TYPED_TEST(SolverTests, SimpleVC){
  using SOLVER=TypeParam;

  std::ifstream file(filename, std::ios::in);
  SOLVER solver(file, cover_size, std::chrono::seconds(cutoff_time));
  solver.cover_LS();

  // let solver check solution
  EXPECT_EQ(solver.check_solution(), true);

  // blackblox check solution
  std::vector<int> solution(solver.get_independent_set());
  unsigned int expected = instance_size - cover_size;
  EXPECT_EQ(solution.size(), expected);
}

TYPED_TEST(SolverTests, EdgeListConstructor){
  using SOLVER=TypeParam;
  using EDGE  = std::pair<int,int>;
  std::ifstream file(filename, std::ios::in);
  std::pair<int, std::vector<EDGE>> edge_list;

  // avoid dependencies
  {
    SOLVER master(file, cover_size, std::chrono::seconds(cutoff_time));
    edge_list = master.get_instance_as_edgelist();
    EXPECT_EQ(edge_list.first, instance_size);
    EXPECT_GT(edge_list.second.size(), 0u);
  }

  SOLVER child(edge_list.second, edge_list.first,
      cover_size, std::chrono::seconds(cutoff_time));
  child.cover_LS();

  // let solver check solution
  EXPECT_EQ(child.check_solution(), true);

  // blackblox check solution
  auto solution(child.get_independent_set());
  unsigned int expected = instance_size - cover_size;
  EXPECT_EQ(solution.size(), expected);
 }

TYPED_TEST(SolverTests, InitialSolution){
  using SOLVER=TypeParam;

  int first_cover = std::max((cover_size + instance_size) / 2, cover_size);
  std::ifstream file(filename, std::ios::in);
  SOLVER master(file, first_cover, std::chrono::seconds(cutoff_time));
  master.cover_LS();

  // let solver check solution
  EXPECT_EQ(master.check_solution(), true);


  auto cover = master.get_cover_as_flaglist();
  EXPECT_EQ(cover.size(), static_cast<unsigned int>(instance_size));

  auto instance = master.get_instance_as_edgelist();

  // construct new sub-solver
  SOLVER child(instance.second, instance.first, cover_size,
      std::chrono::seconds(cutoff_time));
  child.cover_LS();

  EXPECT_EQ(child.check_solution(), true);

  // blackblox check solution
  auto solution(child.get_independent_set());
  unsigned int expected = instance_size - cover_size;
  EXPECT_EQ(solution.size(), expected);
}


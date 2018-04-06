#include "gtest/gtest.h"
#include "../LibNuMVC/numvc.hpp"
#include "../LibFastVC/fastvc.hpp"
#include "../ParallelSolverAdapter/ParallelSolverAdapter.hpp"

#include <fstream>

template<typename SOLVER>
class SolverTests : public ::testing::Test { };

typedef ::testing::Types<NuMVC, FastVC, ParallelSolverAdapter<NuMVC>, ParallelSolverAdapter<FastVC>> SolverTypes;
TYPED_TEST_CASE(SolverTests, SolverTypes);

TYPED_TEST(SolverTests, SimpleVC){
  using SOLVER=TypeParam;
  std::string filename("data/frb30-15-mis/frb30-15-1.mis");
  int cutoff_time   = 10;
  int cover_size    = 420;
  int instance_size = 450;

  std::ifstream file(filename, std::ios::in);
  SOLVER solver(file, cover_size, std::chrono::seconds(cutoff_time));
  solver.cover_LS();

  // let solver check solution
  EXPECT_EQ(solver.check_solution(), true);

  // blackblox check solution
  std::vector<int> solution(solver.get_independent_set());
  EXPECT_EQ(solution.size(), instance_size - cover_size);
}


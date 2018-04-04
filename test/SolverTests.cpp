#include "gtest/gtest.h"
#include "../LibNuMVC/numvc.hpp"
#include "../LibFastVC/fastvc.hpp"

template<typename SOLVER>
class SolverTests : public ::testing::Test { };

typedef ::testing::Types<NuMVC, FastVC> SolverTypes;
TYPED_TEST_CASE(SolverTests, SolverTypes);

TYPED_TEST(SolverTests, SimpleVC){

}

#include "gtest/gtest.h"
#include <string>
#include <fstream>
#include <stdexcept>

using ::testing::UnitTest;

// test instance
std::string filename;
int cutoff_time;
int cover_size;
int instance_size;

int main(int argc, char * argv[]){
  ::testing::InitGoogleTest(&argc, argv);

  // parse arguments
  if(argc == 5){
    filename      = argv[1];
    cover_size    = std::stoi(argv[2]);
    instance_size = std::stoi(argv[3]);
    cutoff_time   = std::stoi(argv[4]);
  } else {
    std::cout << "usage: <filename> <coversize> <instancesize> <cutoff_time[s]>"
              << std::endl;
    return 1;
  }

  // check if file exists
  std::ifstream file(filename, std::ios::in);
  if(!file.good()){
    throw std::runtime_error("File not readable");
  }
  file.close();

  int ret = RUN_ALL_TESTS();
  return ret;
}

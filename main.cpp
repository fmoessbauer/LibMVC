#include "fastvc.hpp"
#include <fstream>
#include <vector>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[])
{
  char title[1024];

  int optimal_size;
  int cutoff_time;

  if(argc != 4){
    cout << "usage: " << argv[0] << " graph-file optimal-size cutoff-time(s)"
         << endl;
    return 1;
  }

  sscanf(argv[2],"%d",&optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
  sscanf(argv[3],"%d",&cutoff_time);

  sprintf(title,"@title file: %s        optimal_size: %d",argv[1], optimal_size);

  std::string filename(argv[1]);
  std::ifstream file(filename, std::ios::in);
  FastVC solver(file, optimal_size, std::chrono::seconds(cutoff_time));

  cout<<"c Improved NuMVC Local Search Solver"<<endl;

  cout<<"c Start local search..."<<endl;

  solver.cover_LS([&solver](){
    // print stats
    auto time_ms = chrono::duration_cast<chrono::milliseconds>(solver.get_best_duration());
    std::cout << "Better MVC found.\tSize: "
              << solver.get_best_cover_size()
              << "\tTime: " << std::fixed << std::setw(4) << std::setprecision(4)
              << time_ms.count() << "ms" << std::endl;
    return false;
  });
  auto solution = solver.get_cover();

    //check solution
    if(solver.check_solution())
    {
      cout<<"c Best found vertex cover size = "<<solution.size()<<endl;
      solver.print_solution();
      auto solver_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(solver.get_best_duration());
      auto total_time_ms  = std::chrono::duration_cast<std::chrono::milliseconds>(solver.get_total_duration());
      double performance = static_cast<double>(solver.get_best_step()) / solver_time_ms.count() / 1000.0;
      cout<<"c searchSteps = "<<solver.get_best_step()<<endl;
      cout<<"c solveTime = "<<solver_time_ms.count()<<"ms"<<endl;
      cout<<"c totalTime = "<<total_time_ms.count()<<"ms"<<endl;
      cout<<"c performance = "<<performance<< "MT/s"<<endl;
    }
    else
    {
      cout<<"the solution is wrong."<<endl;
    }

  return 0;
}

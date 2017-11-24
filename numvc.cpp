#include "tsewf.hpp"
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	char title[1024];
	
  int optimal_size;
  int cutoff_time;

  sscanf(argv[2],"%d",&optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
  sscanf(argv[3],"%d",&cutoff_time);

	sprintf(title,"@title file: %s        optimal_size: %d",argv[1], optimal_size);

  std::string filename(argv[1]);
  std::ifstream file(filename, std::ios::in);
  TSEWF solver(file, optimal_size, cutoff_time);

  cout<<"c Improved NuMVC Local Search Solver"<<endl;

  std::vector<int> solution;

  cout<<"c Start local search..."<<endl;
  solver.cover_LS();
  solution = std::move(solver.get_independent_set());
		
		//check solution
		if(solver.check_solution())
		{
			cout<<"c Best found independent set size = "<<solution.size()<<endl;
			solver.print_solution();
			cout<<"c searchSteps = "<<solver.best_step<<endl;
			cout<<"c solveTime = "<<solver.best_comp_time<<endl;
			cout<<"c performance = "<<(solver.best_step/solver.best_comp_time)/1000000<< "MT/s"<<endl;
		}
		else
		{
			cout<<"the solution is wrong."<<endl;
		}

	return 0;
}

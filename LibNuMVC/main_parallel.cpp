#include "numvc.hpp"
#include <fstream>
#include <chrono>
#include <thread>

using namespace std;

void start_solver(NuMVC * const solver){
  cout<<"Start Solver"<<endl;
  solver->cover_LS(NuMVC::default_stats_printer);
}

int main(int argc, char* argv[])
{
    char title[1024];

    int optimal_size;
    int cutoff_time;
    auto num_solvers = std::max(std::thread::hardware_concurrency()/2, 1u);

    sscanf(argv[2],"%d",&optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
    sscanf(argv[3],"%d",&cutoff_time);

    sprintf(title,"@title file: %s        optimal_size: %d",argv[1], optimal_size);

    std::string filename(argv[1]);
    std::ifstream file(filename, std::ios::in);
    NuMVC solver(file, optimal_size, std::chrono::seconds(3), true);
    file.close();

    // duplicate initial solver
    std::vector<NuMVC> solvers(num_solvers, solver);
    for(unsigned int i=0;i<num_solvers; ++i){
      auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      solvers[i].set_random_seed(seed+i);
    }

    cout<<"c Parallel NuMVC Local Search Solver"<<endl;
    cout<<"c Use "<<num_solvers<<" solvers"<<endl;

    std::vector<std::thread> threads;
    for(unsigned int i=0; i<num_solvers; ++i){
      threads.emplace_back(start_solver, &(solvers[i]));
    }
    for(auto & t : threads){
      t.join();
    }

    // gather results, calculate intersection
    auto intersection = solvers[0].get_cover();
    for(unsigned int i=1; i<num_solvers; ++i){
      const auto & cover = solvers[i].get_cover();
      std::vector<int> target;
      std::set_intersection(
          intersection.begin(), intersection.end(),
          cover.begin(), cover.end(),
          std::back_inserter(target));
      intersection = std::move(target);
    }
    std::cout << "Intersection size: " << intersection.size() << std::endl;
    const auto & instance = solver.get_instance_as_edgelist();
    auto edgelist(instance.second);
    edgelist.erase(std::remove_if(
        edgelist.begin(), edgelist.end(),
        [&intersection](const NuMVC::Edge & x){
          return (std::binary_search(intersection.begin(), intersection.end(), x.first)
               || std::binary_search(intersection.begin(), intersection.end(), x.second));
          }), edgelist.end());

    int inst_size  = instance.first - intersection.size();
    int best_cover = optimal_size - inst_size;
    std::cout << "Final Instance size: " << inst_size << std::endl;

    NuMVC final_solver(edgelist, instance.first, best_cover, std::chrono::seconds(cutoff_time), true);
    final_solver.cover_LS(NuMVC::default_stats_printer);
    // get vertex cover
    auto v_s = final_solver.get_cover();
    std::vector<int> final_cover;
    std::set_union(
        intersection.begin(), intersection.end(),
        v_s.begin(), v_s.end(),
        std::back_inserter(final_cover)
        );

    //check solution
    cout<<"c Best found vertex cover size = "<<final_cover.size()<<endl;
    return 0;
}

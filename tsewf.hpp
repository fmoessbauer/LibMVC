/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/


/************************************************
** Date:	2011.7.1
** TSEWF (Two Stage Exchange and Weighting with Forgetting)
** Author: Shaowei Cai, shaowei_cai@126.com
**       School of EECS, Peking University
**       Beijing, China
**
** Date:	2011.10.28
** Modify: Shaowei Cai
** use dynamic memory for v_adj[][] and v_edge[][], tidy codes.
**
** Date: 2017.11.17
** Modify: Felix Moessbauer
** Object oriented design, C++11 support
************************************************/

#ifndef TSEWF_INCLUDED
#define TSEWF_INCLUDED

#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>

namespace chrono  = std::chrono;

class TSEWF {
private:
  using timepoint_t = chrono::time_point<
                        chrono::system_clock>;
private:

  bool verbose = false;

  timepoint_t start, finish;
  timepoint_t start_time;

  struct Edge{
    int v1;
    int v2;
  };

  const int try_step = 600000;

  /*parameters of algorithm*/
  long long  max_steps;      //step limit
  int        cutoff_time;    //time limit
  long long  step;
  int        optimal_size;   //terminate the algorithm before step limit if it finds a vertex cover of optimal_size

  /*parameters of the instance*/
public:
  int    v_num;//|V|: 1...v
  int    e_num;//|E|: 0...e-1
private:
  /*structures about edge*/
  std::vector<Edge> edge;
  std::vector<int>  edge_weight;

  /*structures about vertex*/
  std::vector<int>       dscore;     //dscore of v
  std::vector<long long> time_stamp;

  //from vertex to it's edges and neighbors
  std::vector<int*> v_edges;      //edges related to v, v_edges[i][k] means vertex v_i's k_th edge
  std::vector<int*> v_adj;        //v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
  std::vector<int>  v_edge_count; //amount of edges (neighbors) related to v


  /* structures about solution */
  //current candidate solution
  int    c_size;            //cardinality of C
  std::vector<int> v_in_c;      //a flag indicates whether a vertex is in C
  std::vector<int> remove_cand; //remove candidates, an array consists of only vertices in C, not including tabu_remove
  std::vector<int> index_in_remove_cand;
  int    remove_cand_size;

  //best solution found
  int              best_c_size;
  std::vector<int> best_v_in_c; //a flag indicates whether a vertex is in best solution

  int best_vertex_improvement;
  int best_count;
  std::vector<int> best_array;

public:
  double  best_comp_time;
  long    best_step;
private:
  //uncovered edge stack
  std::vector<int> uncov_stack;          //store the uncov edge number
  int    uncov_stack_fill_pointer;
  std::vector<int> index_in_uncov_stack; //which position is an edge in the uncov_stack
  std::vector<int> v_edge_count_tmp;


  //CC and taboo
  std::vector<int> conf_change;
  int   tabu_remove=0;

  //smooth
  int    ave_weight=1;
  int    delta_total_weight=0;
  int    threshold;
  float  p_scale=0.3;//w=w*p

public:
  template<typename Is>
  TSEWF(
      Is & str,
      int optimal_size,
      int cutoff_time,
      bool verbose = false)
  : verbose(verbose),
    optimal_size(optimal_size),
    cutoff_time(cutoff_time)
  {
    build_instance(str);
    init_sol();
  }

  ~TSEWF(){
    free_memory();
  }

private:
  void init_internal(int num_vertices, int num_edges){
    ++num_vertices;
    ++num_edges;

    threshold = static_cast<int>(0.5*num_vertices);
    srand(std::time(0));

    edge.resize(num_edges);
    edge_weight.resize(num_edges);
    dscore.resize(num_vertices);
    time_stamp.resize(num_vertices);
    v_edges.resize(num_vertices);
    v_adj.resize(num_vertices);
    v_edge_count.resize(num_vertices);

    v_in_c.resize(num_vertices);
    remove_cand.resize(num_vertices);
    index_in_remove_cand.resize(num_vertices);

    best_v_in_c.resize(num_vertices);

    //uncovered edge stack
    uncov_stack.resize(num_edges);          //store the uncov edge number
    index_in_uncov_stack.resize(num_edges); //which position is an edge in the uncov_stack
    v_edge_count_tmp.resize(num_vertices);

    //CC and taboo
    conf_change.resize(num_vertices);
    best_array.resize(num_vertices);
  }

  // copy v_in_c to best_v_in_c
  void update_best_sol()
  {
    int i;

    for (i=1;i<=v_num;i++)
    {
      best_v_in_c[i] = v_in_c[i];
    }

    best_c_size = c_size;
    finish = chrono::system_clock::now();
    auto diff_in_ms = std::chrono::duration_cast<
                        std::chrono::milliseconds>(finish - start);
    best_comp_time = static_cast<double>(diff_in_ms.count()) / 1000.0;
    best_step = step;
  }

private:
  template<typename Is>
  int build_instance(Is & str)
  {
    char line[1024];
    char tempstr1[10];
    char tempstr2[10];
    int  v,e;

    char  tmp;
    int   v1,v2;

    /*** build problem data structures of the instance ***/
    str.getline(line,1024);
    while (line[0] != 'p')
      str.getline(line,1024);
    sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
    init_internal(v_num, e_num);

    /* read edges and compute v_edge_count */
    for (v=1; v<=v_num; v++)
      v_edge_count[v] = 0;

    for (e=0; e<e_num; e++)
    {
      str>>tmp>>v1>>v2;
      v_edge_count[v1]++;
      v_edge_count[v2]++;

      edge[e].v1 = v1;
      edge[e].v2 = v2;
    }

    /* build v_adj and v_edges arrays */
    for (v=1; v<=v_num; v++)
    {
      v_adj[v] = new int[v_edge_count[v]];
      v_edges[v] = new int[v_edge_count[v]];
    }


    for(v=1; v<=v_num; v++)
      v_edge_count_tmp[v]=0;
    for (e=0; e<e_num; e++)
    {

      v1=edge[e].v1;
      v2=edge[e].v2;

      v_edges[v1][v_edge_count_tmp[v1]] = e;
      v_edges[v2][v_edge_count_tmp[v2]] = e;

      v_adj[v1][v_edge_count_tmp[v1]] = v2;
      v_adj[v2][v_edge_count_tmp[v2]] = v1;

      v_edge_count_tmp[v1]++;
      v_edge_count_tmp[v2]++;
    }

    return 1;
  }

  void free_memory()
  {
    for (int v=1; v<=v_num; v++)
    {
      delete[] v_adj[v];
      delete[] v_edges[v];
      v_adj[v]   = nullptr;
      v_edges[v] = nullptr;
    }
  }

  void reset_remove_cand()
  {
    int v,j;
    j=0;
    for (v=1;v<=v_num;v++)
    {
      if(v_in_c[v]==1)// && v!=tabu_remove)
      {
        remove_cand[j] = v;
        index_in_remove_cand[v]=j;
        j++;
      }
      else index_in_remove_cand[v]=0;
    }

    remove_cand_size = j;
  }

  // kick out the worst vertex in currennt cover
  void update_target_size()
  {
    c_size--;

    int v;
    int max_improvement;
    int max_vertex;//vertex with the highest improvement in C

    max_improvement=-200000000;
    for (int v=1; v<=v_num; ++v)
    {
      if(v_in_c[v]==0)continue;
      if (dscore[v]>max_improvement)
      {
        max_improvement = dscore[v];
        max_vertex = v;
      }
    }
    remove(max_vertex);

    reset_remove_cand();
  }

  inline void uncover(int e)
  {
    index_in_uncov_stack[e] = uncov_stack_fill_pointer;
    uncov_stack[uncov_stack_fill_pointer++] = e;
  }


  inline void cover(int e)
  {
    int index,last_uncov_edge;

    //since the edge is satisfied, its position can be reused to store the last_uncov_edge
    last_uncov_edge = uncov_stack[--uncov_stack_fill_pointer];
    index = index_in_uncov_stack[e];
    uncov_stack[index] = last_uncov_edge;
    index_in_uncov_stack[last_uncov_edge] = index;
  }


  void init_sol()
  {
    int i,v,e;

    /*** build solution data structures of the instance ***/
    //init vertex cover
    for (v=1; v<=v_num; v++)
    {
      v_in_c[v] = 0;
      dscore[v] = 0;

      conf_change[v] = 1;
      time_stamp[v]= 0; // to break ties
    }

    for (e=0; e<e_num; e++)
    {
      edge_weight[e] = 1;
      dscore[edge[e].v1]+=edge_weight[e];
      dscore[edge[e].v2]+=edge_weight[e];
    }

    //init uncovered edge stack and cover_vertrex_count_of_edge array
    uncov_stack_fill_pointer = 0;
    for (e=0; e<e_num; e++)
      uncover(e);



    for (i=0; uncov_stack_fill_pointer>0; )
    {
      best_vertex_improvement = 0;
      best_count = 0;
      for (v=1; v<=v_num; ++v)
      {
        if(v_in_c[v]==1)continue;

        if (dscore[v]>best_vertex_improvement)
        {
          best_vertex_improvement = dscore[v];
          best_array[0] = v;
          best_count = 1;
        }
        else if (dscore[v]==best_vertex_improvement)
        {
          best_array[best_count] = v;
          best_count++;
        }
      }

      if(best_count>0)
      {
        add(best_array[rand()%best_count]);
        ++i;
      }
    }

    if(verbose){
      std::cout << "Initial cover size: " << i << std::endl;
    }

    c_size = i;

    update_best_sol();

    reset_remove_cand();

  }

  // add a vertex to current cover
  void add(int v)
  {
    v_in_c[v] = 1;
    dscore[v] = -dscore[v];

    int i,e,n;

    int edge_count = v_edge_count[v];

    for (i=0; i<edge_count; ++i)
    {
      e = v_edges[v][i];// v's i'th edge
      n = v_adj[v][i];//v's i'th neighbor

      if (v_in_c[n]==0)//this adj isn't in cover set
      {
        dscore[n] -= edge_weight[e];
        conf_change[n] = 1;

        cover(e);
      }
      else
      {
        dscore[n] += edge_weight[e];
      }
    }

  }

  void remove(int v)
  {
    v_in_c[v] = 0;
    dscore[v] = -dscore[v];
    conf_change[v] = 0;

    int i,e,n;

    int edge_count = v_edge_count[v];
    for (i=0; i<edge_count; ++i)
    {
      e = v_edges[v][i];
      n = v_adj[v][i];

      if (v_in_c[n]==0)//this adj isn't in cover set
      {
        dscore[n] += edge_weight[e];
        conf_change[n] = 1;

        uncover(e);
      }
      else
      {
        dscore[n] -= edge_weight[e];
      }
    }

  }


  inline unsigned int rdrand32 ()
  {
    register unsigned int rand;
    asm volatile ("rdrand %0"
      : "=r" (rand)
    );
    return rand;
  }

  void forget_edge_weights()
  {
    int v,e;
    int new_total_weight=0;

    for(v=1; v<=v_num; v++)
      dscore[v]=0;

    //scale_ave=ave_weight*q_scale;
    for (e = 0; e<e_num; e++)
    {
      edge_weight[e] = edge_weight[e]*p_scale;

      new_total_weight+=edge_weight[e];

      //update dscore
      if (v_in_c[edge[e].v1]+v_in_c[edge[e].v2]==0){
        dscore[edge[e].v1]+=edge_weight[e];
        dscore[edge[e].v2]+=edge_weight[e];
        }
      else if(v_in_c[edge[e].v1]+v_in_c[edge[e].v2]==1){
        if(v_in_c[edge[e].v1]==1)dscore[edge[e].v1]-=edge_weight[e];
        else  dscore[edge[e].v2]-=edge_weight[e];
      }
    }
    ave_weight=new_total_weight/e_num;

  }


  void update_edge_weight()
  {
    int i,e;
    for(i=0; i<uncov_stack_fill_pointer; ++i)
    {
      e = uncov_stack[i];

      edge_weight[e]+= 1;
      dscore[edge[e].v1] += 1;
      dscore[edge[e].v2] += 1;
    }


    delta_total_weight += uncov_stack_fill_pointer;
    if(delta_total_weight >= e_num)
    {
      ave_weight += 1;
      delta_total_weight -= e_num;
    }

    //smooth weights
    if(ave_weight >= threshold)
    {
      forget_edge_weights();
    }

  }

public:
  std::vector<int> get_independent_set(){
    std::vector<int> solution;
    solution.reserve(v_num-best_c_size);
    for (int i=1; i<=v_num; i++)
     {
       if (best_v_in_c[i]!=1)//output max independent set
         solution.push_back(i);
     }
    return solution;
  }

  void cover_LS()
  {
    int    best_cov_v;    //the vertex of the highest dscore in C
    int    best_add_v;
    int    e,v1,v2;
    int    i,v;

    step  = 1;
    start = chrono::system_clock::now();

    while(1)// wihin cutoff_time
    //while(step<=max_steps)
    {
      /* ### if there is no uncovered edge ### */
      if (uncov_stack_fill_pointer == 0)
      {
        update_best_sol();      // C* := C

        if(verbose){
          std::cout << "Better MVC found.\tSize: " << v_num-best_c_size
                    << "\tTime: " << std::fixed << std::setw(4) << std::setprecision(4)
                    << best_comp_time << "s" << std::endl;
        }

        if (c_size==optimal_size)
          return;

        update_target_size();    // remove a vertex with the highest dscore from C;

        continue;
      }

      /* monitor status */

      if(step % try_step==0)
      {
        finish = chrono::system_clock::now();
        auto elapsed_ms  = chrono::duration_cast<chrono::milliseconds>(finish - start_time);
        double elap_time = static_cast<double>(elapsed_ms.count()) / 1000;
        if(elap_time >= cutoff_time) return;
      }

      /* choose a vertex u in C with the highest dscore,
        breaking ties in favor of the oldest one; */
      best_cov_v = remove_cand[0];
      for (i=1; i<remove_cand_size; ++i)
      {
        v = remove_cand[i];
        if(v==tabu_remove)
          continue;
        if( dscore[v] < dscore[best_cov_v])
          continue;
        else if( dscore[v] > dscore[best_cov_v] )
          best_cov_v = v;
        else if (time_stamp[v] < time_stamp[best_cov_v])
          best_cov_v = v;
      }

      /* C := C\{u}, confChange(u) := 0 and confChange(z) := 1 for each z in N(u); */
      remove(best_cov_v);


      /* choose an uncovered edge e randomly; */
      //e = uncov_stack[rand()%uncov_stack_fill_pointer];
      e = uncov_stack[rdrand32()%uncov_stack_fill_pointer];

      /* choose a vertex v in e such that confChange(v) = 1 with higher dscore,
        breaking ties in favor of the older one; */
      v1 = edge[e].v1;
      v2 = edge[e].v2;

      if(conf_change[v1]==0 ) best_add_v=v2;
      else if(conf_change[v2]==0) best_add_v=v1;
      else
      {
        if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
          best_add_v=v1;
        else
          best_add_v=v2;
      }


      /* C := C plus {v}, confChange(z) := 1 for each z in N(v); */
      add(best_add_v);

      // update remove candidate index. speed optimization.
      int index = index_in_remove_cand[best_cov_v];
      index_in_remove_cand[best_cov_v] = 0;

      remove_cand[index] = best_add_v;
      index_in_remove_cand[best_add_v] = index;

      // update timestamp used to determine the age of each virtex
      time_stamp[best_add_v] = step;
      time_stamp[best_cov_v] = step;

      tabu_remove = best_add_v;

      /* w(e) := w(e) + 1 for each uncovered edge e; */
      /* if w >= y then w(e) := [p*w(e)] for each edge e; */
      update_edge_weight();

      step++;
    }
    return;
  }

/*On solution*/
  void print_solution()
  {
    int mis_vertex_count=0;

    for (int i=1; i<=v_num; i++)
    {
      if (best_v_in_c[i]!=1)
        mis_vertex_count++;
    }

    if(mis_vertex_count+best_c_size!=v_num){
      std::cout << "The size of independent set + the size of "
                << "vertex cover is not equal to |V(G)|!" << std::endl;
    }

    std::cout << "c Best found independent set size = "
              << mis_vertex_count << std::endl;
    std::cout << "c The following output is the found independent set."
              << std::endl;

    for (int i=1; i<=v_num; i++)
    {
      if (best_v_in_c[i]!=1)//output max independent set
        std::cout << i << "  ";
    }
    std::cout << std::endl;
  }

  //check whether the solution found is a proper solution
  bool check_solution()
  {
    int e;
    for(e=0; e<e_num; ++e)
    {
      if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1)
      {
        std::cout << "uncovered edge " << e << std::endl;
        return false;
      }
    }

    return true;
  }

};

#endif

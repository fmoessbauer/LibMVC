/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/


/************************************************
** Date:  2011.7.1
** NuMVC (Two Stage Exchange and Weighting with Forgetting)
** Author: Shaowei Cai, shaowei_cai@126.com
**       School of EECS, Peking University
**       Beijing, China
**
** Date:  2011.10.28
** Modify: Shaowei Cai
** use dynamic memory for v_adj[][] and v_edge[][], tidy codes.
**
** Date: 2017.11.17
** Modify: Felix Moessbauer
** Object oriented design, C++11 support
** Use heap for construction of initial solution
************************************************/

#ifndef NuMVC_INCLUDED
#define NuMVC_INCLUDED

#include <chrono>
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <functional>
#include <algorithm>
#include <random>

#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "indexed_heap.hpp"

class NuMVC {
private:
    using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
    using duration_ms = std::chrono::milliseconds;

private:

    bool verbose = false;

    timepoint_t start, finish;

    struct Edge {
        int v1;
        int v2;
    };

    const int try_step = 600000;

    /*parameters of algorithm*/
    long long   max_steps;      //step limit
    duration_ms cutoff_time;    //time limit
    long long   step;
    int         optimal_size;   //terminate the algorithm before step limit if it finds a vertex cover of optimal_size

    /*parameters of the instance*/
    int    v_num;//|V|: 1...v
    int    e_num;//|E|: 0...e-1

    /*structures about edge*/
    std::vector<Edge>    edge;
    std::vector<int>     edge_weight;
    const int default_edge_weight = 1;

    /*structures about vertex*/
    std::vector<int>       dscore;     //dscore of v
    std::vector<long long> time_stamp;
    int                    best_cov_v;

    //from vertex to it's edges and neighbors
    std::vector<int> v_beg_idx; // v_edges and v_adj is flattened 2-d array, hence store indices
    std::vector<int> v_edges;   //edges related to v, v_edges[i][k] means vertex v_i's k_th edge
    std::vector<int> v_adj;     //v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
    std::vector<int> v_degree;  //amount of edges (neighbors) related to v


    /* structures about solution */
    //current candidate solution
    int               c_size;      //cardinality of C
    std::vector<bool> v_in_c;      //a flag indicates whether a vertex is in C
    std::vector<int>  remove_cand; //remove candidates, an array consists of only vertices in C, not including tabu_remove
    std::vector<int>  index_in_remove_cand;
    int               remove_cand_size;

    //best solution found
    int               best_c_size;
    std::vector<bool> best_v_in_c; //a flag indicates whether a vertex is in best solution

    duration_ms best_comp_time;
    long        best_step;

    //uncovered edge stack
    std::vector<int> uncov_stack;          //store the uncov edge number
    int    uncov_stack_fill_pointer;
    std::vector<int> index_in_uncov_stack; //which position is an edge in the uncov_stack
    std::vector<int> v_degree_tmp;


    //CC and taboo
    std::vector<int> conf_change;
    int              tabu_remove=0;

    // priority queue for initial construction (max element at first index)
    // build heap by comparing dscores
    std::function<bool (const int &, const int &)>
          dscore_cmp = [&](const int & a, const int & b)
    { return (dscore[a] < dscore[b]); };

    using heap_t = Indexed_Heap<int, decltype(dscore_cmp)>;
    heap_t v_heap;

    //smooth
    static constexpr float  p_scale=0.3;//w=w*p
    int                     delta_total_weight=0;
    int                     ave_weight=1;
    int                     threshold;

    std::mt19937 mt_rand;

public:
    /**
     * Construct solver instance by importing a graph in DIMACS format
     */
    template<
        typename Is,
        typename Duration>
    NuMVC(
        /// Input Stream with graph in DIMACS format
        Is & str,
        /// Size of optimal vertex cover (set to 0 if not known)
        int optimal_size,
        /// Stop calculation after this duration (chrono duration)
        Duration cutoff_time,
        /// Print messages during calculation
        bool verbose = false,
        /// seed for random number generator
        unsigned int rnd_seed =
            std::chrono::high_resolution_clock::now().time_since_epoch().count())
        : verbose(verbose),
          cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
          optimal_size(optimal_size),
          v_heap(dscore_cmp)
    {
        mt_rand.seed(rnd_seed);
        build_instance(str);
    }

private:
    void init_internal(int num_vertices, int num_edges)
    {
        ++num_vertices;
        ++num_edges;

        threshold = static_cast<int>(0.5*num_vertices);

        edge.resize(num_edges);
        edge_weight.resize(num_edges, default_edge_weight);
        dscore.resize(num_vertices);
        time_stamp.resize(num_vertices);
        v_degree.resize(num_vertices);

        v_in_c.resize(num_vertices);
        remove_cand.resize(num_vertices);
        index_in_remove_cand.resize(num_vertices);

        best_v_in_c.resize(num_vertices);

        //uncovered edge stack
        uncov_stack.resize(num_edges);          //store the uncov edge number
        index_in_uncov_stack.resize(num_edges); //which position is an edge in the uncov_stack
        v_degree_tmp.resize(num_vertices);

        v_heap.resize(num_vertices);

        //CC and taboo
        conf_change.resize(num_vertices, 1);
    }

    // copy v_in_c to best_v_in_c
    void update_best_sol()
    {
        int i;

        for (i=1; i<=v_num; i++) {
            best_v_in_c[i] = v_in_c[i];
        }

        best_c_size = c_size;
        finish = std::chrono::system_clock::now();
        best_comp_time = std::chrono::duration_cast<duration_ms>(finish - start);
        best_step = step;
    }

    void update_best_cov_v()
    {
        int i,v;
        best_cov_v = remove_cand[0];
        for (i=1; i<remove_cand_size; ++i) {
            v = remove_cand[i];
            if(v==tabu_remove) continue;
            if( dscore[v] < dscore[best_cov_v])
                continue;
            else if( dscore[v]> dscore[best_cov_v] )
                best_cov_v = v;
            else if (time_stamp[v]<time_stamp[best_cov_v])
                best_cov_v = v;
        }
    }

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

        /* read edges and compute v_degree */
        for (v=1; v<=v_num; v++)
            v_degree[v] = 0;

        for (e=0; e<e_num; e++) {
            str>>tmp>>v1>>v2;
            v_degree[v1]++;
            v_degree[v2]++;

            edge[e].v1 = v1;
            edge[e].v2 = v2;
        }

        /* indices are the partial sums */
        v_beg_idx.reserve(e_num+1);
        v_beg_idx.push_back(0); // shift by one as partial sum calculates end_index
        std::partial_sum(v_degree.begin(), v_degree.end(), std::back_inserter(v_beg_idx));

        v_edges.resize(v_beg_idx.back());
        v_adj.resize(v_beg_idx.back());

        std::vector<int> v_degree_tmp(v_num + 1);
        for (e=0; e<e_num; ++e) {
            v1=edge[e].v1;
            v2=edge[e].v2;

            v_edges[v_beg_idx.at(v1) + v_degree_tmp[v1]] = e;
            v_edges[v_beg_idx.at(v2) + v_degree_tmp[v2]] = e;

            v_adj[v_beg_idx[v1] + v_degree_tmp[v1]] = v2;
            v_adj[v_beg_idx[v2] + v_degree_tmp[v2]] = v1;

            ++(v_degree_tmp[v1]);
            ++(v_degree_tmp[v2]);
        }

        return 1;
    }

    void reset_remove_cand()
    {
        int v,j;
        j=0;
        for (v=1; v<=v_num; ++v) {
            if(v_in_c[v]) { // && v!=tabu_remove)
                remove_cand[j] = v;
                index_in_remove_cand[v]=j;
                ++j;
            } else index_in_remove_cand[v]=0;
        }

        remove_cand_size = j;
    }

    // kick out the worst vertex in current cover
    void update_target_size()
    {
        --c_size;

        int max_improvement = std::numeric_limits<int>::min();
        int max_vertex      = 0;//vertex with the highest improvement in C

        for (int v=1; v<=v_num; ++v) {
            if(!v_in_c[v])continue;
            if (dscore[v]>max_improvement) {
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
        start = std::chrono::system_clock::now();

        /*** build solution data structures of the instance ***/
        //init vertex cover
        //conf_change = 1 (already set in resize call)
        for (int e=0; e<e_num; ++e) {
            const auto weight = edge_weight[e];
            dscore[edge[e].v1]+=weight;
            dscore[edge[e].v2]+=weight;
        }

        //init uncovered edge stack and cover_vertrex_count_of_edge array

#if 0
        uncov_stack_fill_pointer = 0;
        for (int e=0; e<e_num; ++e)
            uncover(e);
#else
        // tweak to uncover all
        std::iota(index_in_uncov_stack.begin(),
                  index_in_uncov_stack.end(),
                  0);
        memcpy(uncov_stack.data(), index_in_uncov_stack.data(), e_num);
        uncov_stack_fill_pointer = e_num;
#endif

        for(int v=1; v<=v_num; ++v) {
            v_heap.push(v);
        }

        int i = 0;
        while(uncov_stack_fill_pointer>0) {
            int best_v = v_heap.top();
            v_heap.pop();
            if(dscore[best_v]>0) {
                add_init(best_v);
                ++i;
            }
        }

        c_size = i;
        update_best_sol();

        finish = std::chrono::system_clock::now();
        auto init_sol_time = std::chrono::duration_cast<duration_ms>(finish-start);
        if(verbose) {
            std::cout << "Initial solution size: " << c_size << std::endl;
            std::cout << "Initial solution time: " << init_sol_time.count()
                      << "ms" << std::endl;
        }

        reset_remove_cand();
        update_best_cov_v();
    }

    // add a vertex to current cover
    inline void add(int v)
    {
        v_in_c[v] = true;
        dscore[v] = -dscore[v];

        int i,e,n;

        int edge_count = v_degree[v];

        for (i=0; i<edge_count; ++i) {
            e = v_edges[v_beg_idx[v]+i]; // v's i'th edge
            n = v_adj[v_beg_idx[v]+i];   //v's i'th neighbor

            if (!v_in_c[n]) { //this adj isn't in cover set
                dscore[n] -= edge_weight[e];
                conf_change[n] = 1;

                cover(e);
            } else {
                dscore[n] += edge_weight[e];
            }
        }
    }

    void add_init(int v)
    {
        v_in_c[v] = true;
        dscore[v] *= (-1);

        int i,e,n;

        const int & degree = v_degree[v];

        for (i=0; i<degree ; ++i) {
            e = v_edges[v_beg_idx[v]+i]; // v's i'th edge
            n = v_adj[v_beg_idx[v]+i];   //v's i'th neighbor

            if (!v_in_c[n]) { //this adj isn't in cover set
                bool inheap = v_heap.count(n) > 0;
                if(inheap){
                  v_heap.erase(v_heap[n]);
                }
                dscore[n] -= edge_weight[e];
                conf_change[n] = 1;
                if(inheap){
                  v_heap.push(n);
                }
                cover(e);
            } else {
                dscore[n] += edge_weight[e];
            }
        }
    }

    inline void remove(int v)
    {
        v_in_c[v] = false;
        dscore[v] = -dscore[v];
        conf_change[v] = 0;

        int i,e,n;

        int edge_count = v_degree[v];
        for (i=0; i<edge_count; ++i) {
            e = v_edges[v_beg_idx[v] + i];
            n = v_adj[v_beg_idx[v] + i];

            if (!v_in_c[n]) { //this adj isn't in cover set
                dscore[n] += edge_weight[e];
                conf_change[n] = 1;

                uncover(e);
            } else {
                dscore[n] -= edge_weight[e];
            }
        }
    }

    void forget_edge_weights()
    {
        int v,e;
        int new_total_weight=0;

        for(v=1; v<=v_num; v++)
            dscore[v]=0;

        //scale_ave=ave_weight*q_scale;
        for (e = 0; e<e_num; e++) {
            edge_weight[e] = edge_weight[e]*p_scale;

            new_total_weight+=edge_weight[e];

            //update dscore
            if (!(v_in_c[edge[e].v1] || v_in_c[edge[e].v2])) {
                dscore[edge[e].v1]+=edge_weight[e];
                dscore[edge[e].v2]+=edge_weight[e];
            } else if(v_in_c[edge[e].v1] != v_in_c[edge[e].v2]) {
                if(v_in_c[edge[e].v1])dscore[edge[e].v1]-=edge_weight[e];
                else  dscore[edge[e].v2]-=edge_weight[e];
            }
        }
        ave_weight=new_total_weight/e_num;

    }


    void update_edge_weight()
    {
        int i,e;
        for(i=0; i<uncov_stack_fill_pointer; ++i) {
            e = uncov_stack[i];

            edge_weight[e]+= 1;
            dscore[edge[e].v1] += 1;
            dscore[edge[e].v2] += 1;
        }


        delta_total_weight += uncov_stack_fill_pointer;
        if(delta_total_weight >= e_num) {
            ave_weight += 1;
            delta_total_weight -= e_num;
        }

        //smooth weights
        if(ave_weight >= threshold) {
            forget_edge_weights();
        }

    }

public:
    /**
     * calculate minimum vertex cover
     */
    void cover_LS()
    {
        cover_LS(nullptr);
    }

    /**
     * calculate minimum vertex cover and call callback after
     * each iteration. If callback returns true, stop calculation.
     */
    void cover_LS(const std::function<bool (const NuMVC&)> & callback_on_update)
    {
        int    best_cov_v;    //the vertex of the highest dscore in C
        int    best_add_v;
        int    e,v1,v2;
        int    i,v;

        // No cover given, calculate inital cover
        if(c_size == 0){
            init_sol();
        }

        step  = 1;
        start = std::chrono::system_clock::now();

        while(true){
            /* ### if there is no uncovered edge ### */
            if (uncov_stack_fill_pointer == 0) {
                update_best_sol();      // C* := C
                if(callback_on_update != nullptr && callback_on_update(*this)) {
                    return;
                }

                if (c_size==optimal_size) return;

                update_target_size();    // remove a vertex with the highest dscore from C;
                continue;
            }

            /* monitor status */

            if(step % try_step==0) {
                finish = std::chrono::system_clock::now();
                auto elapsed_ms  = std::chrono::duration_cast<duration_ms>(finish - start);
                if(elapsed_ms >= cutoff_time) {
                    if(verbose) {
                        std::cout << "Time limit reached:" << elapsed_ms.count() << "ms" << std::endl;
                    }
                    return;
                }
            }

            /* choose a vertex u in C with the highest dscore,
              breaking ties in favor of the oldest one; */
            best_cov_v = remove_cand[0];
            #ifdef _OPENMP
            std::vector<int> best_cov_cands;
            #pragma omp parallel firstprivate(best_cov_v)
            #endif
            {
              #ifdef _OPENMP
              auto num_threads = omp_get_num_threads();
              auto thread_id   = omp_get_thread_num();
              #pragma omp single
              {
                best_cov_cands.resize(num_threads);
              }
              #pragma omp for schedule(static,1024)
              #endif
              for (i=1; i<remove_cand_size; ++i) {
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
              #ifdef _OPENMP
              best_cov_cands[thread_id] = best_cov_v;
              #endif
            }
            #ifdef _OPENMP
            best_cov_v = best_cov_cands[0];
            for(const auto & vert : best_cov_cands){
                if( dscore[vert] < dscore[best_cov_v])
                    continue;
                else if( dscore[vert] > dscore[best_cov_v] )
                    best_cov_v = vert;
                else if (time_stamp[vert] < time_stamp[best_cov_v])
                    best_cov_v = vert;
            }
            #endif

            /* C := C\{u}, confChange(u) := 0 and confChange(z) := 1 for each z in N(u); */
            remove(best_cov_v);


            /* choose an uncovered edge e randomly; */
            e = uncov_stack[mt_rand()%uncov_stack_fill_pointer];

            /* choose a vertex v in e such that confChange(v) = 1 with higher dscore,
              breaking ties in favor of the older one; */
            v1 = edge[e].v1;
            v2 = edge[e].v2;

            if(conf_change[v1]==0) best_add_v=v2;
            else if(conf_change[v2]==0) best_add_v=v1;
            else {
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

            ++step;
        }
        return;
    }

    /**
     * Start solver with this initial cover, given as a list of vertex indices
     */
    void set_initial_cover(const std::vector<int> & cover){
        c_size = cover.size();
        std::fill(v_in_c.begin(), v_in_c.end(), false);
        for(const auto v : cover){
            v_in_c[v] = true;
        }
        reset_remove_cand();
        update_best_cov_v();
    }

    /**
     * Start solver with this initial cover, given as a list of flags which
     * denote if the vertex is in the cover
     */
    void set_initial_cover(const std::vector<bool> & cover){
        v_in_c = cover;
        c_size = std::count(cover.begin(), cover.end(), true);
    }

    /*On solution*/
    void print_solution()
    {
        int mis_vertex_count=0;

        for (int i=1; i<=v_num; i++) {
            if (!best_v_in_c[i])
                mis_vertex_count++;
        }

        if(mis_vertex_count+best_c_size!=v_num) {
            std::cout << "The size of independent set + the size of "
                      << "vertex cover is not equal to |V(G)|!" << std::endl;
        }

        std::cout << "c Best found independent set size = "
                  << mis_vertex_count << std::endl;
        std::cout << "c The following output is the found independent set."
                  << std::endl;

        for (int i=1; i<=v_num; i++) {
            if (!best_v_in_c[i])//output max independent set
                std::cout << i << "  ";
        }
        std::cout << std::endl;
    }

    //check whether the solution found is a proper solution
    bool check_solution()
    {
        int e;
        for(e=0; e<e_num; ++e) {
            if(!(best_v_in_c[edge[e].v1] || best_v_in_c[edge[e].v2])) {
                std::cout << "uncovered edge " << e << std::endl;
                return false;
            }
        }

        return true;
    }

    /**
     * return vertex indices of current best vertex cover
     */
    std::vector<int> get_cover()
    {
        std::vector<int> cover;
        for (int i=1; i<=v_num; i++) {
            if (best_v_in_c[i]) {
                cover.push_back(i);
            }
        }
        return cover;
    }
    std::vector<bool> get_cover_as_flaglist()
    {
        return v_in_c;
    }

    /**
     * return vertex indices of current best independent set
     */
    std::vector<int> get_independent_set()
    {
        std::vector<int> iset;
        for (int i=1; i<=v_num; i++) {
            if (!best_v_in_c[i]) {
                iset.push_back(i);
            }
        }
        return iset;
    }

    /**
     * Number of vertices
     */
    inline int get_vertex_count() const
    {
        return v_num;
    }

    /**
     * Number of edges
     */
    inline int get_edge_count() const
    {
        return e_num;
    }

    /**
     * Size of current best vertex cover
     */
    inline int get_best_cover_size() const
    {
        return best_c_size;
    }

    /**
     * Tries necessary for current best vertex cover
     */
    inline long get_best_step() const
    {
        return best_step;
    }

    /**
     * duration for calculating current best vertex cover
     */
    inline std::chrono::milliseconds get_best_duration() const
    {
        return std::chrono::duration_cast<
               std::chrono::milliseconds>(best_comp_time);
    }

    /**
     * total duration since start of calculation
     */
    inline std::chrono::milliseconds get_total_duration() const
    {
        return std::chrono::duration_cast<
               std::chrono::milliseconds>(
                   (std::chrono::system_clock::now() - start));
    }

    /**
     * Print statistics during calculation
     */
    static bool default_stats_printer(const NuMVC& solver)
    {
        auto time_ms = std::chrono::duration_cast<
                       std::chrono::milliseconds>(solver.get_best_duration());
        std::cout << "Better MVC found.\tSize: "
                  << solver.get_best_cover_size()
                  << "\tTime: " << std::fixed << std::setw(4) << std::setprecision(4)
                  << time_ms.count() << "ms" << std::endl;
        return false;
    }

};

#endif

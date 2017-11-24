/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/


/************************************************
** Date:  2015.2.2
** FastVC
** Author: Shaowei Cai, caisw@ios.ac.cn
**       Key Laboratory of Computer Science,
**       Institute of Software, Chinese Academy of Sciences,
**       Beijing, China
**
** Date: 2017.11.17
** Modify: Felix Moessbauer
** Object oriented design, C++11 support
*************************************************/

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <numeric>
#include <functional>

/**
 * Local search solver for Minimum Vertex Cover
 */
class FastVC {
private:
    using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
    using duration_ms = std::chrono::milliseconds;

    /* print messages while solving */
    bool verbose = false;

    timepoint_t start, finish;

    struct Edge {
        int v1;
        int v2;
    };

    /*parameters of algorithm*/
    long long   max_steps;      //step limit
    duration_ms cutoff_time;    //time limit
    long long   step;
    int         optimal_size;   //terminate the algorithm before step limit if it finds a vertex cover of optimal_size

    /*parameters of the instance*/
    int    v_num; //|V|: 1...v
    int    e_num; //|E|: 0...e-1

    /*structures about edge*/
    std::vector<Edge> edge;

    /*structures about vertex*/
    std::vector<int>       dscore;     //dscore of v
    std::vector<long long> time_stamp;


    //from vertex to it's edges and neighbors
    std::vector<int> v_beg_idx; // v_edges and v_adj is flattened 2-d array, hence store indices
    std::vector<int> v_edges;   //edges related to v, v_edges[i][k] means vertex v_i's k_th edge
    std::vector<int> v_adj;     //v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
    std::vector<int> v_degree;  //amount of edges (neighbors) related to v


    /* structures about solution */
    //current candidate solution
    int               c_size;              //cardinality of C
    std::vector<bool> v_in_c;              //a flag indicates whether a vertex is in C
    //remove candidates, an array consists of only vertices in C, not including tabu_remove
    //use custom stack for perfomance reasons
    std::vector<int>  remove_cand;
    std::vector<int>  index_in_remove_cand;
    int               remove_cand_size;

    //best solution found
    int               best_c_size;
    std::vector<bool> best_v_in_c;        //a flag indicates whether a vertex is in best solution
    duration_ms       best_comp_time;
    long              best_step;


    //uncovered edge stack
    std::vector<int>  uncov_stack;          //store the uncov edge number
    std::vector<int>  index_in_uncov_stack; //which position is an edge in the uncov_stack

    int                  edge_cand;
    static constexpr int try_step = 10;     // check time each number of steps

public:
    /**
     * Construct solver instance by importing a graph in DIMACS format
     */
    template<
        typename Is,
        typename Duration>
    FastVC(
        /// Input Stream with graph in DIMACS format
        Is & str,
        /// Size of optimal vertex cover (set to 0 if not known)
        int optimal_size,
        /// Stop calculation after this duration (chrono duration)
        Duration cutoff_time,
        /// Print messages during calculation
        bool verbose = false)
        : verbose(verbose),
          cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
          optimal_size(optimal_size)
    {
        build_instance(str);
        init_sol();
    }

private:

    void update_best_sol()
    {
        for (int i=1; i<=v_num; ++i)
            best_v_in_c[i] = v_in_c[i];

        best_c_size = c_size;
        finish = std::chrono::system_clock::now();
        best_comp_time = std::chrono::duration_cast<duration_ms>(finish - start);
        best_step = step;
    }

    template<typename Is>
    bool build_instance(Is & str)
    {
        char line[1024];
        char tempstr1[10];
        char tempstr2[10];
        int  e;

        char tmp;
        int  v1,v2;

        /*** build problem data structures of the instance ***/
        str.getline(line,1024);
        while (line[0] != 'p') str.getline(line,1024);
        sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);

        edge.resize(e_num);              //be initialized here
        index_in_uncov_stack.resize(e_num);             //the same as above
        dscore.resize(v_num + 1);                       //be initialized in init_sol()
        time_stamp.resize(v_num + 1);             //be initialized in init_sol()
        v_degree.resize(v_num + 1);                     //the same as above
        v_in_c.resize(v_num + 1);                      //be initialized in init_sol()
        remove_cand.resize(v_num + 1);                  //be initialized in reset_remove_cand() in init_sol()
        index_in_remove_cand.resize(v_num + 1);         //the same as above
        best_v_in_c.resize(v_num + 1);                 //be initialized in update_best_sol() in init_sol()
        /* read edges and compute v_degree */

        for (e=0; e<e_num; ++e) {
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

        //  int v_degree_tmp[MAXV];
        std::vector<int> v_degree_tmp(v_num + 1);
        for (e=0; e<e_num; ++e) {
            v1=edge[e].v1;
            v2=edge[e].v2;

            v_edges[v_beg_idx.at(v1) + v_degree_tmp[v1]] = e;
            v_edges[v_beg_idx.at(v2) + v_degree_tmp[v2]] = e;

            v_adj[v_beg_idx[v1] + v_degree_tmp[v1]] = v2;
            v_adj[v_beg_idx[v2] + v_degree_tmp[v2]] = v1;

            v_degree_tmp[v1]++;
            v_degree_tmp[v2]++;
        }
        return true;
    }

    void reset_remove_cand()
    {
        int v,j;
        j=0;
        for (v=1; v<=v_num; ++v) {
            if(v_in_c[v]==1) {
                remove_cand[j] = v;
                index_in_remove_cand[v]=j;
                ++j;
            } else index_in_remove_cand[v]=0;
        }

        remove_cand_size = j;
    }

    void update_target_size()
    {
        --c_size;

        int v,i;
        //int best_dscore;
        int best_remove_v;//vertex with the highest improvement in C

        best_remove_v = remove_cand[0];
        //best_dscore = dscore[best_remove_v];

        if(dscore[best_remove_v]!=0) {
            for (i=1; i<remove_cand_size; ++i) {
                v = remove_cand[i];

                if(dscore[v]==0) break;

                if (dscore[v] > dscore[best_remove_v])
                    best_remove_v = v;
            }
        }

        remove(best_remove_v);

        //remove best_remove_v from remove_cand, and move the last vertex in remove_cand to the position
        int last_remove_cand_v = remove_cand[--remove_cand_size];
        int index = index_in_remove_cand[best_remove_v];
        remove_cand[index] = last_remove_cand_v;
        index_in_remove_cand[last_remove_cand_v] = index;

        //reset_remove_cand();
    }




    //update the best vertex in C
    int choose_remove_v()
    {
        int cand_count=50;
        int i,v;
        int best_v = remove_cand[rand()%remove_cand_size];

        for (i=1; i<cand_count; ++i) {
            v = remove_cand[rand()%remove_cand_size];

            if( dscore[v] < dscore[best_v])
                continue;
            else if( dscore[v]> dscore[best_v] )
                best_v = v;
            else if (time_stamp[v]<time_stamp[best_v])
                best_v = v;
        }

        return best_v;
    }

    inline void uncover(int e)
    {
        index_in_uncov_stack[e] = uncov_stack.size();
        uncov_stack.push_back(e);
    }

    inline void cover(int e)
    {
        int index,last_uncov_edge;

        //since the edge is satisfied, its position can be reused to store the last_uncov_edge
        last_uncov_edge = uncov_stack.back();

        index = index_in_uncov_stack[e];
        uncov_stack[index] = last_uncov_edge;
        index_in_uncov_stack[last_uncov_edge] = index;
        uncov_stack.pop_back();
    }

    void init_sol()
    {
        int v,e;
        int v1, v2;

        /*** build solution data structures of the instance ***/

        c_size = 0;
        for (e=0; e<e_num; e++) {
            v1=edge[e].v1;
            v2=edge[e].v2;

            if (v_in_c[v1]==0 && v_in_c[v2]==0) { //if uncovered, choose the endpoint with higher degree
                if(v_degree[v1] > v_degree[v2]) {
                    v_in_c[v1]=1;
                } else {
                    v_in_c[v2]=1;
                }
                c_size++;
            }
        }

        //calculate dscores
        for (e=0; e<e_num; e++) {
            v1=edge[e].v1;
            v2=edge[e].v2;

            if (v_in_c[v1]==1 && v_in_c[v2]==0) dscore[v1]--;
            else if (v_in_c[v2]==1 && v_in_c[v1]==0) dscore[v2]--;
        }


        //remove redundent vertices

        for (v=1; v<=v_num; v++) {
            if (v_in_c[v]==1 && dscore[v]==0) {
                remove(v);
                c_size--;
            }
        }


        update_best_sol();//initialize the best found solution

        reset_remove_cand();

    }

    void add(int v)
    {
        v_in_c[v] = 1;
        dscore[v] = -dscore[v];

        int i,e,n;

        int edge_count = v_degree[v];

        for (i=0; i<edge_count; ++i) {
            e = v_edges[v_beg_idx[v]+i];// v's i'th edge
            n = v_adj[v_beg_idx[v]+i];//v's i'th neighbor

            if (v_in_c[n]==0) { //this adj isn't in cover set
                dscore[n]--;
                //conf_change[n] = 1;

                cover(e);
            } else {
                dscore[n]++;
            }
        }
    }

    void remove(int v)
    {
        v_in_c[v] = 0;
        dscore[v] = -dscore[v];
        //conf_change[v] = 0;

        int i,e,n;

        int edge_count = v_degree[v];
        for (i=0; i<edge_count; ++i) {
            e = v_edges[v_beg_idx[v] + i];
            n = v_adj[v_beg_idx[v] + i];

            if (v_in_c[n]==0) { //this adj isn't in cover set
                dscore[n]++;
                //conf_change[n] = 1;

                uncover(e);
            } else {
                dscore[n]--;
            }
        }
    }

public:

    /**
     * print node indices of current best vertex cover
     */
    void print_solution() const
    {
        for (int i=1; i<=v_num; i++) {
            if (best_v_in_c[i]==1)//output vertex cover
                std::cout<<i<<'\t';
        }
        std::cout<<std::endl;
    }

    /**
     * Check if the solution is valid
     */
    bool check_solution()
    {
        for(int e=0; e<e_num; ++e) {
            if(best_v_in_c[edge[e].v1]!=1 && best_v_in_c[edge[e].v2]!=1) {
                if(verbose) {
                    std::cout<<"c error: uncovered edge "<<e<<std::endl;
                }
                return false;
            }
        }

        int verified_vc_size=0;
        for (int i=1; i<=v_num; i++) {
            if (best_v_in_c[i]==1)
                verified_vc_size++;
        }

        if(best_c_size==verified_vc_size) return true;

        else {
            if(verbose) {
                std::cout<<"c error: claimed best found vc size!=verified vc size"<<std::endl;
                std::cout<<"c claimed best found vc size="<<best_c_size<<std::endl;
                std::cout<<"c verified vc size="<<verified_vc_size<<std::endl;
            }
            return false;
        }
    }

    /**
     * calculate minimum vertex cover
     */
    void cover_LS()
    {
        cover_LS([]() {
            return true;
        });
    }

    /**
     * calculate minimum vertex cover and call callback after
     * each iteration. If callback returns true, stop calculation.
     */
    void cover_LS(const std::function<bool (void)> & callback_on_update)
    {
        int    remove_v, add_v;
        int    e,v1,v2;

        step = 1;
        start = std::chrono::system_clock::now();

        while(1) {
            if (uncov_stack.size() == 0) { //update best solution if needed
                update_best_sol();
                if(callback_on_update()) {
                    return;
                }
                if (c_size==optimal_size) return;

                update_target_size();

                continue;
            }

            if(step%try_step==0) { //check cutoff
                finish = std::chrono::system_clock::now();
                auto elap_time = finish - start;
                if(elap_time >= cutoff_time) return;
            }

            remove_v = choose_remove_v();
            //remove_dscore = dscore[remove_v];

            remove(remove_v);

            e = uncov_stack[rand()%uncov_stack.size()];
            v1 = edge[e].v1;
            v2 = edge[e].v2;

            if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
                add_v=v1;
            else add_v=v2;

            add(add_v);

            int index = index_in_remove_cand[remove_v];
            index_in_remove_cand[remove_v] = 0;

            remove_cand[index] = add_v;
            index_in_remove_cand[add_v] = index;

            time_stamp[add_v]=time_stamp[remove_v]=step;

            //tabu_remove = add_v;

            //update_edge_weight();

            ++step;
        }
    }

    /**
     * return vertex indices of current best vertex cover
     */
    std::vector<int> get_cover()
    {
        std::vector<int> cover;
        for (int i=1; i<=v_num; i++) {
            if (best_v_in_c[i]==1) {
                cover.push_back(i);
            }
        }
        return cover;
    }

    /**
     * return vertex indices of current best independent set
     */
    std::vector<int> get_independent_set()
    {
        std::vector<int> iset;
        for (int i=1; i<=v_num; i++) {
            if (best_v_in_c[i]==0) {
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
static bool default_stats_printer(const FastVC & solver){
    auto time_ms = std::chrono::duration_cast<
                     std::chrono::milliseconds>(solver.get_best_duration());
    std::cout << "Better MVC found.\tSize: "
              << solver.get_best_cover_size()
              << "\tTime: " << std::fixed << std::setw(4) << std::setprecision(4)
              << time_ms.count() << "ms" << std::endl;
    return false;
}

};


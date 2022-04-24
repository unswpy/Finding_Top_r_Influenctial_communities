//
// Created by sbian on 2020/7/20.
//

#ifndef WEIGHT_COMMUNITY_TOPR_H
#define WEIGHT_COMMUNITY_TOPR_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <cmath>

#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <map>
#include <hash_map>
#include <set>
#include <stack>
#include <queue>

using namespace std;

class TopRComm{
private:
    int degree_c; // the coreness of the community
    int size_c; // top-r community returned
    int graph_size; // the original graph size
    int edge_size; // edge size
    int c_size; // the current graph size after some nodes removed
    int max_core; // max core of this graph
    double min_thres; // minimum threshold the weight of node must satisfy

    double inf_value;

    std::vector<int> core;
    std::vector<double> weight;
    std::vector< std::vector<int> > AdjList;
    std::vector< std::vector<int> > AdjList_copy;
    // control the storage
//    std::vector<map<int, int>> AdjList_Loc;
    std::vector< __gnu_cxx::hash_map<int, int> > AdjList_Loc;

    const double err = 0.0000000000000000001;
public:
    TopRComm();
    ~TopRComm();
    void init(int numV);
    void load_data(const string & adjpath, const string & corepath, const string & pagepath);

    /// Ronghua's solution
    void naive_max_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void cons_max_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// Ronghua's solution extension
    void naive_min_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void cons_min_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    void non_min_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// sum algorithm
    void naive_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void improved_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void approx_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double eps);
//    void exact_avg_global_topr(const string & adjpath, const string & corepath,
//            const string & pagepath, int k, int r);
//    void improved_avg_global_topr(const string & adjpath, const string & corepath,
//            const string & pagepath, int k, int r);
    /// avg algorithm
    /// delete the smallest weight node, then compute the k core, the order of nodes is fixed
    void improved_greedy_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    /// randomly delete the node, then compute the k core, the order of nodes is fixed
    void improved_random_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    /// each iteration delete the smallest weight node, then compute the k core
    void improved_climb_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    /// each iteration delete node which could improve the influence value most, then compute the k core
    void improved_double_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    /// each iteration delete node which could improve the influence value most, then compute the k core (cut vertex)
    void cut_double_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);

    /// constraint sum algorithm
    void naive_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void improved_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r);
    void approx_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double eps);
    void non_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);

    /// non-overlapping size-constraint sum algorithm
    void climb_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);
    void degree_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);
    void casual_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);

    /// constraint avg algorithm
    /// delete the smallest weight node, then compute the k core, the order of nodes is fixed
    void improved_cons_greedy_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// randomly delete the node, then compute the k core, the order of nodes is fixed
    void improved_cons_random_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// each iteration delete the smallest weight node, then compute the k core
    void improved_cons_climb_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// each iteration delete node which could improve the influence value most, then compute the k core
    void improved_cons_double_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// each iteration delete node which could improve the influence value most, then compute the k core (cut vertex)
    /// it still has bugs?
    void cut_cons_double_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);
    /// exact algorithm
    void improved_cons_exact_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, double p);

    /// non-overlapping avg algorithm
    /// delete the smallest weight node at each iteration, then compute the k core
    void climb_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);
    ///  local search weightest node
    void degree_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);
    /// random delete the node at each iteration, then compute the k core
    void casual_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
            const string & pagepath, int k, int r, int s);

    void split(const string& src, const string& separator, vector<string>& dest);
    void add_edge(int start, int end);
    void remove_edge(int start, int end);
    void compute_max_k_core(int k);
    void compute_cons_max_k_core(int k, double p);
    vector<int> compute_cut_vertex();
    pair<vector<int>, int> compute_minimal_connected_components();
    pair<vector<int>, int> compute_maximal_connected_components();
    vector<vector<int>> compute_connected_components();
    vector<vector<int>> compute_connected_k_core(int k);
    vector<vector<int>> compute_part_connected_components(const vector<int> & vertex_set); // partial connected components
    void dfs_deletion(int v, int k); // after remove the max/min vertex, calculate the new graph
    bool is_connected_k_core(const vector<int> & vertex_set, int k);
    bool is_connected_k_core_graph(int k);
    bool is_element_in_vector(const vector<int> & vec, int element);
    bool is_keep_avg_component(const vector<int> & vec, double threshold, int k);
    int greedy_unsatisfied_node(int k);
    int random_unsatisfied_node(const vector<int> & vec, int k);
    double sum_vector(const vector<int> & vec);
    double avg_vector(const vector<int> & vec);
    vector<pair<int, int>> batch_del_edges(const vector<int> & nodes);
    vector<pair<int, int>> single_del_edges(int node);
    void batch_add_edges(const vector<pair<int, int>> & edges);

    void combine_inner(const vector<int> & data, int start, int n, int m, int depth,
            vector<int> temp, vector<vector<int>> result);
    vector<vector<int>> combine(const vector<int> & data, int m);
    bool determined_k_core(const vector<int> & node_vec, const bool *is_exist, int k);
    int min_degree(const vector<int> & node_vec);
    vector<vector<int>> getSubsets(const vector<int> & A, int k);
    vector<int> peel_vector(vector<int> vec, int k);
};

#endif //WEIGHT_COMMUNITY_TOPR_H

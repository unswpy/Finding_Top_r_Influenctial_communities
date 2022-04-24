//
// Created by sbian on 2020/7/22.
//

#ifndef WEIGHT_COMMUNITY_ENUM_H
#define WEIGHT_COMMUNITY_ENUM_H

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
#include <queue>

using namespace std;

class EnumComm{
private:
    int degree_c; // degree constraint
    double thres_c; // threshold constraint
    int graph_size; // the original graph size
    int edge_size; // edge size
    int c_size; // the current graph size after some nodes removed
    int max_core; // max core of this graph

    double inf_value;

    std::vector<int> core;
    std::vector<double> weight;
    std::vector<std::vector<int>> AdjList;
    std::vector< __gnu_cxx::hash_map<int, int> > AdjList_Loc;

    const double err = 0.0000000001;
public:
    EnumComm();
    ~EnumComm();
    void init(int numV);
    void load_data(const string & adjpath, const string & corepath,
            const string & pagepath);

    void naive_max_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau); // Ronghua's solution extension
    void naive_min_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau); // Ronghua's solution extension
    void naive_sum_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void improved_sum_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void naive_avg_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void improved_avg_global_enum(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);

    void split(const string& src, const string& separator, vector<string>& dest);
    void add_edge(int start, int end);
    void remove_edge(int start, int end);
    void compute_max_k_core(int k);
    vector<vector<int>> compute_connected_components();
    vector<vector<int>> compute_connected_k_core(int k);
    bool is_connected_k_core_graph(int k);
    double sum_vector(const vector<int> & vec);
    vector<pair<int, int>> batch_del_edges(const vector<int> & nodes);
    vector<pair<int, int>> single_del_edges(const int node);
    void batch_add_edges(const vector<pair<int, int>> & edges);
};

#endif //WEIGHT_COMMUNITY_ENUM_H

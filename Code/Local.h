//
// Created by sbian on 2020/7/25.
//

#ifndef WEIGHT_COMMUNITY_LOCAL_H
#define WEIGHT_COMMUNITY_LOCAL_H

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

class LocalComm{
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
    LocalComm();
    ~LocalComm();
    void init(int numV);
    void load_data(const string & adjpath, const string & corepath, const string & pagepath);

    void naive_max_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void naive_min_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void naive_sum_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void improved_sum_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void self_avg_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);
    void neighbor_avg_local(const string & adjpath, const string & corepath,
            const string & pagepath, int k, double tau);

    void split(const string& src, const string& separator, vector<string>& dest);
    void add_edge(int start, const int end);
    void remove_edge(int start, int end);
    vector<vector<int>> preprocess_local_avg(int k, double tau);
    vector<vector<int>> self_greedy(vector<int> & vec, int k, double tau);
    vector<vector<int>> neighbor_greedy(vector<int> & vec, int k, double tau);
};



#endif //WEIGHT_COMMUNITY_LOCAL_H

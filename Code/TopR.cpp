//
// Created by sbian on 2020/7/20.
//

#include "TopR.h"
#include "serialize.h"

TopRComm::TopRComm() {
    degree_c = -1; // the coreness of the community
    size_c = -1; // top-r community returned
    graph_size = -1; // the original graph size
    edge_size = -1;
    c_size = -1; // the current graph size after some nodes removed
    max_core = 0;
    inf_value = 0.0;
}

TopRComm::~TopRComm() {
    core.clear();
    weight.clear();
    for(auto & iter : AdjList) {
        iter.clear();
    }
    AdjList.clear();
    for (auto & iter : AdjList_Loc) {
        iter.clear();
    }
    AdjList_Loc.clear();
}

void TopRComm::init(const int numV) {
    double local;
    graph_size = numV;
    c_size = numV;
    AdjList_Loc.resize(numV);
    for (int i = 0; i < graph_size; i++) {
        for (unsigned int j = 0; j < AdjList[i].size(); j++) {
            AdjList_Loc[i][AdjList[i][j]] = int(j);
        }
    }
}

void TopRComm::load_data(const string & adjpath, const string & corepath, const string & pagepath) {
    load_serialized_graph(adjpath, AdjList);
    load_serialized_graph(corepath, core);
    load_serialized_graph(pagepath, weight);
    graph_size = AdjList[AdjList.size() - 1][0];
    edge_size = AdjList[AdjList.size() - 1][1];
    cout << AdjList[AdjList.size() - 1][0] << " " << AdjList[AdjList.size() - 1][1] << " "
         << AdjList[AdjList.size() - 1][2] << " " << AdjList[AdjList.size() - 1][3] << " "
         << AdjList[AdjList.size() - 1][4] << endl;
    AdjList.pop_back();
    AdjList_copy = AdjList;
    cout << "AdjList load" << endl;
    init(graph_size);
    assert(graph_size == AdjList.size());

    for(int i = 0; i < graph_size; i++) {
        if (core[i] > max_core) {
            max_core = core[i];
        }
    }
    cout << "max_core: " << max_core << endl;
//    /// check load graph
//    int deg_sum = 0;
//    for(int i = 0; i < graph_size; i++) {
//        deg_sum += deg[i];
//    }
//    cout << deg_sum << " " << 2 * edge_size << endl;
//    assert(deg_sum == 2 * edge_size);
}

void TopRComm::naive_max_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    int index = 0;
    int del_id = 0;
    pair<vector<int>, int> res;
    vector<pair<vector<int>, double>> output_res;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    compute_max_k_core(degree_c);
    cout << "c_size: " << c_size << endl;
    while (c_size != 0) {
        cout << "start" << endl;
        res = compute_maximal_connected_components();
        output_res.emplace_back(res.first, weight[res.second]);
        del_id = res.second;
        cout << AdjList[del_id].size() << endl;
        dfs_deletion(del_id, k);
        index = index + 1;
        cout << "index: " << index << endl;
        if (index >= r) {
            break;
        }
        cout << "end" << endl;
    }
//    for (int v = 0; v < graph_size; v++) {
//        assert(deg[v] == AdjList[v].size());
//    }
    for (int i = 0; i < output_res.size(); i++) {
        cout << i << " " << output_res[i].second << endl;
    }
//    res = compute_connected_components();
//    for (int i = 0; i < res.size(); i++) {
//        for (int j = 0; j < res[i].size(); j++) {
//            cout << res[i][j] << " ";
//        }
//        cout << endl;
//    }

//    while (c_size != 0) {
//    }
}

void TopRComm::cons_max_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    int max_idx;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        for (unsigned int j = 0; j < res[i].size(); j++) {
            k_core_vec.push_back(res[i][j]);
        }
    }

    if (res.size() == 0) {
        return;
    }

    L.clear();
    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            int max_node = res[i][0];
            for (unsigned int j = 0; j < res[i].size(); j++) {
                if (weight[max_node] < weight[res[i][j]]) {
                    max_node = res[i][j];
                }
            }
            L.push_back(make_pair(res[i], max_node));
        }
    }

    cout << "cons max top-r begin" << endl;

    batch_del_edges_vec.clear();
    while (output_res.size() < r) {
        add_vec.clear();
        max_idx = 0;
        delete_node = L[0].second;
        for (unsigned int i = 0; i < L.size(); i++) {
            if (weight[delete_node] < weight[L[i].second]) {
                delete_node = L[i].second;
                max_idx = i;
            }
        }
        output_res.push_back(L[max_idx]);
        single_del_edges_vec.clear();
        single_del_edges_vec = single_del_edges(delete_node);
        L.clear();
        vector<vector<int>> divide_com = compute_connected_k_core(k);
        for (unsigned int i = 0; i < divide_com.size(); i++) {
            if (divide_com[i].size() <= k) {
            }
            else {
                int max_node = divide_com[i][0];
                for (unsigned int j = 0; j < divide_com[i].size(); j++) {
                    add_vec.push_back(divide_com[i][j]);
                    if (weight[max_node] < weight[divide_com[i][j]]) {
                        max_node = divide_com[i][j];
                    }
                }
                L.push_back(make_pair(divide_com[i], max_node));
            }
        }
        batch_add_edges(single_del_edges_vec);
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        single_del_edges_vec.clear();
        single_del_edges_vec = batch_del_edges(del_nodes_vec);
        batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                single_del_edges_vec.begin(), single_del_edges_vec.end());
        k_core_vec = add_vec;
    }
    batch_add_edges(batch_del_edges_vec);

    cout << "cons max top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/cons_max_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
    }
    for (unsigned int i = 0; i < output_res[r-1].first.size(); i++) {
        outfile << output_res[r-1].first[i] << " ";
    }
    outfile << endl;
    outfile.close();
}

void TopRComm::naive_min_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    int index = 0;
    int del_id = 0;
    pair<vector<int>, int> res;
    vector<pair<vector<int>, double>> output_res;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    compute_max_k_core(degree_c);
    cout << "c_size: " << c_size << endl;
    while (c_size != 0) {
        cout << "start" << endl;
        res = compute_minimal_connected_components();
        output_res.emplace_back(res.first, weight[res.second]);
        del_id = res.second;
        cout << AdjList[del_id].size() << endl;
        dfs_deletion(del_id, k);
        index = index + 1;
        cout << "index: " << index << endl;
        cout << "end" << endl;
    }
    if (index > r) {
        for(int i = index - 1; i >= index - r; i--) {
            cout << i << " " << output_res[i].second << endl;
        }
    }
    else {
        for(int i = index - 1; i >= 0; i--) {
            cout << i << " " << output_res[i].second << endl;
        }
    }
}

void TopRComm::cons_min_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    int min_idx;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        for (unsigned int j = 0; j < res[i].size(); j++) {
            k_core_vec.push_back(res[i][j]);
        }
    }

    if (res.size() == 0) {
        return;
    }

    L.clear();
    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            int min_node = res[i][0];
            for (unsigned int j = 0; j < res[i].size(); j++) {
                if (weight[min_node] < weight[res[i][j]]) {
                    min_node = res[i][j];
                }
            }
            L.push_back(make_pair(res[i], min_node));
        }
    }

    cout << "cons min top-r begin" << endl;

    batch_del_edges_vec.clear();
    while(k_core_vec.size() != 0) {
        add_vec.clear();
        min_idx = 0;
        delete_node = L[0].second;
        for (unsigned int i = 0; i < L.size(); i++) {
            if (weight[delete_node] > weight[L[i].second]) {
                delete_node = L[i].second;
                min_idx = i;
            }
        }
        output_res.push_back(make_pair(L[min_idx].first, weight[L[min_idx].second]));
        single_del_edges_vec.clear();
        single_del_edges_vec = single_del_edges(delete_node);
        L.clear();
        vector<vector<int>> divide_com = compute_connected_k_core(k);
        for (unsigned int i = 0; i < divide_com.size(); i++) {
            if (divide_com[i].size() <= k) {
            }
            else {
                int min_node = divide_com[i][0];
                for (unsigned int j = 0; j < divide_com[i].size(); j++) {
                    add_vec.push_back(divide_com[i][j]);
                    if (weight[min_node] > weight[divide_com[i][j]]) {
                        min_node = divide_com[i][j];
                    }
                }
                L.push_back(make_pair(divide_com[i], min_node));
            }
        }
        batch_add_edges(single_del_edges_vec);
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        single_del_edges_vec.clear();
        single_del_edges_vec = batch_del_edges(del_nodes_vec);
        batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                single_del_edges_vec.begin(), single_del_edges_vec.end());
        k_core_vec = add_vec;
    }
    batch_add_edges(batch_del_edges_vec);

    cout << "cons min top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/cons_min_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() > r) {
        for (unsigned int i = output_res.size() - 1; i >= output_res.size() - r; i--) {
            outfile << output_res.size() - i << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = output_res.size() - 1; i >= 0; i--) {
            outfile << output_res.size() - i << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
//    if (output_res.size() >= r) {
//        for (unsigned int i = 0; i < output_res[output_res.size() - r].first.size(); i++) {
//            outfile << output_res[output_res.size() - r].first[i] << " ";
//        }
//        outfile << endl;
//    }
//    else {
//        for (unsigned int i = 0; i < output_res[0].first.size(); i++) {
//            outfile << output_res[0].first[i] << " ";
//        }
//        outfile << endl;
//    }
    outfile.close();
}

void TopRComm::non_min_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    queue<pair<vector<int>, double>> Q;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    int min_idx;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        for (unsigned int j = 0; j < res[i].size(); j++) {
            k_core_vec.push_back(res[i][j]);
        }
    }

    if (res.size() == 0) {
        return;
    }

    while (!Q.empty()) {
        Q.pop();
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            int min_node = res[i][0];
            for (unsigned int j = 0; j < res[i].size(); j++) {
                if (weight[min_node] < weight[res[i][j]]) {
                    min_node = res[i][j];
                }
            }
            Q.push(make_pair(res[i], min_node));
        }
    }

    cout << "Q size: " << Q.size() << endl;

    cout << "non min top-r begin" << endl;

    while(Q.size() != 0) {
        add_vec.clear();
        add_vec = Q.front().first;
        delete_node = Q.front().second;
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        batch_del_edges_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        single_del_edges_vec.clear();
        single_del_edges_vec = single_del_edges(delete_node);
        vector<vector<int>> divide_com = compute_connected_k_core(k);
        if (divide_com.size() == 0) {
            output_res.push_back(make_pair(add_vec, weight[delete_node]));
        }
        else {
            for (unsigned int i = 0; i < divide_com.size(); i++) {
                if (divide_com[i].size() <= k) {
                }
                else {
                    int min_node = divide_com[i][0];
                    for (unsigned int j = 0; j < divide_com[i].size(); j++) {
                        if (weight[min_node] > weight[divide_com[i][j]]) {
                            min_node = divide_com[i][j];
                        }
                    }
                    Q.push(make_pair(divide_com[i], min_node));
                }
            }
        }
        Q.pop();
        batch_add_edges(single_del_edges_vec);
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "non min top-r end" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int> ,double> & right){
             return left.second > right.second;
         });

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/non_min_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() > r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                    " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                    " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::naive_sum_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> output_res;
    vector<pair<vector<int>, double>> Lc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    double sum_vec;
    double top_r_weight;
    bool flag = false; // determine whether duplicate
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    bool *is_exist = new bool[graph_size];
//    bool is_exist[graph_size];
    for (int i = 0; i < graph_size; i++) {
        is_exist[i] = false;
    }
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            output_res.push_back(make_pair(res[i], sum_vec));
        }
    }

    sort(output_res.begin(), output_res.end(),
            [](const pair<vector<int>, double> & left, const pair<vector<int> ,double> & right){
                return left.second > right.second;
    });

    if (output_res.size() > r) {
        output_res.assign(output_res.begin(), output_res.begin() + r);
    }
    top_r_weight = output_res.back().second;

    cout << "init: " << output_res.size() << endl;

    for (unsigned int i = 0; i < output_res.size(); i++) {
        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
            is_exist[output_res[i].first[j]] = true;
            k_core_vec.push_back(output_res[i].first[j]);
        }
    }

    cout << "naive sum top-r begin" << endl;

    for (int i = 0; i < graph_size; i++) {
        if (is_exist[i] == true) {
            Lc.clear();
            for (int j = 0; j < output_res.size(); j++) {
                if (output_res[j].first.size() <= k) {
                }
                else if (output_res[j].second - weight[i] <= top_r_weight && output_res.size() + Lc.size() >= r) {
                }
                else {
                    add_vec.assign(output_res[j].first.begin(), output_res[j].first.end());
                    for (unsigned int m = 0; m < add_vec.size(); m++) {
                        if (add_vec[m] == i) {
                            add_vec.erase(add_vec.begin() + m);
                            break;
                        }
                    }
                    if (add_vec.size() == output_res[j].first.size()) {
                    }
                    else {
                        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
                        set<int> add_set(add_vec.begin(), add_vec.end());
                        set<int> del_nodes_set;
                        del_nodes_vec.clear();
                        set_difference(k_core_set.begin(), k_core_set.end(),
                                       add_set.begin(), add_set.end(),
                                       inserter(del_nodes_set, del_nodes_set.begin()));
                        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
                        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
                        vector<vector<int>> divide_com = compute_connected_k_core(k);
                        for (unsigned int m = 0; m < divide_com.size(); m++) {
                            double sum_divide_com = sum_vector(divide_com[m]);
                            if (divide_com[m].size() <= k) {
                            }
                            else {
                                if (output_res.size() + Lc.size() < r) {
                                    Lc.push_back(make_pair(divide_com[m], sum_divide_com));
                                    if (sum_divide_com < top_r_weight) {
                                        top_r_weight = sum_divide_com;
                                    }
                                }
                                else {
                                    if (sum_divide_com > top_r_weight) {
                                        Lc.push_back(make_pair(divide_com[m], sum_divide_com));
                                    }
                                }
                            }
                        }
                        batch_add_edges(batch_del_edges_vec);
                    }
                }
            }
            output_res.insert(output_res.end(), Lc.begin(), Lc.end());
            sort(output_res.begin(), output_res.end(),
                    [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
                        return left.second > right.second;
            });
            if (output_res.size() > r) {
                output_res.assign(output_res.begin(), output_res.begin() + r);
            }
            if (output_res.size() > 0) {
                top_r_weight = output_res.back().second;
            }
        }
    }

    cout << "naive sum top-r end" << endl;

    // output result
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/naive_sum_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                     " " << output_res[i].first.size() << endl;
    }
    outfile.close();

    delete [] is_exist;
}

void TopRComm::improved_sum_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            k_core_vec.push_back(i);
        }
    }

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            if (L.size() < r) {
                L.push_back(make_pair(res[i], sum_vec));
                if (sum_vec < top_r_weight) {
                    top_r_weight = sum_vec;
                    top_r_index = L.size() - 1;
                }
                if (sum_vec > max_weight) {
                    max_weight = sum_vec;
                    max_index = L.size() - 1;
                }
            }
            else {
                if (sum_vec > top_r_weight) {
                    L.erase(L.begin() + top_r_index);
                    L.push_back(make_pair(res[i], sum_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    top_r_index = 0;
                    max_weight = std::numeric_limits<double>::min();
                    max_index = 0;
                    for (unsigned int j = 0; j < L.size(); j++) {
                        if (top_r_weight > L[j].second) {
                            top_r_weight = L[j].second;
                            top_r_index = j;
                        }
                        if (max_weight < L[j].second) {
                            max_weight = L[j].second;
                            max_index = j;
                        }
                    }
                }
            }
        }
    }

    cout << "improve sum top-r begin" << endl;

    while (output_res.size() < r) {
        output_res.push_back(L[max_index]);
        if (output_res.size() == r) {
            break;
        }
        add_vec.assign(L[max_index].first.begin(), L[max_index].first.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sum_vec = sum_vector(add_vec);
        node_weight.clear();
        for (unsigned int i = 0; i < L[max_index].first.size(); i++) {
            node_weight.push_back(make_pair(i, sum_vec - weight[add_vec[i]]));
        }
        sort(node_weight.begin(), node_weight.end(), [](const pair<int, double> &left, const pair<int, double> &right){
            return left.second > right.second;
        });

        for (unsigned int i = 0; i < node_weight.size(); i++) {
            delete_node = add_vec[node_weight[i].first];
            add_vec.erase(add_vec.begin() + node_weight[i].first);
            sum_vec = sum_vector(add_vec);
            if (add_vec.size() <= k) {
            }
            else if (sum_vec <= top_r_weight && L.size() >= r + 1) {
                break;
            }
            else {
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    double sum_divide_com = sum_vector(divide_com[j]);
                    if (divide_com[j].size() <= k) {
                    }
                    else {
                        if (L.size() <= r) {
                            L.push_back(make_pair(divide_com[j], sum_divide_com));
                            if (sum_divide_com < top_r_weight) {
                                top_r_weight = sum_divide_com;
                                top_r_index = L.size() - 1;
                            }
                        }
                        else {
                            if (sum_divide_com > top_r_weight) {
                                L.erase(L.begin() + top_r_index);
                                L.push_back(make_pair(divide_com[j], sum_divide_com));
                                top_r_weight = std::numeric_limits<double>::max();
                                top_r_index = 0;
                                for (unsigned int m = 0; m < L.size(); m++) {
                                    if (top_r_weight > L[m].second) {
                                        top_r_weight = L[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                        }
                    }
                }
                batch_add_edges(single_del_edges_vec);
            }
            add_vec.insert(add_vec.begin() + node_weight[i].first, delete_node);
        }
        batch_add_edges(batch_del_edges_vec);
        L.erase(L.begin() + max_index);
        max_weight = std::numeric_limits<double>::min();
        max_index = 0;
        for (unsigned int i = 0;  i < L.size(); i++) {
            if (L[i].second > max_weight) {
                max_weight = L[i].second;
                max_index = i;
            }
        }
        cout << "current size: " << output_res.size() << endl;
    }

    cout << "improve sum top-r end" << endl;

//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_sum_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
    }
    outfile.close();
}

void TopRComm::approx_sum_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double eps) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double upper_bound_weight;
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            k_core_vec.push_back(i);
        }
    }

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            if (L.size() < r) {
                L.push_back(make_pair(res[i], sum_vec));
                if (sum_vec < top_r_weight) {
                    top_r_weight = sum_vec;
                    top_r_index = L.size() - 1;
                }
                if (sum_vec > max_weight) {
                    max_weight = sum_vec;
                    max_index = L.size() - 1;
                }
            }
            else {
                if (sum_vec > top_r_weight) {
                    L.erase(L.begin() + top_r_index);
                    L.push_back(make_pair(res[i], sum_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    top_r_index = 0;
                    max_weight = std::numeric_limits<double>::min();
                    max_index = 0;
                    for (unsigned int j = 0; j < L.size(); j++) {
                        if (top_r_weight > L[j].second) {
                            top_r_weight = L[j].second;
                            top_r_index = j;
                        }
                        if (max_weight < L[j].second) {
                            max_weight = L[j].second;
                            max_index = j;
                        }
                    }
                }
            }
        }
    }

    upper_bound_weight = max_weight;
    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (sum_vec >= upper_bound_weight * (1 - eps)) {
            output_res.push_back(make_pair(res[i], sum_vec));
        }
    }

    cout << "approx sum top-r begin" << endl;

    while (output_res.size() < r) {
        upper_bound_weight = max_weight;
//        output_res.push_back(L[max_index]);
        if (output_res.size() == r) {
            break;
        }
        add_vec.assign(L[max_index].first.begin(), L[max_index].first.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sum_vec = sum_vector(add_vec);
        node_weight.clear();
        for (unsigned int i = 0; i < L[max_index].first.size(); i++) {
            node_weight.push_back(make_pair(i, sum_vec - weight[add_vec[i]]));
        }
        sort(node_weight.begin(), node_weight.end(), [](const pair<int, double> &left, const pair<int, double> &right){
            return left.second > right.second;
        });

        for (unsigned int i = 0; i < node_weight.size(); i++) {
            if (output_res.size() == r) {
                break;
            }
            delete_node = add_vec[node_weight[i].first];
            add_vec.erase(add_vec.begin() + node_weight[i].first);
            sum_vec = sum_vector(add_vec);
            if (add_vec.size() <= k) {
            }
            else if (sum_vec <= top_r_weight && L.size() >= r + 1) {
                break;
            }
            else {
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    double sum_divide_com = sum_vector(divide_com[j]);
                    if (divide_com[j].size() <= k) {
                    }
                    else {
                        // satisfy approximation constraint output it
                        if (sum_divide_com >= upper_bound_weight * (1 - eps)) {
                            output_res.push_back(make_pair(divide_com[j], sum_divide_com));
                        }
                        if (output_res.size() == r) {
                            break;
                        }

                        if (L.size() <= r) {
                            L.push_back(make_pair(divide_com[j], sum_divide_com));
                            if (sum_divide_com < top_r_weight) {
                                top_r_weight = sum_divide_com;
                                top_r_index = L.size() - 1;
                            }
                        }
                        else {
                            if (sum_divide_com > top_r_weight) {
                                L.erase(L.begin() + top_r_index);
                                L.push_back(make_pair(divide_com[j], sum_divide_com));
                                top_r_weight = std::numeric_limits<double>::max();
                                top_r_index = 0;
                                for (unsigned int m = 0; m < L.size(); m++) {
                                    if (top_r_weight > L[m].second) {
                                        top_r_weight = L[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                        }
                    }
                }
                batch_add_edges(single_del_edges_vec);
            }
            add_vec.insert(add_vec.begin() + node_weight[i].first, delete_node);
        }
        batch_add_edges(batch_del_edges_vec);
        L.erase(L.begin() + max_index);
        max_weight = std::numeric_limits<double>::min();
        max_index = 0;
        for (unsigned int i = 0;  i < L.size(); i++) {
            if (L[i].second > max_weight) {
                max_weight = L[i].second;
                max_index = i;
            }
        }
        cout << "current size: " << output_res.size() << endl;
    }

    cout << "approx sum top-r end" << endl;

//    sort(output_res.begin(), output_res.end(), [](const pair<vector<int>, double> & left,
//            const pair<vector<int>, double> & right){
//        return left.second > right.second;
//    });
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }

    endTime = clock();
    cout << "k: " << k << " r: " << r << " eps: " << eps << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string eps_name = "eps_" + std::to_string(eps) + "_";
    filename = "O3result/approx_sum_topr_" + graphname + "_" + k_name + r_name + eps_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    sort(output_res.begin(), output_res.end(),
            [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
        return left.second > right.second;
    });
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                << " " << output_res[i].first.size() << endl;
    }
    outfile.close();
}

/// sort order delete by decreasing order
void TopRComm::improved_greedy_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    map<int, int> add_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    int top_del_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            L.push_back(make_pair(res[i], avg_vec));
            output_res.push_back(make_pair(res[i], avg_vec));
        }
    }

//    /// test code
//    for (unsigned int i = 0; i < res.size(); i++) {
//        cout << "res size: " << res[i].size() << endl;
//    }
//    /// test code

    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            k_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved greedy avg top-r begin" << endl;

    for (unsigned int i = 0; i < L.size(); i++) {
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        /// store the location of element
        add_vec_loc.clear();
        for (unsigned int j = 0; j < add_vec.size(); j++) {
            add_vec_loc[add_vec[j]] = j;
        }
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });
        batch_iter_del_edges_vec.clear();
        top_del_node = -1;
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            delete_node = exist_weight[i][j].first;
            single_del_edges_vec = single_del_edges(delete_node);
            batch_iter_del_edges_vec.insert(batch_iter_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            int del_elem_loc = add_vec_loc[delete_node];
            add_vec_loc[add_vec.back()] = del_elem_loc;
            swap(add_vec[del_elem_loc], add_vec.back());
            add_vec.pop_back();
            if (top_del_node == -1) {
                top_del_node = greedy_unsatisfied_node(k);
            }
            else if (top_del_node == delete_node){
                top_del_node = greedy_unsatisfied_node(k);
            }
            if (top_del_node == -1 && add_vec.size() > k) {
                update_weight = avg_vector(add_vec);
                output_res.push_back(make_pair(add_vec, update_weight));
            }
//            update_weight = avg_vector(add_vec);
//            if (is_connected_k_core_graph(k)) {
//                output_res.push_back(make_pair(add_vec, update_weight));
//            }
        }
        batch_add_edges(batch_iter_del_edges_vec);
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "improved greedy avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

//    /// test code
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        cout << i + 1 << " " << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }
//    /// test code

    cout << "output res size: " << output_res.size() << endl;

    // output result
//    if (output_res.size() >= r) {
//        for (unsigned int i = 0; i < r; i++) {
//            cout << i + 1 << " " << setprecision(8) << output_res[i].second
//                 << " " << output_res[i].first.size() << endl;
//        }
//    }
//    else {
//        for (unsigned int i = 0; i < output_res.size(); i++) {
//            cout << i + 1 << " " << setprecision(8) << output_res[i].second
//                 << " " << output_res[i].first.size() << endl;
//        }
//    }

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_greedy_avg_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    outfile.close();
}

/// random order delete
void TopRComm::improved_random_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    map<int, int> add_vec_loc;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_order_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    int delete_node;
    int top_del_idx;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            L.push_back(make_pair(res[i], avg_vec));
            output_res.push_back(make_pair(res[i], avg_vec));
        }
    }

//    /// test code
//    for (unsigned int i = 0; i < res.size(); i++) {
//        cout << "res size: " << res[i].size() << endl;
//    }
//    /// test code

    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            k_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved random avg top-r begin" << endl;

    for (unsigned int i = 0; i < L.size(); i++) {
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        del_nodes_order_vec.assign(L[i].first.begin(), L[i].first.end());
        random_shuffle(del_nodes_order_vec.begin(), del_nodes_order_vec.end());
        /// store the location of the delete element
        del_nodes_order_vec_loc.clear();
        for (unsigned int j = 0; j < del_nodes_order_vec.size(); j++) {
            del_nodes_order_vec_loc[del_nodes_order_vec[j]] = j;
        }
        /// store the location of element
        add_vec_loc.clear();
        for (unsigned int j = 0; j < add_vec.size(); j++) {
            add_vec_loc[add_vec[j]] = j;
        }
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        batch_iter_del_edges_vec.clear();
        top_del_idx = -1;
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            delete_node = del_nodes_order_vec[j];
            single_del_edges_vec = single_del_edges(delete_node);
            batch_iter_del_edges_vec.insert(batch_iter_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            int del_elem_loc = add_vec_loc[delete_node];
            add_vec_loc[add_vec.back()] = del_elem_loc;
            swap(add_vec[del_elem_loc], add_vec.back());
            add_vec.pop_back();
            if (top_del_idx == -1) {
                top_del_idx = random_unsatisfied_node(del_nodes_order_vec, k);
            }
            else if (top_del_idx == j){
                top_del_idx = random_unsatisfied_node(del_nodes_order_vec, k);
            }
            if (top_del_idx == -1 && add_vec.size() > k) {
                update_weight = avg_vector(add_vec);
                output_res.push_back(make_pair(add_vec, update_weight));
            }
//            update_weight = avg_vector(add_vec);
//            if (is_connected_k_core_graph(k)) {
//                output_res.push_back(make_pair(add_vec, update_weight));
//            }
        }
        batch_add_edges(batch_iter_del_edges_vec);
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "improved random avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

//    /// test code
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        cout << i + 1 << " " << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }
//    /// test code

    cout << "output res size: " << output_res.size() << endl;

    // output result
//    if (output_res.size() >= r) {
//        for (unsigned int i = 0; i < r; i++) {
//            cout << i + 1 << " " << setprecision(8) << output_res[i].second
//                 << " " << output_res[i].first.size() << endl;
//        }
//    }
//    else {
//        for (unsigned int i = 0; i < output_res.size(); i++) {
//            cout << i + 1 << " " << setprecision(8) << output_res[i].second
//                 << " " << output_res[i].first.size() << endl;
//        }
//    }

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_random_avg_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    outfile.close();
}

/// delete the node which could improve the influence value
void TopRComm::improved_climb_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int curr_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved climb avg top-r begin" << endl;

    while (Q.size() != 0) {
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
//        double curr_inf_value = avg_vector(add_vec);
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        while (add_vec.size() >= k + 1) {
            delete_node = add_vec[0];
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                curr_node = add_vec[i];
                if (weight[delete_node] > weight[curr_node]) {
                    delete_node = curr_node;
                }
            }
            single_del_edges_vec.clear();
            single_del_edges_vec = single_del_edges(delete_node);
            batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            vector<vector<int>> divide_com = compute_connected_k_core(k);
            com_max_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            for (unsigned int i = 0; i < divide_com.size(); i++) {
                avg_vec = avg_vector(divide_com[i]);
                if (avg_vec > com_max_inf_value) {
                    com_max_inf_value = avg_vec;
                    max_inf_com = divide_com[i];
                }
            }
            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com_set.begin(), max_inf_com_set.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                max_inf_value = com_max_inf_value;
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "improved climb avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_climb_avg_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    outfile.close();
}

/// delete the node which could improve the influence value most
void TopRComm::improved_double_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int max_delete_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved double avg top-r begin" << endl;

    while (Q.size() != 0) {
        double curr_inf_value;
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        while (add_vec.size() >= k + 1) {
            max_diff_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            curr_inf_value = avg_vector(add_vec);
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                delete_node = add_vec[i];
                single_del_edges_vec.clear();
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                vector<int> local_max_inf_com;
                com_max_inf_value = std::numeric_limits<double>::min();
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    avg_vec = avg_vector(divide_com[j]);
                    if (avg_vec > com_max_inf_value) {
                        com_max_inf_value = avg_vec;
                        local_max_inf_com = divide_com[j];
                    }
                }
                diff_inf_value = com_max_inf_value - curr_inf_value;
                if (diff_inf_value > max_diff_inf_value) {
                    max_inf_value = com_max_inf_value;
                    max_diff_inf_value = diff_inf_value;
                    max_delete_node = delete_node;
                    max_inf_com = local_max_inf_com;
                }
                batch_add_edges(single_del_edges_vec);
            }
            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com.begin(), max_inf_com.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                single_del_edges_vec = single_del_edges(max_delete_node);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        single_del_edges_vec.begin(), single_del_edges_vec.end());
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "improved double avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_double_avg_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    outfile.close();
}

void TopRComm::cut_double_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    map<int, int> add_vec_loc;
    vector<int> cut_vec;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int max_delete_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "cut double avg top-r begin" << endl;

    while (Q.size() != 0) {
        double curr_inf_value;
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        while (add_vec.size() >= k + 1) {
            curr_inf_value = avg_vector(add_vec);
            /// store the location of element
            add_vec_loc.clear();
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                add_vec_loc[add_vec[i]] = i;
            }

            /// compute cut vertex
            bool *cut_vertex_flag = new bool [graph_size];
            for (unsigned int i = 0; i < graph_size; i++) {
                cut_vertex_flag[i] = false;
            }
            cut_vec = compute_cut_vertex();
            for (unsigned int i = 0; i < cut_vec.size(); i++) {
                cut_vertex_flag[cut_vec[i]] = true;
            }

            max_diff_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                delete_node = add_vec[i];
                if (cut_vertex_flag[delete_node] == false) {
                    bool flag = true;
                    for (unsigned int j = 0; j < AdjList[delete_node].size(); j++) {
                        if (AdjList[AdjList[delete_node][j]].size() <= k) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        if (weight[delete_node] > curr_inf_value) {
                        }
                        else {
                            int del_elem_loc = add_vec_loc[delete_node];
                            add_vec_loc[add_vec.back()] = del_elem_loc;
                            swap(add_vec[del_elem_loc], add_vec.back());
                            add_vec.pop_back();
                            com_max_inf_value = avg_vector(add_vec);
                            diff_inf_value = com_max_inf_value - curr_inf_value;
                            if (diff_inf_value > max_diff_inf_value) {
                                max_inf_value = com_max_inf_value;
                                max_delete_node = delete_node;
                                max_diff_inf_value = diff_inf_value;
                                max_inf_com = add_vec;
                            }
                            add_vec.push_back(delete_node);
                            add_vec_loc[delete_node] = add_vec.size() - 1;
                        }
                    }
                    else {
                        single_del_edges_vec.clear();
                        single_del_edges_vec = single_del_edges(delete_node);
                        vector<vector<int>> divide_com = compute_connected_k_core(k);
                        vector<int> local_max_inf_com;
                        com_max_inf_value = std::numeric_limits<double>::min();
                        for (unsigned int j = 0; j < divide_com.size(); j++) {
                            avg_vec = avg_vector(divide_com[j]);
                            if (avg_vec > com_max_inf_value) {
                                com_max_inf_value = avg_vec;
                                local_max_inf_com = divide_com[j];
                            }
                        }
                        diff_inf_value = com_max_inf_value - curr_inf_value;
                        if (diff_inf_value > max_diff_inf_value) {
                            max_inf_value = com_max_inf_value;
                            max_delete_node = delete_node;
                            max_diff_inf_value = diff_inf_value;
                            max_inf_com = local_max_inf_com;
                        }
                        batch_add_edges(single_del_edges_vec);
                    }
                }
                else {
                    single_del_edges_vec.clear();
                    single_del_edges_vec = single_del_edges(delete_node);
                    vector<vector<int>> divide_com = compute_connected_k_core(k);
                    vector<int> local_max_inf_com;
                    com_max_inf_value = std::numeric_limits<double>::min();
                    for (unsigned int j = 0; j < divide_com.size(); j++) {
                        avg_vec = avg_vector(divide_com[j]);
                        if (avg_vec > com_max_inf_value) {
                            com_max_inf_value = avg_vec;
                            local_max_inf_com = divide_com[j];
                        }
                    }
                    diff_inf_value = com_max_inf_value - curr_inf_value;
                    if (diff_inf_value > max_diff_inf_value) {
                        max_inf_value = com_max_inf_value;
                        max_delete_node = delete_node;
                        max_diff_inf_value = diff_inf_value;
                        max_inf_com = local_max_inf_com;
                    }
                    batch_add_edges(single_del_edges_vec);
                }
            }

            delete [] cut_vertex_flag;

            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com.begin(), max_inf_com.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                                           batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "cut double avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/cut_double_avg_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
    }
    outfile.close();
}

void TopRComm::naive_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> output_res;
    vector<pair<vector<int>, double>> Lc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    double sum_vec;
    double top_r_weight;
    int top_r_index;
    bool flag = false; // determine whether duplicate
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    bool *is_exist = new bool[graph_size];
    for (int i = 0; i < graph_size; i++) {
        is_exist[i] = false;
    }
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            output_res.push_back(make_pair(res[i], sum_vec));
        }
    }

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int> ,double> & right){
             return left.second > right.second;
         });

    if (output_res.size() > r) {
        output_res.assign(output_res.begin(), output_res.begin() + r);
    }
    top_r_weight = output_res.back().second;
    top_r_index = output_res.size() - 1;

    cout << "init: " << output_res.size() << endl;

    for (unsigned int i = 0; i < output_res.size(); i++) {
        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
            is_exist[output_res[i].first[j]] = true;
            k_core_vec.push_back(output_res[i].first[j]);
        }
    }

    Lc.assign(output_res.begin(), output_res.end());
    cout << "naive cons sum top-r begin" << endl;

    for (int i = 0; i < graph_size; i++) {
        if (is_exist[i] == true) {
            for (int j = 0; j < output_res.size(); j++) {
                if (output_res[j].first.size() <= k) {
                }
                else if (output_res[j].second - weight[i] <= top_r_weight && Lc.size() >= r){
                }
                else {
                    add_vec.assign(output_res[j].first.begin(), output_res[j].first.end());
                    bool flag = false;
                    for (unsigned int m = 0; m < add_vec.size(); m++) {
                        if (add_vec[m] == i) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag == false) {
                    }
                    else {
                        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
                        set<int> add_set(add_vec.begin(), add_vec.end());
                        set<int> del_nodes_set;
                        del_nodes_vec.clear();
                        set_difference(k_core_set.begin(), k_core_set.end(),
                                       add_set.begin(), add_set.end(),
                                       inserter(del_nodes_set, del_nodes_set.begin()));
                        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
                        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
                        single_del_edges_vec = single_del_edges(i);
                        vector<vector<int>> divide_com = compute_connected_k_core(k);
                        for (unsigned int m = 0; m < divide_com.size(); m++) {
                            double sum_divide_com = sum_vector(divide_com[m]);
                            if (divide_com[m].size() <= k) {
                            }
                            else {
                                if (Lc.size() < r) {
                                    Lc.push_back(make_pair(divide_com[m], sum_divide_com));
                                    if (sum_divide_com < top_r_weight) {
                                        top_r_weight = sum_divide_com;
                                        top_r_index = Lc.size() - 1;
                                    }
                                }
                                else {
                                    if (sum_divide_com > top_r_weight) {
                                        Lc.push_back(make_pair(divide_com[m], sum_divide_com));
                                        Lc.erase(Lc.begin() + top_r_index);
                                        top_r_weight = std::numeric_limits<double>::max();
                                        top_r_index = 0;
                                        for (unsigned int index = 0; index < Lc.size(); index++) {
                                            if (top_r_weight > Lc[index].second) {
                                                top_r_weight = Lc[index].second;
                                                top_r_index = index;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        batch_add_edges(single_del_edges_vec);
                        batch_add_edges(batch_del_edges_vec);
                    }
                }
            }
            output_res.assign(Lc.begin(), Lc.end());
        }
    }

    output_res.assign(Lc.begin(), Lc.end());
    sort(output_res.begin(), output_res.end(),
            [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
        return left.second > right.second;
    });

    cout << "naive cons sum top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/naive_cons_sum_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
    }
    outfile.close();

    delete [] is_exist;
}

void TopRComm::improved_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
//    compute_max_k_core(degree_c);
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    int num_nodes = 0;
    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            k_core_vec.push_back(i);
            num_nodes = num_nodes + 1;
        }
    }
    cout << "num of nodes: " << num_nodes << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;
    cout << "res size: " << res.size() << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            if (L.size() < r) {
                L.push_back(make_pair(res[i], sum_vec));
                if (sum_vec < top_r_weight) {
                    top_r_weight = sum_vec;
                    top_r_index = L.size() - 1;
                }
                if (sum_vec > max_weight) {
                    max_weight = sum_vec;
                    max_index = L.size() - 1;
                }
            }
            else {
                if (sum_vec > top_r_weight) {
                    L.erase(L.begin() + top_r_index);
                    L.push_back(make_pair(res[i], sum_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    top_r_index = 0;
                    max_weight = std::numeric_limits<double>::min();
                    max_index = 0;
                    for (unsigned int j = 0; j < L.size(); j++) {
                        if (top_r_weight > L[j].second) {
                            top_r_weight = L[j].second;
                            top_r_index = j;
                        }
                        if (max_weight < L[j].second) {
                            max_weight = L[j].second;
                            max_index = j;
                        }
                    }
                }
            }
        }
    }

    cout << "improve cons sum top-r begin" << endl;

    while (output_res.size() < r) {
        output_res.push_back(L[max_index]);
        if (output_res.size() == r) {
            break;
        }
        add_vec.assign(L[max_index].first.begin(), L[max_index].first.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sum_vec = sum_vector(add_vec);
        node_weight.clear();
        for (unsigned int i = 0; i < add_vec.size(); i++) {
            node_weight.push_back(make_pair(i, sum_vec - weight[add_vec[i]]));
        }
        sort(node_weight.begin(), node_weight.end(), [](const pair<int, double> &left, const pair<int, double> &right){
            return left.second > right.second;
        });

        for (unsigned int i = 0; i < node_weight.size(); i++) {
            delete_node = add_vec[node_weight[i].first];
            add_vec.erase(add_vec.begin() + node_weight[i].first);
            sum_vec = sum_vector(add_vec);
            if (add_vec.size() <= k) {
            }
            else if (sum_vec <= top_r_weight && L.size() >= r) {
                break;
            }
            else {
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    double sum_divide_com = sum_vector(divide_com[j]);
                    if (divide_com[j].size() <= k) {
                    }
                    else {
                        if (L.size() < r) {
                            L.push_back(make_pair(divide_com[j], sum_divide_com));
                            if (sum_divide_com < top_r_weight) {
                                top_r_weight = sum_divide_com;
                                top_r_index = L.size() - 1;
                            }
                        }
                        else {
                            if (sum_divide_com > top_r_weight) {
                                L.erase(L.begin() + top_r_index);
                                L.push_back(make_pair(divide_com[j], sum_divide_com));
                                top_r_weight = std::numeric_limits<double>::max();
                                top_r_index = 0;
                                for (unsigned int m = 0; m < L.size(); m++) {
                                    if (top_r_weight > L[m].second) {
                                        top_r_weight = L[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                        }
                    }
                }
                batch_add_edges(single_del_edges_vec);
            }
            add_vec.insert(add_vec.begin() + node_weight[i].first, delete_node);
        }
        batch_add_edges(batch_del_edges_vec);
        L.erase(L.begin() + max_index);
        max_weight = std::numeric_limits<double>::min();
        max_index = 0;
        for (unsigned int i = 0; i < L.size(); i++) {
            if (L[i].second > max_weight) {
                max_weight = L[i].second;
                max_index = i;
            }
        }
        cout << "current size: " << output_res.size() << endl;
    }

    cout << "improve cons sum top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    filename = "O3result/improve_cons_sum_topr_" + graphname + "_" + k_name + r_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
    }
//    for (unsigned int i = 0; i < output_res[output_res.size() - 1].first.size(); i++) {
//        outfile << output_res[output_res.size() - 1].first[i] << " ";
//    }
//    outfile << endl;
    outfile.close();
}

void TopRComm::approx_cons_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const double eps) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double upper_bound_weight;
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            k_core_vec.push_back(i);
        }
    }

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            if (L.size() < r) {
                L.push_back(make_pair(res[i], sum_vec));
                if (sum_vec < top_r_weight) {
                    top_r_weight = sum_vec;
                    top_r_index = L.size() - 1;
                }
                if (sum_vec > max_weight) {
                    max_weight = sum_vec;
                    max_index = L.size() - 1;
                }
            }
            else {
                if (sum_vec > top_r_weight) {
                    L.erase(L.begin() + top_r_index);
                    L.push_back(make_pair(res[i], sum_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    top_r_index = 0;
                    max_weight = std::numeric_limits<double>::min();
                    max_index = 0;
                    for (unsigned int j = 0; j < L.size(); j++) {
                        if (top_r_weight > L[j].second) {
                            top_r_weight = L[j].second;
                            top_r_index = j;
                        }
                        if (max_weight < L[j].second) {
                            max_weight = L[j].second;
                            max_index = j;
                        }
                    }
                }
            }
        }
    }

    upper_bound_weight = max_weight;
    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (sum_vec >= upper_bound_weight * (1 - eps)) {
            output_res.push_back(make_pair(res[i], sum_vec));
        }
    }

    cout << "approx cons sum top-r begin" << endl;

    while (output_res.size() < r) {
        upper_bound_weight = max_weight;
        if (output_res.size() == r) {
            break;
        }
        add_vec.assign(L[max_index].first.begin(), L[max_index].first.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sum_vec = sum_vector(add_vec);
        node_weight.clear();
        for (unsigned int i = 0; i < add_vec.size(); i++) {
            node_weight.push_back(make_pair(i, sum_vec - weight[add_vec[i]]));
        }
        sort(node_weight.begin(), node_weight.end(), [](const pair<int, double> &left, const pair<int, double> &right){
            return left.second > right.second;
        });

        for (unsigned int i = 0; i < node_weight.size(); i++) {
            if (output_res.size() == r) {
                break;
            }
            delete_node = add_vec[node_weight[i].first];
            add_vec.erase(add_vec.begin() + node_weight[i].first);
            sum_vec = sum_vector(add_vec);
            if (add_vec.size() <= k) {
            }
            else if (sum_vec <= top_r_weight && L.size() >= r + 1) {
                break;
            }
            else {
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    double sum_divide_com = sum_vector(divide_com[j]);
                    if (divide_com[j].size() <= k) {
                    }
                    else {
                        // satisfy approximation constraint output it
                        if (sum_divide_com >= upper_bound_weight * (1 - eps)) {
                            output_res.push_back(make_pair(divide_com[j], sum_divide_com));
                        }
                        if (output_res.size() == r) {
                            break;
                        }

                        if (L.size() <= r) {
                            L.push_back(make_pair(divide_com[j], sum_divide_com));
                            if (sum_divide_com < top_r_weight) {
                                top_r_weight = sum_divide_com;
                                top_r_index = L.size() - 1;
                            }
                        }
                        else {
                            if (sum_divide_com > top_r_weight) {
                                L.erase(L.begin() + top_r_index);
                                L.push_back(make_pair(divide_com[j], sum_divide_com));
                                top_r_weight = std::numeric_limits<double>::max();
                                top_r_index = 0;
                                for (unsigned int m = 0; m < L.size(); m++) {
                                    if (top_r_weight > L[m].second) {
                                        top_r_weight = L[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                        }
                    }
                }
                batch_add_edges(single_del_edges_vec);
            }
            add_vec.insert(add_vec.begin() + node_weight[i].first, delete_node);
        }
        batch_add_edges(batch_del_edges_vec);
        L.erase(L.begin() + max_index);
        max_weight = std::numeric_limits<double>::min();
        max_index = 0;
        for (unsigned int i = 0;  i < L.size(); i++) {
            if (L[i].second > max_weight) {
                max_weight = L[i].second;
                max_index = i;
            }
        }
        cout << "current size: " << output_res.size() << endl;
    }

    cout << "approx cons sum top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " eps: " << eps << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string eps_name = "eps_" + std::to_string(eps) + "_";
    filename = "O3result/approx_cons_sum_topr_" + graphname + "_" + k_name + r_name + eps_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                << " " << output_res[i].first.size() << endl;
    }
    outfile.close();
}

void TopRComm::non_cons_sum_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<pair<int, double>> node_weight; // node_id to update_weight
    double update_weight; // updated weight
    double sum_vec;
    double max_weight = std::numeric_limits<double>::min();
    int max_index = 0;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    int delete_node;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
//    compute_max_k_core(degree_c);
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    int num_nodes = 0;
    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            k_core_vec.push_back(i);
            num_nodes = num_nodes + 1;
        }
    }
    cout << "num of nodes: " << num_nodes << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;
    cout << "res size: " << res.size() << endl;

    if (res.size() == 0) {
        return;
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            if (L.size() < r) {
                L.push_back(make_pair(res[i], sum_vec));
                if (sum_vec < top_r_weight) {
                    top_r_weight = sum_vec;
                    top_r_index = L.size() - 1;
                }
                if (sum_vec > max_weight) {
                    max_weight = sum_vec;
                    max_index = L.size() - 1;
                }
            }
            else {
                if (sum_vec > top_r_weight) {
                    L.erase(L.begin() + top_r_index);
                    L.push_back(make_pair(res[i], sum_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    top_r_index = 0;
                    max_weight = std::numeric_limits<double>::min();
                    max_index = 0;
                    for (unsigned int j = 0; j < L.size(); j++) {
                        if (top_r_weight > L[j].second) {
                            top_r_weight = L[j].second;
                            top_r_index = j;
                        }
                        if (max_weight < L[j].second) {
                            max_weight = L[j].second;
                            max_index = j;
                        }
                    }
                }
            }
        }
    }

    cout << "non cons sum top-r begin" << endl;

    output_res = L;
    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "non cons sum top-r end" << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/non_cons_sum_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    for (unsigned int i = 0; i < output_res.size(); i++) {
        outfile << i + 1 << " " << setprecision(8) << output_res[i].second <<
                " " << output_res[i].first.size() << endl;
        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
            outfile << output_res[i].first[j] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void TopRComm::climb_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    map<int, int> add_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<int> local_Q;
    vector<pair<int, double>> local_candidate;
    vector<int> C;
    double top_r_weight;
    int top_r_index;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    batch_del_edges_vec.clear();
    top_r_index = 0;
    top_r_weight = std::numeric_limits<double>::min();
    for (unsigned int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(i);
            in_vec[i] = true;
            while (local_Q.size() != 0 && local_candidate.size() <= s) {
                u = local_Q.front();
                local_candidate.push_back(make_pair(u, weight[u]));
                if (local_Q.size() + local_candidate.size() < s) {
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(v);
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
	    //int temp_add = s <  local_candidate.size() ? s : local_candidate.size();
	    //cout << "vertex : "  << i << "\t with local size" << local_candidate.size() << endl;
            sort(local_candidate.begin(), local_candidate.end(),//begin()+temp_add,//end(),
                    [](const pair<int, double> & left, const pair<int, double> & right){
                return left.second > right.second;
            });
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(local_candidate[j].first);
                }
                in_vec[local_candidate[j].first] = false;
            }
            for (unsigned int j = 0; j < C.size(); j++) {
                core_exist[C[j]] = true;
            }
            while (C.size() >= k + 1 && sum_vector(C) > top_r_weight) {
                if (determined_k_core(C, core_exist, k)) {
//                    batch_del_edges_local_vec.clear();
//                    batch_del_edges_local_vec = batch_del_edges(C);
//                    batch_del_edges_vec.insert(batch_del_edges_vec.end(),
//                            batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                    if (output_res.size() >= r) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(C, sum_vector(C)));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (output_res[j].second < top_r_weight) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                    else if (output_res.size() == r - 1) {
                        output_res.push_back(make_pair(C, sum_vector(C)));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (output_res[j].second < top_r_weight) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                    else {
                        output_res.push_back(make_pair(C, sum_vector(C)));
                    }
                    for (unsigned int j = 0; j < C.size(); j++) {
                        core_exist[C[j]] = false;
                    }
                    break;
                }
                else {
                    core_exist[C.back()] = false;
                    C.pop_back();
                }
            }
        }
    }
//    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "climb size sum top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/climb_size_sum_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::degree_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    map<int, int> add_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<pair<int, int>> local_Q;
    vector<tuple<int, int, double>> local_candidate;
    vector<int> C;
    double com_max_inf_value;
    int com_max_idx;
    int delete_node;
    int curr_node;
    double sum_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    int hier;
    int curr_hier;
    vector<pair<int, double>> new_hier;
    batch_del_edges_vec.clear();
    for (unsigned int i = 0; i < graph_size; i++) {
        curr_hier = 0;
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(make_pair(i, 0));
            in_vec[i] = true;
            while (local_Q.size() != 0) {
                u = local_Q.front().first;
                hier = local_Q.front().second;
                if (hier > curr_hier) {
                    if (local_candidate.size() < s) {
                        curr_hier = hier;
                        local_candidate.push_back(make_tuple(u, hier, weight[u]));
                        for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                            v = AdjList[u][j];
                            if (AdjList[v].size() >= k) {
                                if (in_vec[v] == false) {
                                    local_Q.push(make_pair(v, hier + 1));
                                    in_vec[v] = true;
                                }
                            }
                        }
                    }
                }
                else {
                    local_candidate.push_back(make_tuple(u, hier, weight[u]));
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(make_pair(v, hier + 1));
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
            sort(local_candidate.begin(), local_candidate.end(),
                 [](const tuple<int, int, double> & left, const tuple<int, int, double> & right){
                     if (get<1>(left) == get<1>(right)) {
                         return get<2>(left) > get<2>(right);
                     }
                     return get<1>(left) < get<1>(right);
                 });
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(get<0>(local_candidate[j]));
                }
                in_vec[get<0>(local_candidate[j])] = false;
            }
            for (unsigned int j = 0; j < C.size(); j++) {
                core_exist[C[j]] = true;
            }
            while (C.size() >= k + 1) {
                if (determined_k_core(C, core_exist, k)) {
                    batch_del_edges_local_vec.clear();
                    batch_del_edges_local_vec = batch_del_edges(C);
                    batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                            batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                    output_res.push_back(make_pair(C, sum_vector(C)));
                    for (unsigned int j = 0; j < C.size(); j++) {
                        core_exist[C[j]] = false;
                    }
                    break;
                }
                core_exist[C.back()] = false;
                C.pop_back();
            }
        }
    }
    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "degree size sum top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/degree_size_sum_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::casual_size_sum_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<int> local_Q;
    vector<pair<int, double>> local_candidate;
    vector<int> C;
    double top_r_weight;
    int top_r_index;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    batch_del_edges_vec.clear();
    top_r_index = 0;
    top_r_weight = std::numeric_limits<double>::min();
    for (unsigned int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(i);
            in_vec[i] = true;
            while (local_Q.size() != 0){
                u = local_Q.front();
                local_candidate.push_back(make_pair(u, weight[u]));
                if (local_Q.size() + local_candidate.size() < s) {
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(v);
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
//            srand ( unsigned ( time(0) ) );
//            random_shuffle(local_candidate.begin(), local_candidate.end());
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(local_candidate[j].first);
                }
                in_vec[local_candidate[j].first] = false;
            }
            for (unsigned int j = 0; j < C.size(); j++) {
                core_exist[C[j]] = true;
            }
            while (C.size() >= k + 1 && sum_vector(C) > top_r_weight) {
                if (determined_k_core(C, core_exist, k)) {
//                    batch_del_edges_local_vec.clear();
//                    batch_del_edges_local_vec = batch_del_edges(C);
//                    batch_del_edges_vec.insert(batch_del_edges_vec.end(),
//                            batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                    if (output_res.size() >= r) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(C, sum_vector(C)));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (output_res[j].second < top_r_weight) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                    else if (output_res.size() == r - 1) {
                        output_res.push_back(make_pair(C, sum_vector(C)));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (output_res[j].second < top_r_weight) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                    else {
                        output_res.push_back(make_pair(C, sum_vector(C)));
                    }
                    for (unsigned int j = 0; j < C.size(); j++) {
                        core_exist[C[j]] = false;
                    }
                    break;
                }
                else {
                    core_exist[C.back()] = false;
                    C.pop_back();
                }
            }
        }
    }
//    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "casual size sum top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/casual_size_sum_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::improved_cons_greedy_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    map<int, int> add_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    int top_del_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            L.push_back(make_pair(res[i], avg_vec));
            output_res.push_back(make_pair(res[i], avg_vec));
        }
    }

    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            k_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "improved cons greedy avg top-r begin" << endl;

    for (unsigned int i = 0; i < L.size(); i++) {
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        /// store the location of element
        add_vec_loc.clear();
        for (unsigned int j = 0; j < add_vec.size(); j++) {
            add_vec_loc[add_vec[j]] = j;
        }
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });
        batch_iter_del_edges_vec.clear();
        top_del_node = -1;
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            delete_node = exist_weight[i][j].first;
            single_del_edges_vec = single_del_edges(delete_node);
            batch_iter_del_edges_vec.insert(batch_iter_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            int del_elem_loc = add_vec_loc[delete_node];
            add_vec_loc[add_vec.back()] = del_elem_loc;
            swap(add_vec[del_elem_loc], add_vec.back());
            add_vec.pop_back();
            if (top_del_node == -1) {
                top_del_node = greedy_unsatisfied_node(k);
            }
            else if (top_del_node == delete_node){
                top_del_node = greedy_unsatisfied_node(k);
            }
            if (top_del_node == -1 && add_vec.size() > k) {
                update_weight = avg_vector(add_vec);
                output_res.push_back(make_pair(add_vec, update_weight));
            }
        }
        batch_add_edges(batch_iter_del_edges_vec);
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "improved cons greedy avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(0);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/improve_greedy_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[r-1].first.size(); i++) {
            outfile << output_res[r-1].first[i] << " ";
        }
        outfile << endl;
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[output_res.size() - 1].first.size(); i++) {
            outfile << output_res[output_res.size() - 1].first[i] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void TopRComm::improved_cons_random_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> output_res;
    map<int, int> add_vec_loc;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_order_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    int delete_node;
    int top_del_idx;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            L.push_back(make_pair(res[i], avg_vec));
            output_res.push_back(make_pair(res[i], avg_vec));
        }
    }

    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            k_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved cons random avg top-r begin" << endl;

    for (unsigned int i = 0; i < L.size(); i++) {
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        del_nodes_order_vec.assign(L[i].first.begin(), L[i].first.end());
        random_shuffle(del_nodes_order_vec.begin(), del_nodes_order_vec.end());
        /// store the location of the delete element
        del_nodes_order_vec_loc.clear();
        for (unsigned int j = 0; j < del_nodes_order_vec.size(); j++) {
            del_nodes_order_vec_loc[del_nodes_order_vec[j]] = j;
        }
        /// store the location of element
        add_vec_loc.clear();
        for (unsigned int j = 0; j < add_vec.size(); j++) {
            add_vec_loc[add_vec[j]] = j;
        }
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        batch_iter_del_edges_vec.clear();
        top_del_idx = -1;
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            delete_node = del_nodes_order_vec[j];
            single_del_edges_vec = single_del_edges(delete_node);
            batch_iter_del_edges_vec.insert(batch_iter_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            int del_elem_loc = add_vec_loc[delete_node];
            add_vec_loc[add_vec.back()] = del_elem_loc;
            swap(add_vec[del_elem_loc], add_vec.back());
            add_vec.pop_back();
            if (top_del_idx == -1) {
                top_del_idx = random_unsatisfied_node(del_nodes_order_vec, k);
            }
            else if (top_del_idx == j){
                top_del_idx = random_unsatisfied_node(del_nodes_order_vec, k);
            }
            if (top_del_idx == -1 && add_vec.size() > k) {
                update_weight = avg_vector(add_vec);
                output_res.push_back(make_pair(add_vec, update_weight));
            }
        }
        batch_add_edges(batch_iter_del_edges_vec);
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "improved cons random avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/improve_cons_random_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::improved_cons_climb_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int curr_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved cons climb avg top-r begin" << endl;

    while (Q.size() != 0) {
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        while (add_vec.size() >= k + 1) {
            delete_node = add_vec[0];
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                curr_node = add_vec[i];
                if (weight[delete_node] > weight[curr_node]) {
                    delete_node = curr_node;
                }
            }
            single_del_edges_vec.clear();
            single_del_edges_vec = single_del_edges(delete_node);
            batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                    single_del_edges_vec.begin(), single_del_edges_vec.end());
            vector<vector<int>> divide_com = compute_connected_k_core(k);
            com_max_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            for (unsigned int i = 0; i < divide_com.size(); i++) {
                avg_vec = avg_vector(divide_com[i]);
                if (avg_vec > com_max_inf_value) {
                    com_max_inf_value = avg_vec;
                    max_inf_com = divide_com[i];
                }
            }
            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com_set.begin(), max_inf_com_set.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                max_inf_value = com_max_inf_value;
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "improved cons climb avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/improve_cons_climb_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::improved_cons_double_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int max_delete_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved cons double avg top-r begin" << endl;

    while (Q.size() != 0) {
        double curr_inf_value;
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        while (add_vec.size() >= k + 1) {
            max_diff_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            curr_inf_value = avg_vector(add_vec);
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                delete_node = add_vec[i];
                single_del_edges_vec.clear();
                single_del_edges_vec = single_del_edges(delete_node);
                vector<vector<int>> divide_com = compute_connected_k_core(k);
                vector<int> local_max_inf_com;
                com_max_inf_value = std::numeric_limits<double>::min();
                for (unsigned int j = 0; j < divide_com.size(); j++) {
                    avg_vec = avg_vector(divide_com[j]);
                    if (avg_vec > com_max_inf_value) {
                        com_max_inf_value = avg_vec;
                        local_max_inf_com = divide_com[j];
                    }
                }
                diff_inf_value = com_max_inf_value - curr_inf_value;
                if (diff_inf_value > max_diff_inf_value) {
                    max_inf_value = com_max_inf_value;
                    max_diff_inf_value = diff_inf_value;
                    max_delete_node = delete_node;
                    max_inf_com = local_max_inf_com;
                }
                batch_add_edges(single_del_edges_vec);
            }
            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com_set.begin(), max_inf_com_set.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                                           batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                single_del_edges_vec = single_del_edges(max_delete_node);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        single_del_edges_vec.begin(), single_del_edges_vec.end());
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "improved cons double avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/improve_cons_double_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::cut_cons_double_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, const int k, const int r, const double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    map<int, int> add_vec_loc;
    vector<int> cut_vec;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double diff_inf_value;
    double max_diff_inf_value;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int max_delete_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() <= k) {
        }
        else {
            Q.push(make_pair(res[i], avg_vec));
            if (output_res.size() < r) {
                output_res.push_back(make_pair(res[i], avg_vec));
                if (avg_vec < top_r_weight) {
                    top_r_weight = avg_vec;
                    top_r_index = output_res.size() - 1;
                }
            }
            else {
                if (avg_vec > top_r_weight) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(res[i], avg_vec));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (top_r_weight > output_res[j].second) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
            }
        }
    }
    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "cut cons double avg top-r begin" << endl;

    while (Q.size() != 0) {
        double curr_inf_value;
        vector<int> curr_res = Q.front().first;
        add_vec.assign(curr_res.begin(), curr_res.end());
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        while (add_vec.size() >= k + 1) {
            curr_inf_value = avg_vector(add_vec);
            /// store the location of element
            add_vec_loc.clear();
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                add_vec_loc[add_vec[i]] = i;
            }

            /// compute cut vertex
            bool *cut_vertex_flag = new bool [graph_size];
            for (unsigned int i = 0; i < graph_size; i++) {
                cut_vertex_flag[i] = false;
            }
            cut_vec = compute_cut_vertex();
            for (unsigned int i = 0; i < cut_vec.size(); i++) {
                cut_vertex_flag[cut_vec[i]] = true;
            }

            max_diff_inf_value = std::numeric_limits<double>::min();
            max_inf_com.clear();
            for (unsigned int i = 0; i < add_vec.size(); i++) {
                delete_node = add_vec[i];
                if (cut_vertex_flag[delete_node] == false) {
                    bool flag = true;
                    for (unsigned int j = 0; j < AdjList[delete_node].size(); j++) {
                        if (AdjList[AdjList[delete_node][j]].size() <= k) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        if (weight[delete_node] > curr_inf_value) {
                        }
                        else {
                            int del_elem_loc = add_vec_loc[delete_node];
                            add_vec_loc[add_vec.back()] = del_elem_loc;
                            swap(add_vec[del_elem_loc], add_vec.back());
                            add_vec.pop_back();
                            com_max_inf_value = avg_vector(add_vec);
                            diff_inf_value = com_max_inf_value - curr_inf_value;
                            if (diff_inf_value > max_diff_inf_value) {
                                max_inf_value = com_max_inf_value;
                                max_delete_node = delete_node;
                                max_diff_inf_value = diff_inf_value;
                                max_inf_com = add_vec;
                            }
                            add_vec.push_back(delete_node);
                            add_vec_loc[delete_node] = add_vec.size() - 1;
                        }
                    }
                    else {
                        single_del_edges_vec.clear();
                        single_del_edges_vec = single_del_edges(delete_node);
                        vector<vector<int>> divide_com = compute_connected_k_core(k);
                        vector<int> local_max_inf_com;
                        com_max_inf_value = std::numeric_limits<double>::min();
                        for (unsigned int j = 0; j < divide_com.size(); j++) {
                            avg_vec = avg_vector(divide_com[j]);
                            if (avg_vec > com_max_inf_value) {
                                com_max_inf_value = avg_vec;
                                local_max_inf_com = divide_com[j];
                            }
                        }
                        diff_inf_value = com_max_inf_value - curr_inf_value;
                        if (diff_inf_value > max_diff_inf_value) {
                            max_inf_value = com_max_inf_value;
                            max_delete_node = delete_node;
                            max_diff_inf_value = diff_inf_value;
                            max_inf_com = local_max_inf_com;
                        }
                        batch_add_edges(single_del_edges_vec);
                    }
                }
                else {
                    single_del_edges_vec.clear();
                    single_del_edges_vec = single_del_edges(delete_node);
                    vector<vector<int>> divide_com = compute_connected_k_core(k);
                    vector<int> local_max_inf_com;
                    com_max_inf_value = std::numeric_limits<double>::min();
                    for (unsigned int j = 0; j < divide_com.size(); j++) {
                        avg_vec = avg_vector(divide_com[j]);
                        if (avg_vec > com_max_inf_value) {
                            com_max_inf_value = avg_vec;
                            local_max_inf_com = divide_com[j];
                        }
                    }
                    diff_inf_value = com_max_inf_value - curr_inf_value;
                    if (diff_inf_value > max_diff_inf_value) {
                        max_inf_value = com_max_inf_value;
                        max_delete_node = delete_node;
                        max_diff_inf_value = diff_inf_value;
                        max_inf_com = local_max_inf_com;
                    }
                    batch_add_edges(single_del_edges_vec);
                }
            }

            delete [] cut_vertex_flag;

            if (max_inf_com.size() >= k + 1) {
                /// delete edges for next iteration
                set<int> add_local_set(add_vec.begin(), add_vec.end());
                set<int> max_inf_com_set(max_inf_com.begin(), max_inf_com.end());
                set<int> del_nodes_local_set;
                del_nodes_local_vec.clear();
                set_difference(add_local_set.begin(), add_local_set.end(),
                               max_inf_com.begin(), max_inf_com.end(),
                               inserter(del_nodes_local_set, del_nodes_local_set.begin()));
                del_nodes_local_vec.assign(del_nodes_local_set.begin(), del_nodes_local_set.end());
                batch_del_edges_local_vec = batch_del_edges(del_nodes_local_vec);
                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                        batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                /// delete edges for next iteration
                if (output_res.size() < r) {
                    output_res.push_back(make_pair(max_inf_com, max_inf_value));
                    if (max_inf_value < top_r_weight) {
                        top_r_weight = max_inf_value;
                        top_r_index = output_res.size() - 1;
                    }
                }
                else {
                    if (max_inf_value > top_r_weight) {
                        output_res.erase(output_res.begin() + top_r_index);
                        output_res.push_back(make_pair(max_inf_com, max_inf_value));
                        top_r_weight = std::numeric_limits<double>::max();
                        for (unsigned int j = 0; j < output_res.size(); j++) {
                            if (top_r_weight > output_res[j].second) {
                                top_r_weight = output_res[j].second;
                                top_r_index = j;
                            }
                        }
                    }
                }
            }
            add_vec = max_inf_com;
        }
        batch_add_edges(batch_del_edges_vec);
        Q.pop();
    }

    cout << "cut cons double avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/cut_cons_double_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[r-1].first.size(); i++) {
            outfile << output_res[r-1].first[i] << " ";
        }
        outfile << endl;
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[output_res.size() - 1].first.size(); i++) {
            outfile << output_res[output_res.size() - 1].first[i] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void TopRComm::improved_cons_exact_avg_global_topr(const string & graphname, const string & adjpath,
        const string & corepath, const string & pagepath, int k, int r, double p) {
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<vector<int>> candidate;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    queue<pair<vector<int>, double>> Q;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_cons_max_k_core(degree_c, p);
    coreEndTime = clock();
    cout << "compute_cons_max_k_core" << endl;

    res = compute_connected_components();
    cout << "compute connected components" << endl;

    int num_nodes = 0;
    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() <= k) {
        }
        else {
            for (unsigned int j = 0; j < res[i].size(); j++) {
                k_core_vec.push_back(res[i][j]);
                num_nodes = num_nodes + 1;
            }
        }
    }
    cout << "num of nodes: " << num_nodes << endl;

    cout << "exist weight and k core vector" << endl;

    cout << "init: " << output_res.size() << endl;

    cout << "improved cons exact avg top-r begin" << endl;

    candidate = getSubsets(k_core_vec, k);
    for (unsigned int i = 0; i < candidate.size(); i++) {
//        cout << "i: " << i << " total: " << candidate.size() << endl;
        add_vec = candidate[i];
        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        vector<vector<int>> divide_com = compute_connected_k_core(k);
        if (divide_com.size() == 1) {
            output_res.push_back(make_pair(add_vec, avg_vector(add_vec)));
        }
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "improved cons exact avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " p: " << p << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    filename = "O3result/improve_cons_exact_avg_topr_" + graphname + "_" + k_name + r_name + p_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[r-1].first.size(); i++) {
            outfile << output_res[r-1].first[i] << " ";
        }
        outfile << endl;
    }
    else if (output_res.size() != 0) {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
        }
        for (unsigned int i = 0; i < output_res[output_res.size() - 1].first.size(); i++) {
            outfile << output_res[output_res.size() - 1].first[i] << " ";
        }
        outfile << endl;
    }
    outfile.close();
//    vector<int> vec = {1,2,3,4,5,6,7,8,9,10,11,12};
//    vector<vector<int>> res = getSubsets(vec, 3);
//    for (unsigned int i = 0; i < res.size(); i++) {
//        for (unsigned int j = 0; j < res[i].size(); j++) {
//            cout << res[i][j] << " ";
//        }
//        cout << endl;
//    }
}

void TopRComm::climb_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> candidate_res;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<int> local_Q;
    vector<pair<int, double>> local_candidate;
    vector<int> C;
    double top_r_weight;
    int top_r_index;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    batch_del_edges_vec.clear();
    top_r_index = 0;
    top_r_weight = std::numeric_limits<double>::min();
    for (unsigned int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(i);
            in_vec[i] = true;
            while (local_Q.size() != 0 && local_candidate.size()<=s) {
                u = local_Q.front();
                local_candidate.push_back(make_pair(u, weight[u]));
                if (local_Q.size() + local_candidate.size() < s) {
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(v);
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
            sort(local_candidate.begin(), local_candidate.end(),
                 [](const pair<int, double> & left, const pair<int, double> & right){
                     return left.second > right.second;
                 });
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                in_vec[local_candidate[j].first] = false;
            }
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(local_candidate[j].first);
                    core_exist[local_candidate[j].first] = true;
                    if (C.size() >= k + 1 && avg_vector(C) > top_r_weight) {
                        if (determined_k_core(C, core_exist, k)) {
//                            batch_del_edges_local_vec.clear();
//                            batch_del_edges_local_vec = batch_del_edges(C);
//                            batch_del_edges_vec.insert(batch_del_edges_vec.end(),
//                                    batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                            if (output_res.size() >= r) {
                                output_res.erase(output_res.begin() + top_r_index);
                                output_res.push_back(make_pair(C, avg_vector(C)));
                                top_r_weight = std::numeric_limits<double>::max();
                                for (unsigned int m = 0; m < output_res.size(); m++) {
                                    if (output_res[m].second < top_r_weight) {
                                        top_r_weight = output_res[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                            else if (output_res.size() == r - 1) {
                                output_res.push_back(make_pair(C, avg_vector(C)));
                                top_r_weight = std::numeric_limits<double>::max();
                                for (unsigned int m = 0; m < output_res.size(); m++) {
                                    if (output_res[m].second < top_r_weight) {
                                        top_r_weight = output_res[m].second;
                                        top_r_index = m;
                                    }
                                }
                            }
                            else {
                                output_res.push_back(make_pair(C, avg_vector(C)));
                            }
                            for (unsigned int m = 0; m < C.size(); m++) {
                                core_exist[C[m]] = false;
                            }
                            break;
                        }
                    }
                }
                else if (C.size() == s) {
                    for (unsigned int m = 0; m < C.size(); m++) {
                        core_exist[C[m]] = false;
                    }
                }
            }
        }
    }
//    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "climb size avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/climb_size_avg_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << weight[output_res[i].first[j]] << " ";
//            }
//            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::degree_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> candidate_res;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<pair<int, int>> local_Q;
    vector<int> local_vec;
    vector<tuple<int, int, double>> local_candidate;  // (node_id, hierarchy, weight)
    vector<int> C;
    double max_inf_value;
    double com_max_inf_value;
    int delete_node;
    int curr_node;
    double avg_vec;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    int hier;
    int curr_hier;
    vector<pair<int, double>> new_hier;
    batch_del_edges_vec.clear();
    for (unsigned int i = 0; i < graph_size; i++) {
        curr_hier = 0;
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(make_pair(i, 0));
            in_vec[i] = true;
            while (local_Q.size() != 0) {
                u = local_Q.front().first;
                hier = local_Q.front().second;
                if (hier > curr_hier) {
                    if (local_candidate.size() < s) {
                        curr_hier = hier;
                        local_candidate.push_back(make_tuple(u, hier, weight[u]));
                        for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                            v = AdjList[u][j];
                            if (AdjList[v].size() >= k) {
                                if (in_vec[v] == false) {
                                    local_Q.push(make_pair(v, hier + 1));
                                    in_vec[u] = true;
                                }
                            }
                        }
                    }
                }
                else {
                    local_candidate.push_back(make_tuple(u, hier, weight[u]));
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(make_pair(v, hier + 1));
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
            sort(local_candidate.begin(), local_candidate.end(),
                 [](const tuple<int, int, double> & left, const tuple<int, int, double> & right){
                    if (get<1>(left) == get<1>(right)) {
                        return get<2>(left) > get<2>(right);
                    }
                    return get<1>(left) < get<1>(right);
                 });
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                in_vec[get<0>(local_candidate[j])] = false;
            }
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(get<0>(local_candidate[j]));
                    core_exist[get<0>(local_candidate[j])] = true;
                    if (C.size() >= k + 1 && determined_k_core(C, core_exist, k)) {
                        batch_del_edges_local_vec.clear();
                        batch_del_edges_local_vec = batch_del_edges(C);
                        batch_del_edges_vec.insert(batch_del_edges_vec.end(),
                                batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                        output_res.push_back(make_pair(C, avg_vector(C)));
                        for (unsigned int m = 0; m < C.size(); m++) {
                            core_exist[C[m]] = false;
                        }
                        break;
                    }
                }
                else if (C.size() == s) {
                    for (unsigned int m = 0; m < C.size(); m++) {
                        core_exist[C[m]] = false;
                    }
                }
            }
        }
    }
    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "degree size avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/degree_size_avg_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << output_res[i].first[j] << " ";
            }
            outfile << endl;
            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
                outfile << weight[output_res[i].first[j]] << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::casual_size_avg_global_topr(const string & graphname, const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const int r, const int s) {
    if (s < k + 1) {
        return;
    }
    clock_t startTime, endTime;
    clock_t coreStartTime, coreEndTime;
    clock_t localStartTime, localEndTime;
    vector<vector<int>> res;
    vector<int> max_inf_com;
    vector<pair<vector<int>, double>> output_res;
    vector<pair<vector<int>, double>> L;
    vector<pair<vector<int>, double>> candidate_res;
    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<int> del_nodes_local_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_local_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    queue<int> local_Q;
    vector<pair<int, double>> local_candidate;
    vector<pair<vector<int>, double>> k_candidate;
    vector<int> C;
    double top_r_weight;
    int top_r_index;
    double k_top_weight;
    int k_top_index;
    degree_c = k;
    size_c = r;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    startTime = clock();
    coreStartTime = clock();
    compute_max_k_core(degree_c);
    coreEndTime = clock();
    cout << "compute_max_k_core" << endl;

    bool *in_vec = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        in_vec[i] = false;
    }
    bool *core_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        core_exist[i] = false;
    }
    localStartTime = clock();
    int u, v;
    batch_del_edges_vec.clear();
    top_r_index = 0;
    top_r_weight = std::numeric_limits<double>::min();
    for (unsigned int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            local_candidate.clear();
            C.clear();
            while (!local_Q.empty()) {
                local_Q.pop();
            }
            local_Q.push(i);
            in_vec[i] = true;
            while (local_Q.size() != 0) {
                u = local_Q.front();
                local_candidate.push_back(make_pair(u, weight[u]));
                if (local_Q.size() + local_candidate.size() < s) {
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        if (AdjList[v].size() >= k) {
                            if (in_vec[v] == false) {
                                local_Q.push(v);
                                in_vec[v] = true;
                            }
                        }
                    }
                }
                local_Q.pop();
            }
//            srand ( unsigned ( time(0) ) );
//            random_shuffle(local_candidate.begin(), local_candidate.end());
            k_candidate.clear();
            for (unsigned int j = 0; j < local_candidate.size(); j++) {
                if (C.size() < s) {
                    C.push_back(local_candidate[j].first);
                    core_exist[local_candidate[j].first] = true;
                    if (C.size() >= k + 1 && avg_vector(C) > top_r_weight) {
                        if (determined_k_core(C, core_exist, k)) {
                            k_candidate.push_back(make_pair(C, avg_vector(C)));
                        }
                    }
                }
                else if (C.size() == s) {
                    for (unsigned int m = 0; m < C.size(); m++) {
                        core_exist[C[m]] = false;
                    }
                }
                in_vec[local_candidate[j].first] = false;
            }
            if (k_candidate.size() != 0) {
                k_top_index = 0;
                k_top_weight = std::numeric_limits<double>::min();
                for (unsigned int j = 0; j < k_candidate.size(); j++) {
                    if (k_candidate[j].second > k_top_weight) {
                        k_top_weight = k_candidate[j].second;
                        k_top_index = j;
                    }
                }
//                batch_del_edges_local_vec.clear();
//                batch_del_edges_local_vec = batch_del_edges(k_candidate[k_top_index].first);
//                batch_del_edges_vec.insert(batch_del_edges_vec.end(),
//                        batch_del_edges_local_vec.begin(), batch_del_edges_local_vec.end());
                if (output_res.size() >= r) {
                    output_res.erase(output_res.begin() + top_r_index);
                    output_res.push_back(make_pair(k_candidate[k_top_index].first,
                            k_candidate[k_top_index].second));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (output_res[j].second < top_r_weight) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
                else if (output_res.size() == r - 1) {
                    output_res.push_back(make_pair(k_candidate[k_top_index].first,
                            k_candidate[k_top_index].second));
                    top_r_weight = std::numeric_limits<double>::max();
                    for (unsigned int j = 0; j < output_res.size(); j++) {
                        if (output_res[j].second < top_r_weight) {
                            top_r_weight = output_res[j].second;
                            top_r_index = j;
                        }
                    }
                }
                else {
                    output_res.push_back(make_pair(k_candidate[k_top_index].first,
                            k_candidate[k_top_index].second));
                }
            }
        }
    }
//    batch_add_edges(batch_del_edges_vec);
    localEndTime = clock();
    delete [] in_vec;
    delete [] core_exist;
    cout << "Local Search Time : " << (double)(localEndTime - localStartTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "output_res size: " << output_res.size() << endl;

    cout << "casual size avg top-r end" << endl;

    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    cout << "output res size: " << output_res.size() << endl;

    endTime = clock();
    cout << "k: " << k << " r: " << r << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << "s" << endl;

    // write result to file
    char timename[128];
    struct tm *newtime;
    time_t now = time(nullptr);
    newtime = localtime(&now);
    strftime(timename, 128, "%Y%m%d_%H%M%S.txt", newtime);
    string filename = timename;
    string k_name = "k_" + std::to_string(k) + "_";
    string r_name = "r_" + std::to_string(r) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    filename = "O3result/casual_size_avg_topr_" + graphname + "_" + k_name + r_name + s_name + filename;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime) / CLOCKS_PER_SEC << endl;
    if (output_res.size() >= r) {
        for (unsigned int i = 0; i < r; i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
        }
    }
    else {
        for (unsigned int i = 0; i < output_res.size(); i++) {
            outfile << i + 1 << " " << setprecision(8) << output_res[i].second
                    << " " << output_res[i].first.size() << endl;
//            for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//                outfile << output_res[i].first[j] << " ";
//            }
//            outfile << endl;
        }
    }
    outfile.close();
}

void TopRComm::split(const string& src, const string& separator, vector<string>& dest)
{
    string str = src;
    string substring;
    string::size_type start = 0, index;

    do
    {
        index = str.find_first_of(separator,start);
        if (index != string::npos)
        {
            substring = str.substr(start,index-start);
            dest.push_back(substring);
            start = str.find_first_not_of(separator,index);
            if (start == string::npos) return;
        }
    }while(index != string::npos);

    //the last token
    substring = str.substr(start);
    dest.push_back(substring);
}

void TopRComm::add_edge(const int start, const int end) {
    if (start < end) {
        AdjList[start].push_back(end);
        AdjList_Loc[start][end] = AdjList[start].size() - 1;
        AdjList[end].push_back(start);
        AdjList_Loc[end][start] = AdjList[end].size() - 1;
    }
}

// one direction which means for edge (u,v), we only remove u, v but remain v, u
void TopRComm::remove_edge(const int start, const int end) {
    int idx = AdjList_Loc[start][end];
    AdjList_Loc[start][AdjList[start].back()] = idx;
    swap(AdjList[start][idx], AdjList[start].back());
    AdjList[start].pop_back();
}

void TopRComm::compute_max_k_core(const int k) {
    for (int i = 0; i < graph_size; i++) {
        if (core[i] < k) {
            for (unsigned int j = 0; j < AdjList[i].size(); j++) {
                remove_edge(AdjList[i][j], i);
            }
            AdjList[i].clear();
            c_size = c_size - 1;
        }
    }
}

void TopRComm::compute_cons_max_k_core(const int k, const double p) {
    int u, v;
    queue<int> Q;
    /// initial threshold
    double max_weight = 0.0;
    for (int i = 0; i < graph_size; i++) {
        if (max_weight < weight[i]) {
            max_weight = weight[i];
        }
    }
    min_thres = max_weight * p;
    for (int i = 0; i < graph_size; i++) {
        if (core[i] < k) {
            for (unsigned int j = 0; j < AdjList[i].size(); j++) {
                remove_edge(AdjList[i][j], i);
            }
            AdjList[i].clear();
            c_size = c_size - 1;
        }
    }
    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() != 0) {
            if (AdjList[i].size() < k || weight[i] < min_thres) {
                Q.push(i);
                while (Q.size() != 0) {
                    u = Q.front();
                    if (AdjList[u].size() < k || weight[u] < min_thres) {
                        for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                            v = AdjList[u][j];
                            remove_edge(v, u);
                            Q.push(v);
                        }
                        AdjList[u].clear();
                        c_size = c_size - 1;
                    }
                    Q.pop();
                }
            }
        }
    }
}

vector<int> TopRComm::compute_cut_vertex() {
    vector<int> res;
    int inde = 0;
    int u = 0;
    int v = 0;
    int father = 0;
    int neigh;
    int *num = new int [graph_size];
    int *low = new int [graph_size];
    int *children = new int [graph_size];
    int *flag = new int [graph_size];
    int *neighbor = new int [graph_size];
    bool *vis = new bool [graph_size];
    stack<int> S;
    stack<int> S_father;

    // init
    for (unsigned int i = 0; i < graph_size; i++) {
        children[i] = 0;
        vis[i] = false;
        flag[i] = false;
        neighbor[i] = 0;
    }

    for (u = 0; u < graph_size; u++) {
        if (!vis[u]) {
            S.push(u);
            S_father.push(u);
            while(S.size() != 0) {
                u = S.top();
                father = S_father.top();
                if (vis[u] == false) {
                    vis[u] = true;
                    low[u] = num[u] = ++inde;
                }
                for (neigh = neighbor[u]; neigh < AdjList[u].size(); neigh++) {
                    v = AdjList[u][neigh];
                    neighbor[u]++;
                    if (!vis[v]) {
                        children[u]++;
                        S.push(v);
                        S_father.push(u);
                        break;
                    }
                    else if (v != father) {
                        low[u] = min(low[u], num[v]);
                    }
                }
                if (neigh == AdjList[u].size()) {
                    if (father == u && children[u] >= 2 && !flag[u]) {
                        flag[u] = true;
                    }
                    if (S.size() > 1) {
                        v = u;
                        u = father;
                        S.pop();
                        S_father.pop();
                        father = S_father.top();
                        low[u] = min(low[u], low[v]);
                        if (father != u && low[v] >= num[u] && !flag[u]) {
                            flag[u] = true;
                        }
                    }
                    else {
                        S.pop();
                        S_father.pop();
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < graph_size; i++) {
        if (flag[i] == true) {
            res.push_back(i);
        }
    }

    delete [] num;
    delete [] low;
    delete [] children;
    delete [] flag;
    delete [] neighbor;
    delete [] vis;

    return res;
}

pair<vector<int>, int> TopRComm::compute_minimal_connected_components() {
    bool *visited = new bool [graph_size];
    int u;
    vector<int> res;
    queue<int> q;
    double min_weight = std::numeric_limits<double>::max();
    int min_idx = 0;
    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
        if (weight[v] < min_weight && AdjList[v].size() != 0) {
            min_weight = weight[v];
            min_idx = v;
        }
    }

    q.push(min_idx);
    visited[min_idx] = true;
    while(q.size() != 0) {
        u = q.front();
        res.push_back(u);
        q.pop();
        for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
            if (visited[*iter] == false && AdjList[*iter].size() != 0) {
                q.push(*iter);
                visited[*iter] = true;
            }
        }
    }

    delete [] visited;
    return make_pair(res, min_idx);
}

pair<vector<int>, int> TopRComm::compute_maximal_connected_components() {
    bool *visited = new bool [graph_size];
    int u;
    vector<int> res;
    queue<int> q;
    double max_weight = std::numeric_limits<double>::min();
    int max_idx = 0;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
        if (weight[v] > max_weight && AdjList[v].size() != 0) {
            max_weight = weight[v];
            max_idx = v;
        }
    }

    q.push(max_idx);
    visited[max_idx] = true;
    while (q.size() != 0) {
        u = q.front();
        res.push_back(u);
        q.pop();
        for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
            if (visited[*iter] == false && AdjList[*iter].size() != 0) {
                q.push(*iter);
                visited[*iter] = true;
            }
        }
    }

    delete [] visited;
    return make_pair(res, max_idx);
}

vector<vector<int>> TopRComm::compute_connected_components() {
    bool *visited = new bool [graph_size];
    int u;
    vector<int> id_vec;
    queue<int> q;
    vector<vector<int>> res;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    for (int v = 0; v < graph_size; v++) {
        id_vec.clear();
        if (visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (unsigned int i = 0; i < AdjList[u].size(); i++) {
                    if (visited[AdjList[u][i]] == false && AdjList[AdjList[u][i]].size() != 0) {
                        q.push(AdjList[u][i]);
                        visited[AdjList[u][i]] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

//    /// check connected components # vertex
//    int sum = 0;
//    for(auto & re : res) {
//        sum += re.size();
//    }
//    assert(sum == c_size);

    delete [] visited;
    return res;
}

// if remove vertex, the component does not satisfy k-core and the component is divided into several components
vector<vector<int>> TopRComm::compute_connected_k_core(const int k) {
    vector<vector<int>> res;
    int u;
    queue<int> q;
    bool *visited = new bool [graph_size];
    vector<int> id_vec;
    vector<pair<int, int>> batch_del_edges_vec;

    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        if (visited[i] == false && AdjList[i].size() < k) {
            q.push(i);
            visited[i] = true;
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                    if (u < AdjList[u][j]) {
                        batch_del_edges_vec.push_back(make_pair(u, AdjList[u][j]));
                    }
                    else {
                        batch_del_edges_vec.push_back(make_pair(AdjList[u][j], u));
                    }
                    remove_edge(AdjList[u][j], u);
                    if (visited[AdjList[u][j]] == false && AdjList[AdjList[u][j]].size() < k) {
                        q.push(AdjList[u][j]);
                        visited[AdjList[u][j]] = true;
                    }
                }
                AdjList[u].clear();
            }
        }
    }

    assert(q.size() == 0);

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    for (int v = 0; v < graph_size; v++) {
        id_vec.clear();
        if (visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (unsigned int i = 0; i < AdjList[u].size(); i++) {
                    if (visited[AdjList[u][i]] == false && AdjList[AdjList[u][i]].size() != 0) {
                        q.push(AdjList[u][i]);
                        visited[AdjList[u][i]] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

    batch_add_edges(batch_del_edges_vec);
    delete [] visited;
    return res;
}

vector<pair<int, int>> TopRComm::batch_del_edges(const vector<int> & nodes){
    int u;
    vector<pair<int, int>> res;
    for(unsigned int i = 0; i < nodes.size(); i++) {
        u = nodes[i];
        for (unsigned int j = 0; j < AdjList[u].size(); j++) {
            if (u < AdjList[u][j]) {
                res.push_back(make_pair(u, AdjList[u][j]));
            }
            else {
                res.push_back(make_pair(AdjList[u][j], u));
            }
            remove_edge(AdjList[u][j], u);
        }
        AdjList[u].clear();
    }
    return res;
}

vector<pair<int, int>> TopRComm::single_del_edges(const int node) {
    vector<pair<int, int>> res;
    for (unsigned int i = 0; i < AdjList[node].size(); i++) {
        if (node < AdjList[node][i]) {
            res.push_back(make_pair(node, AdjList[node][i]));
        }
        else {
            res.push_back(make_pair(AdjList[node][i], node));
        }
        remove_edge(AdjList[node][i], node);
    }
    AdjList[node].clear();
    return res;
}

void TopRComm::batch_add_edges(const vector<pair<int, int>> & edges){
    for(unsigned int i = 0; i < edges.size(); i++) {
        if (edges[i].first < edges[i].second) {
            add_edge(edges[i].first, edges[i].second);
        }
    }
}

bool TopRComm::is_connected_k_core(const vector<int> & vertex_set, const int k) {
    bool *visited = new bool [graph_size];
    int count = 0;
    int u;
    int k_neigh = 0;
    bool is_in_vec = false;
    queue<int> q;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    for (int v = 0; v < graph_size; v++) {
        if (visited[v] == false && AdjList[v].size() != 0 && is_element_in_vector(vertex_set, v)) {
            q.push(v);
            visited[v] = true;
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                k_neigh = 0;
                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
                    is_in_vec = is_element_in_vector(vertex_set, v);
                    if (is_in_vec && AdjList[*iter].size() != 0) {
                        if (visited[*iter] == false) {
                            q.push(*iter);
                            visited[*iter] = true;
                        }
                        k_neigh += 1;
                    }
                }
                if (k_neigh < k) {
                    return false;
                }
            }
            count += 1;
        }
        if (count > 1) {
            delete [] visited;
            return false;
        }
    }
    return true;
}

bool TopRComm::is_connected_k_core_graph(const int k) {
    bool *visited = new bool [graph_size];
    int count = 0;
    int u;
    int k_neigh = 0;
    queue<int> q;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    for (int v = 0; v < graph_size; v++) {
        if (visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                k_neigh = 0;
                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
                    if (AdjList[*iter].size() != 0) {
                        if (visited[*iter] == false) {
                            q.push(*iter);
                            visited[*iter] = true;
                        }
                        k_neigh += 1;
                    }
                }
                if (k_neigh < k) {
                    delete [] visited;
                    return false;
                }
            }
            count += 1;
        }
        if (count > 1) {
            delete [] visited;
            return false;
        }
    }
    delete [] visited;
    if (count == 0) {
        return false;
    }
    return true;
}

int TopRComm::greedy_unsatisfied_node(const int k) {
    bool *visited = new bool [graph_size];
    int count = 0;
    int u;
    int k_neigh = 0;
    queue<int> q;
    vector<pair<int, double>> uns_node_vec;
    vector<pair<int, double>> comp_max_node_vec;
    double max_node_weight;
    int max_node_id;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    uns_node_vec.clear();
    comp_max_node_vec.clear();
    for (int v = 0; v < graph_size; v++) {
        if (visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            max_node_weight = weight[v];
            max_node_id = v;
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                k_neigh = 0;
                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
                    if (AdjList[*iter].size() != 0) {
                        if (visited[*iter] == false) {
                            q.push(*iter);
                            visited[*iter] = true;
                        }
                        k_neigh += 1;
                    }
                }
                if (k_neigh < k) {
                    uns_node_vec.push_back(make_pair(u, weight[u]));
                }
                if (weight[u] > max_node_weight) {
                    max_node_weight = weight[u];
                    max_node_id = u;
                }
            }
            count += 1;
            comp_max_node_vec.push_back(make_pair(max_node_id, max_node_weight));
        }
    }

    delete [] visited;
    double max_uns_node_weight, max_comp_node_weight;
    int max_uns_node_id, max_comp_node_id;
    if (comp_max_node_vec.size() > 1) {
        if (uns_node_vec.size() == 0) {
            max_comp_node_weight = comp_max_node_vec[0].second;
            max_comp_node_id = comp_max_node_vec[0].first;
            for (unsigned int i = 0; i < comp_max_node_vec.size(); i++) {
                if (comp_max_node_vec[i].second > max_comp_node_weight) {
                    max_comp_node_weight = comp_max_node_vec[i].second;
                    max_comp_node_id = comp_max_node_vec[i].first;
                }
            }
            return max_comp_node_id;
        }
        else {
            max_uns_node_weight = uns_node_vec[0].second;
            max_uns_node_id = uns_node_vec[0].first;
            for (unsigned int i = 0; i < uns_node_vec.size(); i++) {
                if (uns_node_vec[i].second > max_uns_node_weight) {
                    max_uns_node_weight = uns_node_vec[i].second;
                    max_uns_node_id = uns_node_vec[i].first;
                }
            }
            max_comp_node_weight = comp_max_node_vec[0].second;
            max_comp_node_id = comp_max_node_vec[0].first;
            for (unsigned int i = 0; i < comp_max_node_vec.size(); i++) {
                if (comp_max_node_vec[i].second > max_comp_node_weight) {
                    max_comp_node_weight = comp_max_node_vec[i].second;
                    max_comp_node_id = comp_max_node_vec[i].first;
                }
            }
            if (max_comp_node_weight > max_uns_node_weight) {
                return max_comp_node_id;
            }
            else {
                return max_uns_node_id;
            }
        }
    }
    else {
        if (uns_node_vec.size() == 0) {
            return -1;
        }
        else {
            max_uns_node_weight = uns_node_vec[0].second;
            max_uns_node_id = uns_node_vec[0].first;
            for (unsigned int i = 0; i < uns_node_vec.size(); i++) {
                if (uns_node_vec[i].second > max_uns_node_weight) {
                    max_uns_node_weight = uns_node_vec[i].second;
                    max_uns_node_id = uns_node_vec[i].first;
                }
            }
            return max_uns_node_id;
        }
    }
}

int TopRComm::random_unsatisfied_node(const vector<int> & vec, const int k) {
    bool *visited = new bool [graph_size];
    int count = 0;
    int u;
    int k_neigh = 0;
    queue<int> q;
    vector<int> uns_node_idx_vec;
    vector<int> comp_max_node_idx_vec;
    map<int, int> vec_loc;
    int max_node_idx;

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec_loc[vec[i]] = i;
    }

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    uns_node_idx_vec.clear();
    comp_max_node_idx_vec.clear();
    for (int v = 0; v < graph_size; v++) {
        if (visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            max_node_idx = vec_loc[v];
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                k_neigh = 0;
                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
                    if (AdjList[*iter].size() != 0) {
                        if (visited[*iter] == false) {
                            q.push(*iter);
                            visited[*iter] = true;
                        }
                        k_neigh += 1;
                    }
                }
                if (k_neigh < k) {
                    uns_node_idx_vec.push_back(vec_loc[u]);
                }
                if (vec_loc[u] > max_node_idx) {
                    max_node_idx = vec_loc[u];
                }
            }
            count += 1;
            comp_max_node_idx_vec.push_back(max_node_idx);
        }
    }

    delete [] visited;
    int max_uns_node_idx, max_comp_node_idx;
    if (comp_max_node_idx_vec.size() > 1) {
        if (uns_node_idx_vec.size() == 0) {
            max_comp_node_idx = comp_max_node_idx_vec[0];
            for (unsigned int i = 0; i < comp_max_node_idx_vec.size(); i++) {
                if (comp_max_node_idx_vec[i] > max_comp_node_idx) {
                    max_comp_node_idx = comp_max_node_idx_vec[i];
                }
            }
            return max_comp_node_idx;
        }
        else {
            max_uns_node_idx = uns_node_idx_vec[0];
            for (unsigned int i = 0; i < uns_node_idx_vec.size(); i++) {
                if (uns_node_idx_vec[i] > max_uns_node_idx) {
                    max_uns_node_idx = uns_node_idx_vec[i];
                }
            }
            max_comp_node_idx = comp_max_node_idx_vec[0];
            for (unsigned int i = 0; i < comp_max_node_idx_vec.size(); i++) {
                if (comp_max_node_idx_vec[i] > max_comp_node_idx) {
                    max_comp_node_idx = comp_max_node_idx_vec[i];
                }
            }
            if (max_comp_node_idx > max_uns_node_idx) {
                return max_comp_node_idx;
            }
            else {
                return max_uns_node_idx;
            }
        }
    }
    else {
        if (uns_node_idx_vec.size() == 0) {
            return -1;
        }
        else {
            max_uns_node_idx = uns_node_idx_vec[0];
            for (unsigned int i = 0; i < uns_node_idx_vec.size(); i++) {
                if (uns_node_idx_vec[i] > max_uns_node_idx) {
                    max_uns_node_idx = uns_node_idx_vec[i];
                }
            }
            return max_uns_node_idx;
        }
    }
}

vector<vector<int>> TopRComm::compute_part_connected_components(const vector<int> & vertex_set) {
    bool *visited = new bool [graph_size];
    int u;
    queue<int> q;
    vector<int> id_vec;
    vector<vector<int>> res;

    for (int v = 0; v < graph_size; v++) {
        visited[v] = false;
    }

    for (int v = 0; v < graph_size; v++) {
        id_vec.clear();
        if (is_element_in_vector(vertex_set, v) && visited[v] == false && AdjList[v].size() != 0) {
            q.push(v);
            visited[v] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
                    if (is_element_in_vector(vertex_set, v) && visited[*iter] == false && AdjList[*iter].size() != 0) {
                        q.push(*iter);
                        visited[*iter] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

    delete [] visited;
    return res;
}

void TopRComm::dfs_deletion(const int v, const int k) {
    vector<int> neighs = AdjList[v];
    AdjList[v].clear();
    c_size = c_size - 1;
    for (int i = 0; i < neighs.size(); i++) {
        remove_edge(neighs[i], v);
        if (AdjList[neighs[i]].size() == 0) {
            AdjList[neighs[i]].clear();
        }
        else if (AdjList[neighs[i]].size() < k && AdjList[neighs[i]].size() != 0) {
            dfs_deletion(neighs[i], k);
        }
    }
}

bool TopRComm::determined_k_core(const vector<int> & node_vec, const bool *is_exist, const int k) {
    int count = 0;
    int node = 0;

    for (unsigned int i = 0; i < node_vec.size(); i++) {
        count = 0;
        node = node_vec[i];
        for (auto iter = AdjList[node].begin(); iter != AdjList[node].end(); iter++) {
            if (is_exist[*iter] == true) {
                count += 1;
            }
        }
        if (count < k) {
            return false;
        }
    }
    return true;
}

int TopRComm::min_degree(const vector<int> & node_vec) {
    int count = 0;
    int node = 0;
    int min_deg = max_core;

    bool *is_exist = new bool [graph_size];
    for (unsigned int i = 0; i < graph_size; i++) {
        is_exist[i] = false;
    }
    for (unsigned int i = 0; i < node_vec.size(); i++) {
        is_exist[node_vec[i]] = true;
    }
    for (unsigned int i = 0; i < node_vec.size(); i++) {
        count = 0;
        node = node_vec[i];
        for (auto iter = AdjList[node].begin(); iter != AdjList[node].end(); iter++) {
            if (is_exist[*iter] == true) {
                count += 1;
            }
        }
        if (count < min_deg) {
            min_deg = count;
        }
    }
    delete [] is_exist;
    return min_deg;
}

double TopRComm::sum_vector(const vector<int> & vec) {
    double sum = 0.0;
    for (int i : vec) {
        sum += weight[i];
    }
    return sum;
}

double TopRComm::avg_vector(const vector<int> & vec) {
    double sum = 0.0;
    for (int i : vec) {
        sum += weight[i];
    }
    return sum * 1.0 / vec.size();
}

vector<vector<int>> TopRComm::getSubsets(const vector<int> & A, const int k) {
    int n = A.size();
    long long int size = 1 << n;
    vector<vector<int>> subsets;
    for (int i = size - 1; i > 0; --i) {
        vector<int> subset;
        for (int j = n - 1; j >= 0; --j) {
            if ((i >> j) & 1) {
                subset.push_back(A[j]);
            }
        }
        if (subset.size() > k) {
            subsets.emplace_back(subset);
        }
    }
    return subsets;
}

bool TopRComm::is_element_in_vector(const vector<int> & vec, const int element) {
    bool flag = false;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] == element) {
            flag = true;
            break;
        }
    }
    return flag;
}

bool TopRComm::is_keep_avg_component(const vector<int> & vec, const double threshold, const int k) {
    double weight_sum = 0.0;
    double weight_avg;
    vector<double> weight_vec;
    weight_vec.resize(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++) {
        weight_vec[i] = weight[vec[i]];
    }
    sort(weight_vec.begin(), weight_vec.end(), [](const double & left, const double & right){
        return left > right;
    });
    for (unsigned int i = 0; i < k + 1; i++) {
        weight_sum += weight_vec[i];
    }
    weight_avg = weight_sum * 1.0 / (k + 1);
    if (weight_avg > threshold) {
        return true;
    }
    else {
        return false;
    }
}

void TopRComm::combine_inner(const vector<int> & data, int start, int n, int m, int depth,
        vector<int> temp, vector<vector<int>> result) {
    if (depth == m - 1) {
        for (int i = start; i < n - (m - depth - 1); ++i) {
            temp[depth] = data[i];
            result.push_back(temp);
        }
    }
    else {
        for (int i = start; i < n - (m - depth - 1); ++i) {
            temp[depth] = data[i];
            combine_inner(data,i + 1, n, m, depth+1,temp,result);
        }
    }
}
vector<vector<int>> TopRComm::combine(const vector<int> & data, int m) {
    if (m <= 0)
        return{};
    int depth = 0;
    vector<vector<int>> result;
    vector<int> temp(m,0);
    combine_inner(data,0, data.size(), m, depth,temp,result);
    return result;
}

vector<int> TopRComm::peel_vector(vector<int> vec, const int k) {
    vector<list<int>> AdjList_copy;
    vector<pair<int, int>> vec_vertex;
    vector<int> res;
    vector<int> remove_item;
    AdjList_copy.resize(graph_size);

    for (unsigned int i = 0; i < vec.size(); i++) {
        AdjList_copy[vec[i]].assign(AdjList[vec[i]].begin(), AdjList[vec[i]].end());
    }

    for (unsigned int i = 0; i < vec.size(); i++) {
        remove_item.clear();
        for (auto j = AdjList_copy[i].begin(); j != AdjList_copy[i].end(); j++) {
//            cout << *j << endl;
            if (is_element_in_vector(vec, *j)) {
            }
            else {
                remove_item.push_back(*j);
//                AdjList_copy[i].remove(*j);
            }
        }
        for (unsigned int j = 0; j < remove_item.size(); j++) {
            AdjList_copy[i].remove(remove_item[j]);
        }
    }

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec_vertex.push_back(make_pair(vec[i], AdjList_copy[i].size()));
    }

    sort(vec_vertex.begin(), vec_vertex.end(), [](const pair<int, int> & left, const pair<int, int> & right){
        return left.second < right.second;
    });

    for (unsigned int i = 0; i < vec_vertex.size(); i++) {
        sort(vec_vertex.begin() + i, vec_vertex.end(), [](const pair<int, int> & left, const pair<int, int> & right){
            return left.second < right.second;
        });
        if (vec_vertex[i].second >= k) {
            res.push_back(vec_vertex[i].first);
        }
        for (auto j = AdjList_copy[vec_vertex[i].first].begin(); j != AdjList_copy[vec_vertex[i].first].end(); j++) {
            if (AdjList_copy[*j].size() > AdjList_copy[vec_vertex[i].first].size()) {
                AdjList_copy[*j].remove(vec_vertex[i].first);
            }
        }
        AdjList_copy[vec_vertex[i].first].clear();
    }
    return res;
//    int u;
//    queue<int> q;
//    bool visited[graph_size];
//
//    for (int i = 0; i < graph_size; i++) {
//        visited[i] = false;
//    }
//
//    for (int i = 0; i < graph_size; i++) {
//        if (visited[i] == false && AdjList[i].size() < k) {
//            q.push(i);
//            visited[i] = true;
//            while (q.size() != 0) {
//                u = q.front();
//                q.pop();
//                for (auto iter = AdjList[u].begin(); iter != AdjList[u].end(); iter++) {
//                    AdjList[*iter].remove(u);
//                    if (visited[*iter] == false && AdjList[*iter].size() < k) {
//                        q.push(*iter);
//                        visited[*iter] = true;
//                    }
//                }
//                AdjList[u].clear();
//                c_size = c_size - 1;
//            }
//        }
//    }
}

//void TopRComm::naive_avg_global_topr(const string & adjpath, const string & corepath,
//                                     const string & pagepath, const int k, const int r) {
//    clock_t startTime, endTime;
//    vector<vector<int>> res;
//    vector<pair<vector<int>, double>> temp;
//    vector<pair<vector<int>, double>> output_res;
//    vector<int> add_vec; // after delete one node, new vector
//    vector<int> k_core_vec;
//    vector<int> del_nodes_vec;
//    vector<pair<int, int>> batch_del_edges_vec;
//    vector<pair<int, int>> single_del_edges_vec;
//    vector<double> weight_vec;
//    vector<double> exist_weight;
//    double update_weight; // updated weight
//    double upper_bound = -1.0;
//    double lower_bound = -1.0;
//    double bound = -1.0;
//    double eps = 0.1;
//    double avg_vec;
//    degree_c = k;
//    size_c = r;
//    load_data(adjpath, corepath, pagepath);
//    cout << "successfully load data" << endl;
//    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
//    bool is_exist[graph_size];
//    for (int i = 0; i < graph_size; i++) {
//        is_exist[i] = false;
//    }
//    startTime = clock();
//    compute_max_k_core(degree_c);
//    cout << "compute_max_k_core" << endl;
//
//    res = compute_connected_components();
//    cout << "compute connected components" << endl;
//
//    if (res.size() == 0) {
//        return;
//    }
//
//    for (unsigned int i = 0; i < res.size(); i++) {
//        for (unsigned int j = 0; j < res[i].size(); j++) {
//            is_exist[res[i][j]] = true;
//            exist_weight.push_back(weight[res[i][j]]);
//            k_core_vec.push_back(res[i][j]);
//        }
//    }
//
//    sort(exist_weight.begin(), exist_weight.end(), [](double & left, double & right) {
//        return left > right;
//    });
//
//    upper_bound = 0.0;
//    for (unsigned int i = 0; i < (k + 1); i++) {
//        upper_bound += exist_weight[i];
//    }
//    upper_bound = upper_bound * 1.0 / (k + 1);
//    cout << "upper bound: " << upper_bound << endl;
//
//    bound = upper_bound * (1 - eps);
//    for (unsigned int i = 0; i < res.size(); i++) {
//        avg_vec = avg_vector(res[i]);
//        if (res[i].size() <= k) {
//        }
//        else {
//            if (avg_vec >= bound) {
//                output_res.push_back(make_pair(res[i], avg_vec));
//            }
//        }
//    }
//    upper_bound = -1.0;
//
//    cout << "init: " << output_res.size() << endl;
//
//    cout << "naive avg top-r begin" << endl;
//
//    while (true) {
//        for (int i = 0; i < graph_size; i++) {
//            if (is_exist[i] == true) {
//                for (int j = 0; j < output_res.size(); j++) {
//                    if (output_res[j].first.size() <= k + 1) {
//                    }
//                    else {
//                        add_vec.assign(output_res[j].first.begin(), output_res[j].first.end());
//                        for (unsigned int m = 0; m < add_vec.size(); m++) {
//                            if (add_vec[m] == i) {
//                                add_vec.erase(add_vec.begin() + m);
//                                break;
//                            }
//                        }
//                        update_weight = avg_vector(add_vec);
//                        if (add_vec.size() == output_res[j].first.size()) {
//                        }
//                        else if (is_keep_avg_component(add_vec, bound, k)) {
//                            set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
//                            set<int> add_set(add_vec.begin(), add_vec.end());
//                            set<int> del_nodes_set;
//                            del_nodes_vec.clear();
//                            set_difference(k_core_set.begin(), k_core_set.end(),
//                                           add_set.begin(), add_set.end(),
//                                           inserter(del_nodes_set, del_nodes_set.begin()));
//                            del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
//                            batch_del_edges_vec = batch_del_edges(del_nodes_vec);
//                            if (is_connected_k_core_graph(k)) {
//                                output_res.push_back(make_pair(add_vec, update_weight));
//                                if (output_res.size() == r) {
//                                    break;
//                                }
//                            }
//                            else {
//                                vector<vector<int>> divide_com = compute_connected_k_core(k);
//                                for (unsigned int m = 0; m < divide_com.size(); m++) {
//                                    double avg_divide_com = avg_vector(divide_com[m]);
//                                    if (divide_com[m].size() <= k) {
//                                    }
//                                    else if (is_keep_avg_component(divide_com[m], bound, k)) {
//                                        output_res.push_back(make_pair(divide_com[m], avg_divide_com));
//                                        if (output_res.size() == r) {
//                                            break;
//                                        }
//                                    }
//                                }
//                            }
//                            batch_add_edges(batch_del_edges_vec);
//                        }
//                        if (output_res.size() == r) {
//                            break;
//                        }
//                    }
//                }
//            }
//            if (output_res.size() == r) {
//                break;
//            }
//        }
//        cout << "output_res size: " << output_res.size() << " bound value: " << bound << endl;
//        cout << "upper bound value: " << upper_bound << " lower bound value: " << lower_bound << endl;
//        if (output_res.size() == 0) {
//            bound = bound * (1 - eps);
//        }
//        else if (output_res.size() < r && output_res.size() > 0) {
//            upper_bound = bound;
//            if (lower_bound < 0) {
//                bound = (bound * (1 - eps) + bound) / 2;
//            }
//            else {
//                bound = (upper_bound + lower_bound) / 2;
//            }
//        }
//        else if (output_res.size() == r) {
//            lower_bound = bound;
//            if (upper_bound < 0) {
//                bound = (bound / (1 - eps) + bound) / 2;
//            }
//            else {
//                bound = (upper_bound + lower_bound) / 2;
//            }
//        }
//        sort(output_res.begin(), output_res.end(), [](const pair<vector<int>, double> & left,
//                                                      const pair<vector<int>, double> & right){
//            return left.second > right.second;
//        });
//        if (output_res.size() != 0 && fabs(output_res[r - 1].second - bound) < err) {
//            break;
//        }
//        if (output_res.size() != 0) {
//            cout << "top-r value: " << output_res[r - 1].second << endl;
//        }
//        output_res.clear();
//        for (unsigned int i = 0; i < res.size(); i++) {
//            avg_vec = avg_vector(res[i]);
//            if (res[i].size() <= k) {
//            }
//            else {
//                if (avg_vec >= bound) {
//                    output_res.push_back(make_pair(res[i], avg_vec));
//                }
//            }
//        }
//    }
//
//    cout << "naive avg top-r end" << endl;
//
////    sort(output_res.begin(), output_res.end(), [](pair<vector<int>, double> & left, pair<vector<int>, double> & right){
////        return left.second > right.second;
////    });
//
//    // output result
//    for (unsigned int i = 0; i < r; i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }
//
//    endTime = clock();
//    cout << "k: " << k << " r: " << r << endl;
//    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//}

//void TopRComm::naive_avg_global_topr(const string & adjpath, const string & corepath,
//        const string & pagepath, const int k, const int r) {
//    clock_t startTime, endTime;
//    vector<vector<int>> res;
//    vector<pair<vector<int>, double>> temp;
//    vector<pair<vector<int>, double>> output_res;
//    vector<int> add_vec; // after delete one node, new vector
//    vector<int> k_core_vec;
//    vector<int> del_nodes_vec;
//    vector<pair<int, int>> batch_del_edges_vec;
//    vector<pair<int, int>> single_del_edges_vec;
//    vector<double> weight_vec;
//    double update_weight; // updated weight
//    double avg_vec;
//    degree_c = k;
//    size_c = r;
//    load_data(adjpath, corepath, pagepath);
//    cout << "successfully load data" << endl;
//    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
//    bool is_exist[graph_size];
//    for (int i = 0; i < graph_size; i++) {
//        is_exist[i] = false;
//    }
//    startTime = clock();
//    compute_max_k_core(degree_c);
//    cout << "compute_max_k_core" << endl;
//
//    res = compute_connected_components();
//    cout << "compute connected components" << endl;
//
//    if (res.size() == 0) {
//        return;
//    }
//
//    for (unsigned int i = 0; i < res.size(); i++) {
//        avg_vec = avg_vector(res[i]);
//        if (res[i].size() <= k) {
//        }
//        else {
//            output_res.push_back(make_pair(res[i], avg_vec));
//        }
//    }
//
//    cout << "init: " << output_res.size() << endl;
//
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//            is_exist[output_res[i].first[j]] = true;
//            k_core_vec.push_back(output_res[i].first[j]);
//        }
//    }
//
//    cout << "naive avg top-r begin" << endl;
//
//    for (int i = 0; i < graph_size; i++) {
//        if (is_exist[i] == true) {
//            for (int j = 0; j < output_res.size(); j++) {
//                if (output_res[j].first.size() <= k + 1) {
//                }
//                else {
//                    add_vec.assign(output_res[j].first.begin(), output_res[j].first.end());
//                    for (unsigned int m = 0; m < add_vec.size(); m++) {
//                        if (add_vec[m] == i) {
//                            add_vec.erase(add_vec.begin() + m);
//                            break;
//                        }
//                    }
//                    if (add_vec.size() == output_res[j].first.size()) {
//                    }
//                    else {
//                        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
//                        set<int> add_set(add_vec.begin(), add_vec.end());
//                        set<int> del_nodes_set;
//                        del_nodes_vec.clear();
//                        set_difference(k_core_set.begin(), k_core_set.end(),
//                                       add_set.begin(), add_set.end(),
//                                       inserter(del_nodes_set, del_nodes_set.begin()));
//                        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
//                        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
//                        if (is_connected_k_core_graph(k)) {
//                            update_weight = avg_vector(add_vec);
//                            output_res.push_back(make_pair(add_vec, update_weight));
//                        }
//                        else {
//                            vector<vector<int>> divide_com = compute_connected_k_core(k);
//                            for (unsigned int m = 0; m < divide_com.size(); m++) {
//                                double avg_divide_com = avg_vector(divide_com[m]);
//                                if (divide_com[m].size() <= k) {
//                                }
//                                else {
//                                    output_res.push_back(make_pair(divide_com[m], avg_divide_com));
//                                }
//                            }
//                        }
//                        batch_add_edges(batch_del_edges_vec);
//                    }
//                }
//            }
//        }
//    }
//
//    cout << "naive avg top-r end" << endl;
//
//    sort(output_res.begin(), output_res.end(), [](pair<vector<int>, double> & left, pair<vector<int>, double> & right){
//        return left.second > right.second;
//    });
//
//    // output result
//    for (unsigned int i = 0; i < r; i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }
//
//    endTime = clock();
//    cout << "k: " << k << " r: " << r << endl;
//    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//}

//void TopRComm::improved_avg_global_topr(const string & adjpath, const string & corepath,
//                                        const string & pagepath, const int k, const int r) {
//    clock_t startTime, endTime;
//    vector<vector<int>> res;
//    vector<pair<vector<int>, double>> temp;
//    vector<pair<vector<int>, double>> output_res;
//    vector<int> add_vec; // after delete one node, new vector
//    vector<int> k_core_vec;
//    vector<int> del_nodes_vec;
//    vector<pair<int, int>> batch_del_edges_vec;
//    vector<pair<int, int>> single_del_edges_vec;
//    vector<double> weight_vec;
//    double update_weight; // updated weight
//    double avg_vec;
//    double top_r_weight = std::numeric_limits<double>::max();
//    vector<double> sort_weight_vec;
//    degree_c = k;
//    size_c = r;
//    load_data(adjpath, corepath, pagepath);
//    cout << "successfully load data" << endl;
//    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
//    bool is_exist[graph_size];
//    for (int i = 0; i < graph_size; i++) {
//        is_exist[i] = false;
//    }
//    startTime = clock();
//    compute_max_k_core(degree_c);
//    cout << "compute_max_k_core" << endl;
//
//    res = compute_connected_components();
//    cout << "compute connected components" << endl;
//
//    if (res.size() == 0) {
//        return;
//    }
//
//    for (unsigned int i = 0; i < res.size(); i++) {
//        avg_vec = avg_vector(res[i]);
//        if (res[i].size() <= k) {
//        }
//        else {
//            if (output_res.size() < r) {
//                output_res.push_back(make_pair(res[i], avg_vec));
//                if (avg_vec < top_r_weight) {
//                    top_r_weight = avg_vec;
//                }
//            }
//            else {
//                if (avg_vec > top_r_weight) {
//                    output_res.push_back(make_pair(res[i], avg_vec));
//                    sort_weight_vec.clear();
//                    sort_weight_vec.resize(output_res.size());
//                    for (unsigned int j = 0; j < output_res.size(); j++) {
//                        sort_weight_vec[j] = output_res[j].second;
//                    }
//                    sort(sort_weight_vec.begin(), sort_weight_vec.end(), [](const double & left, const double & right){
//                        return left > right;
//                    });
//                    top_r_weight = sort_weight_vec[r];
//                    for (unsigned int j = 0; j < output_res.size(); ) {
//                        if (is_keep_avg_component(output_res[j].first, top_r_weight, k)) {
//                            j++;
//                        }
//                        else {
//                            output_res.erase(output_res.begin() + j);
//                        }
//                    }
//                }
//                else {
//                    if (is_keep_avg_component(res[i], top_r_weight, k)) {
//                        output_res.push_back(make_pair(res[i], avg_vec));
//                    }
//                }
//            }
//        }
//    }
//
//    cout << "init: " << output_res.size() << endl;
//
//    for (unsigned int i = 0; i < output_res.size(); i++) {
//        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
//            is_exist[output_res[i].first[j]] = true;
//            k_core_vec.push_back(output_res[i].first[j]);
//        }
//    }
//
//    cout << "improve avg top-r begin" << endl;
//
//    for (int i = 0; i < graph_size; i++) {
//        if (is_exist[i] == true) {
//            for (int j = 0; j < output_res.size(); j++) {
//                if (output_res[j].first.size() <= k + 1) {
//                }
//                else {
//                    add_vec.assign(output_res[j].first.begin(), output_res[j].first.end());
//                    for (unsigned int m = 0; m < add_vec.size(); m++) {
//                        if (add_vec[m] == i) {
//                            add_vec.erase(add_vec.begin() + m);
//                            break;
//                        }
//                    }
//                    if (add_vec.size() == output_res[j].first.size()) {
//                    }
//                    else {
//                        set<int> k_core_set(k_core_vec.begin(), k_core_vec.end());
//                        set<int> add_set(add_vec.begin(), add_vec.end());
//                        set<int> del_nodes_set;
//                        del_nodes_vec.clear();
//                        set_difference(k_core_set.begin(), k_core_set.end(),
//                                       add_set.begin(), add_set.end(),
//                                       inserter(del_nodes_set, del_nodes_set.begin()));
//                        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
//                        batch_del_edges_vec = batch_del_edges(del_nodes_vec);
//                        if (is_connected_k_core_graph(k)) {
//                            update_weight = avg_vector(add_vec);
//                            if (output_res.size() < r) {
//                                output_res.push_back(make_pair(add_vec, update_weight));
//                                if (update_weight < top_r_weight) {
//                                    top_r_weight = update_weight;
//                                }
//                            }
//                            else {
//                                if (update_weight > top_r_weight) {
//                                    output_res.push_back(make_pair(add_vec, update_weight));
//                                    sort_weight_vec.clear();
//                                    sort_weight_vec.resize(output_res.size());
//                                    for (unsigned int m = 0; m < output_res.size(); m++) {
//                                        sort_weight_vec[m] = output_res[m].second;
//                                    }
//                                    sort(sort_weight_vec.begin(), sort_weight_vec.end(),
//                                         [](const double & left, const double & right){
//                                             return left > right;
//                                         });
//                                    top_r_weight = sort_weight_vec[r];
//                                    for (unsigned int m = 0; m < output_res.size(); ) {
//                                        if (is_keep_avg_component(output_res[m].first, top_r_weight, k)) {
//                                            m++;
//                                        }
//                                        else {
//                                            output_res.erase(output_res.begin() + m);
//                                        }
//                                    }
//                                }
//                                else {
//                                    if (is_keep_avg_component(add_vec, top_r_weight, k)) {
//                                        output_res.push_back(make_pair(add_vec, update_weight));
//                                    }
//                                }
//                            }
//                        }
//                        else {
//                            vector<vector<int>> divide_com = compute_connected_k_core(k);
//                            for (unsigned int m = 0; m < divide_com.size(); m++) {
//                                double avg_divide_com = avg_vector(divide_com[m]);
//                                if (divide_com[m].size() <= k) {
//                                }
//                                else {
//                                    if (output_res.size() < r) {
//                                        output_res.push_back(make_pair(divide_com[m], avg_divide_com));
//                                        if (avg_divide_com < top_r_weight) {
//                                            top_r_weight = avg_divide_com;
//                                        }
//                                    }
//                                    else {
//                                        if (avg_divide_com > top_r_weight) {
//                                            output_res.push_back(make_pair(divide_com[m], avg_divide_com));
//                                            sort_weight_vec.clear();
//                                            sort_weight_vec.resize(output_res.size());
//                                            for (unsigned int n = 0; n < output_res.size(); n++) {
//                                                sort_weight_vec[n] = output_res[n].second;
//                                            }
//                                            sort(sort_weight_vec.begin(), sort_weight_vec.end(),
//                                                 [](const double & left, const double & right){
//                                                     return left > right;
//                                                 });
//                                            top_r_weight = sort_weight_vec[r];
//                                            for (unsigned int n = 0; n < output_res.size(); ) {
//                                                if (is_keep_avg_component(output_res[n].first, top_r_weight, k)) {
//                                                    n++;
//                                                }
//                                                else {
//                                                    output_res.erase(output_res.begin() + n);
//                                                }
//                                            }
//                                        }
//                                        else {
//                                            if (is_keep_avg_component(divide_com[m], top_r_weight, k)) {
//                                                output_res.push_back(make_pair(divide_com[m], avg_divide_com));
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        batch_add_edges(batch_del_edges_vec);
//                    }
//                }
//            }
//        }
//    }
//
//    cout << "improve avg top-r end" << endl;
//
//    sort(output_res.begin(), output_res.end(), [](pair<vector<int>, double> & left, pair<vector<int>, double> & right){
//        return left.second > right.second;
//    });
//
//    // output result
//    for (unsigned int i = 0; i < r; i++) {
//        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
//    }
//
//    endTime = clock();
//    cout << "k: " << k << " r: " << r << endl;
//    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//}

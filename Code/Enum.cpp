//
// Created by sbian on 2020/7/22.
//

#include "Enum.h"
#include "serialize.h"

EnumComm::EnumComm() {
    degree_c = -1; // the coreness of the community
    thres_c = -1.0; // the threshold constraint of the community
    graph_size = -1; // the original graph size
    edge_size = -1;
    c_size = -1; // the current graph size after some nodes removed
    max_core = 0;
    inf_value = 0.0;
}

EnumComm::~EnumComm() {
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

void EnumComm::init(const int numV) {
    graph_size = numV;
    c_size = numV;
    AdjList_Loc.resize(numV);
    for (int i = 0; i < graph_size; i++) {
        for (unsigned int j = 0; j < AdjList[i].size(); j++) {
            AdjList_Loc[i][AdjList[i][j]] = int(j);
        }
    }
}

void EnumComm::load_data(const string & adjpath, const string & corepath, const string & pagepath) {
    load_serialized_graph(adjpath, AdjList);
    load_serialized_graph(corepath, core);
    load_serialized_graph(pagepath, weight);
    graph_size = AdjList[AdjList.size() - 1][0];
    edge_size = AdjList[AdjList.size() - 1][1];
    cout << AdjList[AdjList.size() - 1][0] << " " << AdjList[AdjList.size() - 1][1] << " "
         << AdjList[AdjList.size() - 1][2] << " " << AdjList[AdjList.size() - 1][3] << " "
         << AdjList[AdjList.size() - 1][4] << endl;
    AdjList.pop_back();
    init(graph_size);
    assert(graph_size == AdjList.size());

    for(int i = 0; i < graph_size; i++) {
        if (core[i] > max_core) {
            max_core = core[i];
        }
    }
    cout << "max_core: " << max_core << endl;
}

void EnumComm::naive_max_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void EnumComm::naive_min_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void EnumComm::naive_sum_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {
    clock_t startTime, endTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> k_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    double sum_vec;
    double top_r_weight = std::numeric_limits<double>::max();
    int top_r_index = 0;
    bool flag = false; // determine whether duplicate
    degree_c = k;
    thres_c = tau;
    load_data(adjpath, corepath, pagepath);
    cout << graph_size << endl;
    bool is_exist[graph_size];
    for (int i = 0; i < graph_size; i++) {
        is_exist[i] = false;
    }
    startTime = clock();
    compute_max_k_core(degree_c);
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
        else if (sum_vec < tau) {
        }
        else {
            output_res.push_back(make_pair(res[i], sum_vec));
        }
    }

    cout << "init: " << output_res.size() << endl;

    for (unsigned int i = 0; i < output_res.size(); i++) {
        for (unsigned int j = 0; j < output_res[i].first.size(); j++) {
            is_exist[output_res[i].first[j]] = true;
            k_core_vec.push_back(output_res[i].first[j]);
        }
    }

    for (int i = 0; i < graph_size; i++) {
        if (is_exist[i] == true) {
            for (int j = 0; j < output_res.size(); j++) {
                if (output_res[j].first.size() <= k + 1) {
                }
                else if (output_res[j].second - weight[i] < tau) {
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
                        if (is_connected_k_core_graph(k)) {
                            update_weight = sum_vector(add_vec);
                            if (update_weight >= tau) {
                                output_res.push_back(make_pair(add_vec, update_weight));
                            }
                            cout << "update_weight: " << update_weight << endl;
                            cout << "final size: " << output_res.size() << endl;
                            cout << "top_r_weight: " << top_r_weight << endl;
                        }
                        else {
                            vector<vector<int>> divide_com = compute_connected_k_core(k);
                            for (unsigned int m = 0; m < divide_com.size(); m++) {
                                double sum_divide_com = sum_vector(divide_com[m]);
                                if (divide_com[m].size() <= k) {
                                }
                                else {
                                    if (sum_divide_com >= tau) {
                                        output_res.push_back(make_pair(divide_com[m], sum_divide_com));
                                    }
                                }
                            }
                        }
                        batch_add_edges(batch_del_edges_vec);
                    }
                }
            }
        }
    }

    // output result
    for (unsigned int i = 0; i < output_res.size(); i++) {
        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
    }

    endTime = clock();
    cout << "k: " << k << " tau: " << tau << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void EnumComm::improved_sum_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {
    clock_t startTime, endTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> temp;
    queue<pair<vector<int>, double>> L;
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
    thres_c = tau;
    load_data(adjpath, corepath, pagepath);
    cout << "successfully load data" << endl;
    cout << graph_size << endl;
    startTime = clock();
    compute_max_k_core(degree_c);
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
            if (sum_vec >= tau) {
                L.push(make_pair(res[i], sum_vec));
            }
        }
    }

    cout << "enumeration begin" << endl;

    while (L.size() != 0) {
        output_res.push_back(L.front());
        add_vec.assign(L.front().first.begin(), L.front().first.end());
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
        for (unsigned int i = 0; i < L.front().first.size(); i++) {
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
            else if (sum_vec < tau) {
                break;
            }
            else {
                single_del_edges_vec = single_del_edges(delete_node);
                if (is_connected_k_core_graph(k)) {
                    if (sum_vec >= tau) {
                        L.push(make_pair(add_vec, sum_vec));
                    }
                }
                else {
                    vector<vector<int>> divide_com = compute_connected_k_core(k);
                    for (unsigned int j = 0; j < divide_com.size(); j++) {
                        double sum_divide_com = sum_vector(divide_com[j]);
                        if (divide_com[j].size() <= k) {
                        }
                        else {
                            if (sum_divide_com >= tau) {
                                L.push(make_pair(divide_com[j], sum_divide_com));
                            }
                        }
                    }
                }
                batch_add_edges(single_del_edges_vec);
            }
            add_vec.insert(add_vec.begin() + node_weight[i].first, delete_node);
        }
        batch_add_edges(batch_del_edges_vec);
        L.pop();
        cout << "current size: " << output_res.size() << endl;
    }

    cout << "enumeration end" << endl;

    for (unsigned int i = 0; i < output_res.size(); i++) {
        cout << i + 1 << " " << setprecision(8) << output_res[i].second << " " << output_res[i].first.size() << endl;
    }

    endTime = clock();
    cout << "k: " << k << " tau: " << tau << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void EnumComm::naive_avg_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void EnumComm::improved_avg_global_enum(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void EnumComm::split(const string &src, const string &separator, vector<string> &dest) {
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

void EnumComm::add_edge(const int start, const int end) {
    if (start < end) {
        AdjList[start].push_back(end);
        AdjList_Loc[start][end] = AdjList[start].size() - 1;
        AdjList[end].push_back(start);
        AdjList_Loc[end][start] = AdjList[end].size() - 1;
    }
}

// one direction which means for edge (u,v), we only remove u, v but remain v, u
void EnumComm::remove_edge(const int start, const int end) {
    int idx = AdjList_Loc[start][end];
    AdjList_Loc[start][AdjList[start].back()] = idx;
    swap(AdjList[start][idx], AdjList[start].back());
    AdjList[start].pop_back();
}

void EnumComm::compute_max_k_core(const int k) {
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

vector<vector<int>> EnumComm::compute_connected_components() {
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

    delete [] visited;
    return res;
}

vector<vector<int>> EnumComm::compute_connected_k_core(const int k) {
    vector<vector<int>> res;
    int u;
    queue<int> q;
    bool visited [graph_size];
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
    return res;
}

bool EnumComm::is_connected_k_core_graph(const int k) {
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

double EnumComm::sum_vector(const vector<int> & vec) {
    double sum = 0.0;
    for (int i : vec) {
        sum += weight[i];
    }
    return sum;
}

vector<pair<int, int>> EnumComm::batch_del_edges(const vector<int> & nodes) {
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

vector<pair<int, int>> EnumComm::single_del_edges(const int node) {
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

void EnumComm::batch_add_edges(const vector<pair<int, int>> & edges) {
    for(unsigned int i = 0; i < edges.size(); i++) {
        if (edges[i].first < edges[i].second) {
            add_edge(edges[i].first, edges[i].second);
        }
    }
}


//
// Created by sbian on 2020/7/25.
//

#include "Local.h"
#include "serialize.h"

LocalComm::LocalComm() {
    degree_c = -1; // the coreness of the community
    thres_c = -1.0; // the threshold constraint of the community
    graph_size = -1; // the original graph size
    edge_size = -1;
    c_size = -1; // the current graph size after some nodes removed
    max_core = 0;
    inf_value = 0.0;
}

LocalComm::~LocalComm() {
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

void LocalComm::init(const int numV) {
    graph_size = numV;
    c_size = numV;
    AdjList_Loc.resize(numV);
    for (int i = 0; i < graph_size; i++) {
        for (unsigned int j = 0; j < AdjList[i].size(); j++) {
            AdjList_Loc[i][AdjList[i][j]] = int(j);
        }
    }
}

void LocalComm::load_data(const string & adjpath, const string & corepath, const string & pagepath) {
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

void LocalComm::naive_max_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void LocalComm::naive_min_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void LocalComm::naive_sum_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void LocalComm::improved_sum_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {

}

void LocalComm::self_avg_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {
    clock_t startTime, endTime;
    vector<vector<int>> res;
    vector<int> node_vec;
    int curr_max_size = 0;
    int curr_max_idx = 0;
    int next_max_size = 0;
    int next_max_idx = 0;
    degree_c = k;
    thres_c = tau;
    load_data(adjpath, corepath, pagepath);
    cout << graph_size << endl;
    startTime = clock();
    res = preprocess_local_avg(k, tau);

    if (res.size() == 0) {
        return;
    }

    curr_max_size = 0;
    curr_max_idx = 0;
    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() > curr_max_size) {
            curr_max_size = res[i].size();
            curr_max_idx = i;
        }
    }

    vector<vector<int>> divide_components;
    while (true) {
        node_vec = res[curr_max_idx];
        assert(node_vec.size() == curr_max_size);
        res.erase(res.begin() + curr_max_idx);
        divide_components = self_greedy(node_vec, k, tau);
        curr_max_idx = 0;
        curr_max_size = 0;
        for (unsigned int i = 0; i < divide_components.size(); i++) {
            if (divide_components[i].size() > curr_max_size) {
                curr_max_size = divide_components[i].size();
                curr_max_idx = i;
            }
        }
        next_max_idx = 0;
        next_max_size = 0;
        for (unsigned int i = 0; i < res.size(); i++) {
            if (res[i].size() > next_max_size) {
                next_max_size = res[i].size();
                next_max_idx = i;
            }
        }
        if (curr_max_size > next_max_size) {
            break;
        }
        else {
            res.insert(res.end(), divide_components.begin(), divide_components.end());
        }
        curr_max_size = next_max_size;
        curr_max_idx = next_max_idx;
    }

    endTime = clock();
    cout << "k: " << k << " tau: " << tau << endl;
    cout << "max size: " << curr_max_size << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void LocalComm::neighbor_avg_local(const string & adjpath, const string & corepath,
        const string & pagepath, const int k, const double tau) {
    clock_t startTime, endTime;
    vector<vector<int>> res;
    vector<int> node_vec;
    int curr_max_size = 0;
    int curr_max_idx = 0;
    int next_max_size = 0;
    int next_max_idx = 0;
    degree_c = k;
    thres_c = tau;
    load_data(adjpath, corepath, pagepath);
    cout << graph_size << endl;

    startTime = clock();
    res = preprocess_local_avg(k, tau);

    if (res.size() == 0) {
        return;
    }

    curr_max_size = 0;
    curr_max_idx = 0;
    for (unsigned int i = 0; i < res.size(); i++) {
        if (res[i].size() > curr_max_size) {
            curr_max_size = res[i].size();
            curr_max_idx = i;
        }
    }

    vector<vector<int>> divide_components;
    while (true) {
        node_vec = res[curr_max_idx];
        assert(node_vec.size() == curr_max_size);
        res.erase(res.begin() + curr_max_idx);
        divide_components = neighbor_greedy(node_vec, k, tau);
        curr_max_idx = 0;
        curr_max_size = 0;
        for (unsigned int i = 0; i < divide_components.size(); i++) {
            if (divide_components[i].size() > curr_max_size) {
                curr_max_size = divide_components[i].size();
                curr_max_idx = i;
            }
        }
        next_max_idx = 0;
        next_max_size = 0;
        for (unsigned int i = 0; i < res.size(); i++) {
            if (res[i].size() > next_max_size) {
                next_max_size = res[i].size();
                next_max_idx = i;
            }
        }
        if (curr_max_size > next_max_size) {
            break;
        }
        else {
            res.insert(res.end(), divide_components.begin(), divide_components.end());
        }
        curr_max_size = next_max_size;
        curr_max_idx = next_max_idx;
    }

    endTime = clock();
    cout << "k: " << k << " tau: " << tau << endl;
    cout << "max size: " << curr_max_size << endl;
    cout << "Total Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void LocalComm::split(const string & src, const string & separator, vector<string> & dest) {
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

void LocalComm::add_edge(const int start, const int end) {
    if (start < end) {
        AdjList[start].push_back(end);
        AdjList_Loc[start][end] = AdjList[start].size() - 1;
        AdjList[end].push_back(start);
        AdjList_Loc[end][start] = AdjList[end].size() - 1;
    }
}

void LocalComm::remove_edge(const int start, const int end) {
    int idx = AdjList_Loc[start][end];
    AdjList_Loc[start][AdjList[start].back()] = idx;
    swap(AdjList[start][idx], AdjList[start].back());
    AdjList[start].pop_back();
}

vector<vector<int>> LocalComm::preprocess_local_avg(const int k, const double tau) {
    vector<vector<int>> res;
    int u, v;
    double neigh_k_sum;
    double neigh_k_avg;
    queue<int> q;
    vector<double> neighbor_weight_vec;
    bool visited [graph_size];
    vector<int> id_vec;

    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        neighbor_weight_vec.clear();
        for (unsigned j = 0; j < AdjList[i].size(); j++) {
            v = AdjList[i][j];
            neighbor_weight_vec.push_back(weight[v]);
        }
        sort(neighbor_weight_vec.begin(), neighbor_weight_vec.end(), greater<double>());
        neigh_k_sum = 0.0;
        for (unsigned int j = 0; j < k; j++) {
            neigh_k_sum += neighbor_weight_vec[j];
        }
        neigh_k_avg = neigh_k_sum * 1.0 / k;
        if (visited[i] == false && (AdjList[i].size() < k || neigh_k_avg < tau)) {
            q.push(i);
            visited[i] = true;
            while (q.size() != 0) {
                u = q.front();
                q.pop();
                for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                    v = AdjList[u][j];
                    remove_edge(v, u);
                    neighbor_weight_vec.clear();
                    for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                        neighbor_weight_vec.push_back(weight[AdjList[v][m]]);
                    }
                    sort(neighbor_weight_vec.begin(), neighbor_weight_vec.end(), greater<double>());
                    neigh_k_sum = 0.0;
                    for (unsigned int m = 0; m < k; m++) {
                        neigh_k_sum += neighbor_weight_vec[m];
                    }
                    neigh_k_avg = neigh_k_sum * 1.0 / k;
                    if (visited[v] == false && (AdjList[v].size() < k || neigh_k_avg < tau)) {
                        q.push(v);
                        visited[v] = true;
                    }
                }
                AdjList[u].clear();
            }
        }
    }

    assert(q.size() == 0);

    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        id_vec.clear();
        if (visited[i] == false && AdjList[i].size() != 0) {
            q.push(i);
            visited[i] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                    v = AdjList[u][j];
                    if (visited[v] == false && AdjList[v].size() != 0) {
                        q.push(v);
                        visited[v] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

    return res;
}

vector<vector<int>> LocalComm::self_greedy(vector<int> & vec, const int k, const double tau) {
    vector<vector<int>> res;
    int u, v;
    queue<int> q;
    double neigh_sum;
    double neigh_avg;
    vector<int> id_vec;
    bool visited[graph_size];

    for (int i = 0; i < graph_size; i++) {
        visited[i] = true;
    }

    for (unsigned int i = 0; i < vec.size(); i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        if (visited[i] == false) {
            neigh_sum = 0.0;
            for (unsigned int j = 0; j < AdjList[i].size(); j++) {
                v = AdjList[i][j];
                neigh_sum += weight[v];
            }
            neigh_avg = neigh_sum * 1.0 / AdjList[i].size();
            if (AdjList[i].size() < k || neigh_avg < tau) {
                q.push(i);
                visited[i] = true;
                while (q.size() != 0) {
                    u = q.front();
                    q.pop();
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        remove_edge(v, u);
                        neigh_sum = 0.0;
                        for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                            neigh_sum += weight[AdjList[v][m]];
                        }
                        neigh_avg = neigh_sum * 1.0 / AdjList[v].size();
                        if (visited[v] == false && (AdjList[v].size() < k || neigh_avg < tau)) {
                            q.push(v);
                            visited[v] = true;
                        }
                    }
                    AdjList[u].clear();
                }
            }
        }
    }

    assert(q.size() == 0);

    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        id_vec.clear();
        if (visited[i] == false && AdjList[i].size() != 0) {
            q.push(i);
            visited[i] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                    v = AdjList[u][j];
                    if (visited[v] == false && AdjList[v].size() != 0) {
                        q.push(v);
                        visited[v] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

    return res;
}

vector<vector<int>> LocalComm::neighbor_greedy(vector<int> & vec, const int k, const double tau) {
    vector<vector<int>> res;
    int u, v;
    queue<int> q;
    int neigh_min_idx;
    double neigh_min_weight;
    double neigh_sum;
    double neigh_avg;
    vector<int> id_vec;
    bool visited[graph_size];

    for (int i = 0; i < graph_size; i++) {
        visited[i] = true;
    }

    for (unsigned int i = 0; i < vec.size(); i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        if (visited[i] == false) {
            neigh_sum = 0.0;
            for (unsigned int j = 0; j < AdjList[i].size(); j++) {
                v = AdjList[i][j];
                neigh_sum += weight[v];
            }
            neigh_avg = neigh_sum * 1.0 / AdjList[i].size();
            if (AdjList[i].size() < k) {
                q.push(i);
                visited[i] = true;
                while (q.size() != 0) {
                    u = q.front();
                    q.pop();
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        remove_edge(v, u);
                        neigh_sum = 0.0;
                        for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                            neigh_sum += weight[AdjList[v][m]];
                        }
                        neigh_avg = neigh_sum * 1.0 / AdjList[v].size();
                        if (visited[v] == false) {
                            if (AdjList[v].size() < k) {
                                q.push(v);
                                visited[v] = true;
                            }
                            else if (neigh_avg < tau) {
                                neigh_min_idx = 0;
                                neigh_min_weight = std::numeric_limits<double>::max();
                                for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                                    if (weight[AdjList[v][m]] < neigh_min_weight) {
                                        neigh_min_weight = weight[AdjList[v][m]];
                                        neigh_min_idx = AdjList[v][m];
                                    }
                                }
                                q.push(neigh_min_idx);
                                visited[neigh_min_idx] = true;
                            }
                        }
                    }
                    AdjList[u].clear();
                }
            }
            else if (neigh_avg < tau) {
                neigh_min_idx = 0;
                neigh_min_weight = std::numeric_limits<double>::max();
                for (unsigned int j = 0; j < AdjList[i].size(); j++) {
                    v = AdjList[i][j];
                    if (weight[v] < neigh_min_weight) {
                        neigh_min_weight = weight[v];
                        neigh_min_idx = v;
                    }
                }
                q.push(neigh_min_idx);
                visited[neigh_min_idx] = true;
                while (q.size() != 0) {
                    u = q.front();
                    q.pop();
                    for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                        v = AdjList[u][j];
                        remove_edge(v, u);
                        neigh_sum = 0.0;
                        for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                            neigh_sum += weight[AdjList[v][m]];
                        }
                        neigh_avg = neigh_sum * 1.0 / AdjList[v].size();
                        if (visited[v] == false) {
                            if (AdjList[v].size() < k) {
                                q.push(v);
                                visited[v] = true;
                            }
                            else if (neigh_avg < tau) {
                                neigh_min_idx = 0;
                                neigh_min_weight = std::numeric_limits<double>::max();
                                for (unsigned int m = 0; m < AdjList[v].size(); m++) {
                                    if (weight[AdjList[v][m]] < neigh_min_weight) {
                                        neigh_min_weight = weight[AdjList[v][m]];
                                        neigh_min_idx = AdjList[v][m];
                                    }
                                }
                                q.push(neigh_min_idx);
                                visited[neigh_min_idx] = true;
                            }
                        }
                    }
                    AdjList[u].clear();
                }
            }
        }
    }

    assert(q.size() == 0);

    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {
        id_vec.clear();
        if (visited[i] == false && AdjList[i].size() != 0) {
            q.push(i);
            visited[i] = true;
            while (q.size() != 0) {
                u = q.front();
                id_vec.push_back(u);
                q.pop();
                for (unsigned int j = 0; j < AdjList[u].size(); j++) {
                    v = AdjList[u][j];
                    if (visited[v] == false && AdjList[v].size() != 0) {
                        q.push(v);
                        visited[v] = true;
                    }
                }
            }
            res.push_back(id_vec);
        }
    }

    return res;
}

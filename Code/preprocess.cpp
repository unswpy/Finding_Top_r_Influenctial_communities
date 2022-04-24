//
// Created by sbian on 2020/7/25.
//

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <list>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;

void preprocess(const string & file);
void split(const string& src, const string& separator, vector<string>& dest);
vector<int> compute_core(int graph_size, vector<list<int>> AdjList);
vector<double> compute_pagerank(int graph_size, double alpha, const vector<list<int>> & AdjList);

int main() {
//    string path = "D:\\Git\\Weight_Community\\origin\\com-dblp.ungraph.txt";
    vector<string> file_vec = {"dblp", "friendster", "lj", "orkut", "youtube", "email"};
    for (unsigned int i = 0; i < file_vec.size(); i++) {
        preprocess(file_vec[i]);
    }
}

void preprocess(const string & file) {
    string inpath, outpath;
    if (file == "dblp") {
        inpath = "origin/com-dblp.ungraph.txt";
        outpath = "cpreprocess/com-dblp.txt";
    }
    else if (file == "lj") {
        inpath = "origin/com-lj.ungraph.txt";
        outpath = "cpreprocess/com-lj.txt";
    }
    else if (file == "orkut") {
        inpath = "origin/com-orkut.ungraph.txt";
        outpath = "cpreprocess/com-orkut.txt";
    }
    else if (file == "email") {
        inpath = "origin/Email-Enron.txt";
        outpath = "cpreprocess/email.txt";
    }
    else if (file == "youtube") {
        inpath = "origin/com-youtube.ungraph.txt";
        outpath = "cpreprocess/com-youtube.txt";
    }
    else if (file == "friendster") {
        inpath = "origin/com-friendster.ungraph.txt";
        outpath = "cpreprocess/com-friendster.txt";
    }
    int count = 0;
    int graph_size = 0;
    int edge_size = 0;
    int start = 0;
    int end = 0;
    int count_node_id = 0;
    double d_avg;
    int d_sum = 0;
    int d_max = 0;
    int k_max = 0;
    map<int, int> node_id_map;
    vector<list<int>> AdjList;
    vector<string> svec;
    vector<int> core_vec;
    vector<double> pagerank_vec;
    char buf[1024];
    string message;
    ifstream infile;
    infile.open(inpath);
    if (infile.is_open()) {
        while (infile.good() && !infile.eof()) {
            memset(buf, 0, 1024);
            infile.getline(buf, 1024);
            message = buf;
            if (count == 2) {
                split(message, " ", svec);
                graph_size = stoi(svec[2]);
                edge_size = stoi(svec[4]);
                svec.clear();
                AdjList.resize(graph_size);
                cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
            }
            else if (count > 3) {
                split(message, "\t", svec);
                if (svec.size() == 2) {
                    start = stoi(svec[0]);
                    end = stoi(svec[1]);
                    if (node_id_map.count(start) == 0) {
                        node_id_map[start] = count_node_id;
                        start = node_id_map[start];
                        count_node_id += 1;
                    }
                    else {
                        start = node_id_map[start];
                    }
                    if (node_id_map.count(end) == 0) {
                        node_id_map[end] = count_node_id;
                        end = node_id_map[end];
                        count_node_id += 1;
                    }
                    else {
                        end = node_id_map[end];
                    }
                    AdjList[start].push_back(end);
                    AdjList[end].push_back(start);
                }
                else {
                    cout << svec.size() << endl;
                }
                svec.clear();
            }
            count += 1;
        }
        infile.close();
    }

    for (int i = 0; i < graph_size; i++) {
        d_sum += AdjList[i].size();
        if (d_max < AdjList[i].size()) {
            d_max = AdjList[i].size();
        }
    }
    d_avg = d_sum * 1.0 / graph_size;

    core_vec = compute_core(graph_size, AdjList);
    for (unsigned int i = 0; i < core_vec.size(); i++) {
        if (k_max < core_vec[i]) {
            k_max = core_vec[i];
        }
    }
    pagerank_vec = compute_pagerank(graph_size, 0.85, AdjList);

    // write to preprocess_file
    ofstream outfile;
    outfile.open(outpath);
    outfile << graph_size << "  " << edge_size << "  " << d_max << "  " << d_avg << "  " << k_max << "  " << endl;
    for (int i = 0; i < graph_size; i++) {
        outfile << i << "  " << setprecision(10) << std::fixed << pagerank_vec[i] << "  " << core_vec[i] << endl;
    }
    for (int i = 0; i < graph_size; i++) {
        for (auto iter = AdjList[i].begin(); iter != AdjList[i].end(); iter++) {
            outfile << i << "  " << *iter << endl;
        }
    }
    outfile.close();
}

void split(const string& src, const string& separator, vector<string>& dest)
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

vector<int> compute_core(int graph_size, vector<list<int>> AdjList) {
    vector<int> res;
    res.resize(graph_size);
    vector<pair<int, int>> nodes;
    pair<int, int> temp;
    for (int i = 0; i < graph_size; i++) {
        nodes.push_back(make_pair(i, AdjList[i].size()));
    }
    sort(nodes.begin(), nodes.end(), [](const pair<int, int> & left, const pair<int, int> & right){
        return left.second < right.second;
    });
    vector<int> bin_boundaries;
    bin_boundaries.push_back(0);
    int curr_degree = 0;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        if (nodes[i].second > curr_degree) {
            for (int j = 0; j < nodes[i].second - curr_degree; j++) {
                bin_boundaries.push_back(i);
            }
            curr_degree = nodes[i].second;
        }
    }
    map<int, int> node_pos;
    map<int, int> core;
    int pos;
    int bin_start;
    int v;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        node_pos[nodes[i].first] = i;
        core[nodes[i].first] = nodes[i].second;
    }
    for (unsigned int i = 0; i < nodes.size(); i++) {
        v = nodes[i].first;
        for (auto u = AdjList[v].begin(); u != AdjList[v].end(); u++) {
            if (core[*u] > core[v]) {
                AdjList[*u].remove(v);
                pos = node_pos[*u];
                bin_start = bin_boundaries[core[*u]];
                node_pos[*u] = bin_start;
                node_pos[nodes[bin_start].first] = pos;
                temp = nodes[pos];
                nodes[pos] = nodes[bin_start];
                nodes[bin_start] = temp;
                bin_boundaries[core[*u]] += 1;
                core[*u] -= 1;
            }
        }
    }
    for (int i = 0; i < graph_size; i++) {
        res[i] = core[i];
    }
    return res;
}

vector<double> compute_pagerank(int graph_size, double alpha, const vector<list<int>> & AdjList) {
    vector<double> res;
    res.resize(graph_size);
    vector<int> dangling_nodes;
    vector<double> x;
    x.resize(graph_size);
    vector<double> p;
    p.resize(graph_size);
    vector<double> xlast;
    vector<double> dangling_weights;
    double danglesum = 0.0;
    double xlast_sum;
    double err;
    for (int i = 0; i < graph_size; i++) {
        x[i] = 1.0 / graph_size;
        p[i] = 1.0 / graph_size;
    }
    dangling_weights.assign(p.begin(), p.end());

    for (int i = 0; i < graph_size; i++) {
        if (AdjList[i].size() == 0) {
            dangling_nodes.push_back(i);
        }
    }

    for (int i = 0; i < 100; i++) {
        xlast.assign(x.begin(), x.end());
        for (unsigned int j = 0; j < x.size(); j++) {
            x[j] = 0.0;
        }
        xlast_sum = 0.0;
        for (unsigned int j = 0; j < dangling_nodes.size(); j++) {
             xlast_sum += xlast[dangling_nodes[j]];
        }
        danglesum = alpha * xlast_sum;
        for (unsigned int j = 0; j < x.size(); j++) {
            for (auto iter = AdjList[j].begin(); iter != AdjList[j].end(); iter++) {
                x[*iter] += alpha * xlast[j] * 1.0 / AdjList[j].size();
            }
            x[j] = x[j] + danglesum * dangling_weights[j] + (1.0 - alpha) * p[j];
        }
        err = 0.0;
        for (unsigned int j = 0; j < x.size(); j++) {
            err += fabs(x[j] - xlast[j]);
        }
        if (err < graph_size * 0.0000001) {
            for (unsigned int j = 0; j < x.size(); j++) {
                res[j] = x[j];
            }
        }
    }
    return res;
}


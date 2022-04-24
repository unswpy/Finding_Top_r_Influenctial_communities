//
// Created by sbian on 2020/7/27.
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
#include "file_ctrl.h"

using namespace std;

void preprocess(const string & file);
void split(const string& src, const string& separator, vector<string>& dest);
vector<int> compute_core(int graph_size, vector<vector<int>> AdjList);
vector<double> compute_pagerank(int graph_size, double alpha, const vector<vector<int>> & AdjList);

int main() {
//    string path = "D:\\Git\\Weight_Community\\origin\\com-dblp.ungraph.txt";
//    vector<string> file_vec = {"dblp", "lj", "orkut", "youtube", "email"};
//    vector<string> file_vec = {"friendster"};
//    vector<string> file_vec = {"domain_case", "domain_label", "domain_pub"};
    vector<string> file_vec = {"topic_case", "topic_label", "topic_pub"};
    for (unsigned int i = 0; i < file_vec.size(); i++) {
        preprocess(file_vec[i]);
    }
}

void preprocess(const string & file) {
    string inpath, adjpath, corepath, pagepath;
    string pubpath, labelpath;
    if (file == "dblp") {
        inpath = "origin/com-dblp.ungraph.txt";
        adjpath = "spreprocess/com-dblp-adj.file";
        corepath = "spreprocess/com-dblp-core.file";
        pagepath = "spreprocess/com-dblp-page.file";
    }
    else if (file == "dblp_case") {
        inpath = "origin/dblp_case.txt";
        adjpath = "spreprocess/dblp-case-adj.file";
        corepath = "spreprocess/dblp-case-core.file";
        pagepath = "spreprocess/dblp-case-page.file";
    }
    else if (file == "domain_case") {
        inpath = "origin/domain_case.txt";
        adjpath = "spreprocess/domain-case-adj.file";
        corepath = "spreprocess/domain-case-core.file";
        pagepath = "spreprocess/domain-case-page.file";
    }
    else if (file == "domain_pub") {
        inpath = "origin/domain_case.txt";
        pubpath = "origin/domain_pub_list.txt";
        adjpath = "spreprocess/domain-pub-adj.file";
        corepath = "spreprocess/domain-pub-core.file";
        pagepath = "spreprocess/domain-pub-page.file";
    }
    else if (file == "domain_label") {
        inpath = "origin/domain_case.txt";
        labelpath = "origin/domain_label_list.txt";
        adjpath = "spreprocess/domain-label-adj.file";
        corepath = "spreprocess/domain-label-core.file";
        pagepath = "spreprocess/domain-label-page.file";
    }
    else if (file == "topic_case") {
        inpath = "origin/topic_case.txt";
        adjpath = "spreprocess/topic-case-adj.file";
        corepath = "spreprocess/topic-case-core.file";
        pagepath = "spreprocess/topic-case-page.file";
    }
    else if (file == "topic_pub") {
        inpath = "origin/topic_case.txt";
        pubpath = "origin/topic_pub_list.txt";
        adjpath = "spreprocess/topic-pub-adj.file";
        corepath = "spreprocess/topic-pub-core.file";
        pagepath = "spreprocess/topic-pub-page.file";
    }
    else if (file == "topic_label") {
        inpath = "origin/topic_case.txt";
        labelpath = "origin/topic_label_list.txt";
        adjpath = "spreprocess/topic-label-adj.file";
        corepath = "spreprocess/topic-label-core.file";
        pagepath = "spreprocess/topic-label-page.file";
    }
    else if (file == "lj") {
        inpath = "origin/com-lj.ungraph.txt";
        adjpath = "spreprocess/com-lj-adj.file";
        corepath = "spreprocess/com-lj-core.file";
        pagepath = "spreprocess/com-lj-page.file";
    }
    else if (file == "orkut") {
        inpath = "origin/com-orkut.ungraph.txt";
        adjpath = "spreprocess/com-orkut-adj.file";
        corepath = "spreprocess/com-orkut-core.file";
        pagepath = "spreprocess/com-orkut-page.file";
    }
    else if (file == "email") {
        inpath = "origin/Email-Enron.txt";
        adjpath = "spreprocess/email-adj.file";
        corepath = "spreprocess/email-core.file";
        pagepath = "spreprocess/email-page.file";
    }
    else if (file == "youtube") {
        inpath = "origin/com-youtube.ungraph.txt";
        adjpath = "spreprocess/com-youtube-adj.file";
        corepath = "spreprocess/com-youtube-core.file";
        pagepath = "spreprocess/com-youtube-page.file";
    }
    else if (file == "friendster") {
        inpath = "origin/com-friendster.ungraph.txt";
        adjpath = "spreprocess/com-friendster-adj.file";
        corepath = "spreprocess/com-friendster-core.file";
        pagepath = "spreprocess/com-friendster-page.file";
    }
    int count = 0;
    int graph_size = 0;
    int edge_size = 0;
    int edge_count = 0;
    int start = 0;
    int end = 0;
    int count_node_id = 0;
    double d_avg;
    int d_sum = 0;
    int d_max = 0;
    int k_max = 0;
    map<int, int> node_id_map;
    vector<vector<int>> AdjList;
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
                AdjList.resize(graph_size + 1);
                cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
            }
            else if (count > 3) {
                split(message, "\t", svec);
                if (svec.size() == 2) {
                    start = stoi(svec[0]);
                    end = stoi(svec[1]);
                    if (file != "dblp_case" &&
                        file != "domain_case" &&
                        file != "domain_pub" &&
                        file != "domain_label") {
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
                    }
                    else {
                        if (node_id_map.count(start) == 0) {
                            node_id_map[start] = count_node_id;
                            count_node_id += 1;
                        }
                        if (node_id_map.count(end) == 0) {
                            node_id_map[end] = count_node_id;
                            count_node_id += 1;
                        }
                    }
                    if (file == "email") {
                        if (start < end) {
                            AdjList[start].push_back(end);
                            AdjList[end].push_back(start);
                            edge_count += 1;
                        }
                    }
                    else {
                        AdjList[start].push_back(end);
                        AdjList[end].push_back(start);
                        edge_count += 1;
                    }
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
    cout << edge_size << " " << edge_count << endl;
    cout << "line number: " << count << endl;
    cout << "node id: " << count_node_id << endl;
    if (file == "email") {
        assert(graph_size == count_node_id);
        assert(edge_size == edge_count * 2);
    }
    else {
        assert(graph_size == count_node_id);
        assert(edge_size == edge_count);
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
    AdjList[graph_size].push_back(graph_size);
    AdjList[graph_size].push_back(edge_size);
    AdjList[graph_size].push_back(d_max);
    AdjList[graph_size].push_back(d_avg);
    AdjList[graph_size].push_back(k_max);
    save_serialized_graph(adjpath, AdjList);
    save_serialized_graph(corepath, core_vec);

    if (file != "domain_pub" && file != "domain_label") {
        pagerank_vec = compute_pagerank(graph_size, 0.85, AdjList);
    }
    else if (file == "domain_pub"){
        infile.open(pubpath);
        if (infile.is_open()) {
            while (infile.good() && !infile.eof()) {
                memset(buf, 0, 1024);
                infile.getline(buf, 1024);
                message = buf;
                split(message, " ", svec);
                if (svec.size() == 2) {
                    pagerank_vec.push_back(stod(svec[1]));
                }
                svec.clear();
            }
            infile.close();
        }
    }
    else if (file == "domain_label") {
        infile.open(labelpath);
        if (infile.is_open()) {
            while (infile.good() && !infile.eof()) {
                memset(buf, 0, 1024);
                infile.getline(buf, 1024);
                message = buf;
                split(message, " ", svec);
                if (svec.size() == 2) {
                    pagerank_vec.push_back(stod(svec[1]));
                }
                svec.clear();
            }
            infile.close();
        }
    }
    save_serialized_graph(pagepath, pagerank_vec);

    // write to serial preprocess_file


    // write to preprocess_file
//    ofstream outfile;
//    outfile.open(outpath);
//    outfile << graph_size << "  " << edge_size << "  " << d_max << "  " << d_avg << "  " << k_max << "  " << endl;
//    for (int i = 0; i < graph_size; i++) {
//        outfile << i << "  " << setprecision(10) << std::fixed << pagerank_vec[i] << "  " << core_vec[i] << endl;
//    }
//    for (int i = 0; i < graph_size; i++) {
//        for (auto iter = AdjList[i].begin(); iter != AdjList[i].end(); iter++) {
//            outfile << i << "  " << *iter << endl;
//        }
//    }
//    outfile.close();
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

vector<int> compute_core(int graph_size, vector<vector<int>> AdjList) {
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
    int u, v;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        node_pos[nodes[i].first] = i;
        core[nodes[i].first] = nodes[i].second;
    }
    for (unsigned int i = 0; i < nodes.size(); i++) {
        v = nodes[i].first;
        for (unsigned int j = 0; j < AdjList[v].size(); j++) {
            u = AdjList[v][j];
            if (core[u] > core[v]) {
                for (unsigned int k = 0; k < AdjList[u].size(); k++) {
                    if (AdjList[u][k] == v) {
                        AdjList[u].erase(AdjList[u].begin() + k);
                        break;
                    }
                }
                pos = node_pos[u];
                bin_start = bin_boundaries[core[u]];
                node_pos[u] = bin_start;
                node_pos[nodes[bin_start].first] = pos;
                temp = nodes[pos];
                nodes[pos] = nodes[bin_start];
                nodes[bin_start] = temp;
                bin_boundaries[core[u]] += 1;
                core[u] -= 1;
            }
        }
    }
    for (int i = 0; i < graph_size; i++) {
        res[i] = core[i];
    }
    return res;
}

vector<double> compute_pagerank(int graph_size, double alpha, const vector<vector<int>> & AdjList) {
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


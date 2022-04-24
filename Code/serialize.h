//
// Created by sbian on 2020/7/27.
//

#ifndef WEIGHT_COMMUNITY_SERIALIZE_H
#define WEIGHT_COMMUNITY_SERIALIZE_H

#include <vector>
#include <tuple>
#include <numeric>

#include <cassert>

#include <fstream>
#include <cstdlib>

#include <iostream>

using namespace std;

typedef std::vector<uint8_t> StreamType;

template <class T>
void serialize(const T&, StreamType&);


template <size_t>
struct int_
{
};


// get_size
template <class T>
size_t get_size(const T& obj);


template <class T>
struct get_size_helper;

template <class tuple_type>
size_t get_tuple_size(const tuple_type& obj, int_<0>);

template <class tuple_type, size_t pos>
size_t get_tuple_size(const tuple_type& obj, int_<pos>);




template <class T>
struct serialize_helper;

template <class T>
void serializer(const T& obj, StreamType::iterator&);

template <class tuple_type>
void serialize_tuple(const tuple_type& obj, StreamType::iterator& res, int_<0>);

template <class tuple_type, size_t pos>
void serialize_tuple(const tuple_type& obj, StreamType::iterator& res, int_<pos>);



template <class T>
struct deserialize_helper;

template <class tuple_type>
void deserialize_tuple(tuple_type& obj,
                       StreamType::const_iterator& begin,
                       StreamType::const_iterator end, int_<0>);

template <class tuple_type, size_t pos>
void deserialize_tuple(tuple_type& obj,
                       StreamType::const_iterator& begin,
                       StreamType::const_iterator end, int_<pos>);


template <class T>
T deserialize(StreamType::const_iterator& begin, const StreamType::const_iterator& end);

template <class T>
T deserialize(const StreamType& res);

void split_filename(const std::string& filepath, std::string &dir, std::string &filename);
int make_dir(const std::string& filepath);

template <class T>
static void save_file(const std::string filename, const T& output);

template <class T>
static void load_file(const std::string filename, T& input);

template <class T>
void save_serialized_graph(const std::string file_name, const T& graph);

template <class T>
void load_serialized_graph(const std::string file_name, T& graph);

bool seraizlied_graph_exist(const std::string file_name);

#endif //WEIGHT_COMMUNITY_SERIALIZE_H

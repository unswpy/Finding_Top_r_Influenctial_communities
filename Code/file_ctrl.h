//
// Created by sbian on 2020/7/27.
//

#ifndef WEIGHT_COMMUNITY_FILE_CTRL_H
#define WEIGHT_COMMUNITY_FILE_CTRL_H

#include <fstream>
#include <cstdlib>
#include "serialize_qt.h"

void split_filename(const std::string& filepath, std::string &dir, std::string &filename)
{
    std::size_t found = filepath.find_last_of("/");
    if (found ==  std::string::npos)
    {
        dir = "";
        filename = filepath;
    }
    else
    {
        dir = filepath.substr(0, found);
        filename = filepath.substr(found + 1);
    }

    std::cout << " path: " << dir << '\n';
    std::cout << " file: " << filename << '\n';
}

int make_dir(const std::string& filepath)
{
    std::string dir;
    std::string filename;
    split_filename(filepath, dir, filename);
    if (dir == "")
    {
        return 0;
    }

    std::string cmd = "mkdir -p " + dir;
    return std::system(cmd.c_str());
}

/// Save a serialized file
template <class T>
static void save_file(const std::string filename, const T& output)
{
    std::ofstream outfile(filename, std::ios::binary);

    if (!outfile.eof() && !outfile.fail())
    {
        StreamType res;
        serialize(output, res);
        outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
        outfile.close();
        res.clear();
        std::cout << "Save file successfully: " << filename << '\n';
    }
    else
    {
        std::cout << "Save file failed: " + filename << '\n';
        exit(1);
    }
}

/// Load a serialized file
template <class T>
static void load_file(const std::string filename, T& input)
{
    std::ifstream infile(filename, std::ios::binary);

    if (!infile.eof() && !infile.fail())
    {
        infile.seekg(0, std::ios_base::end);
        const std::streampos fileSize = infile.tellg();
        infile.seekg(0, std::ios_base::beg);
        std::vector<uint8_t> res(fileSize);
        infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
        infile.close();
        input.clear();
        auto it = res.cbegin();
        input = deserialize<T>(it, res.cend());
        res.clear();
    }
    else
    {
        std::cout << "Cannot open file: " + filename << '\n';
        exit(1);
    }
}

/// Save graph structure to a file
template <class T>
void save_serialized_graph(const std::string file_name, const T& graph)
{
    make_dir(file_name);
    save_file(file_name, graph);
}

/// Load graph structure from a file
template <class T>
void load_serialized_graph(const std::string file_name, T& graph)
{
    load_file(file_name, graph);
}

bool seraizlied_graph_exist(const std::string file_name)
{
    std::ifstream infile(file_name, std::ios::binary);

    if (!infile.is_open())
    {
        return false;
    }

    infile.close();
    return true;
}

#endif //WEIGHT_COMMUNITY_FILE_CTRL_H

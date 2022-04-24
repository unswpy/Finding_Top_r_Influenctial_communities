//
// Created by sbian on 2020/7/27.
//

#include "serialize.h"

template <class tuple_type>
size_t get_tuple_size(const tuple_type& obj, int_<0>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - 1;
    return get_size(std::get<idx>(obj));
}

template <class tuple_type, size_t pos>
size_t get_tuple_size(const tuple_type& obj, int_<pos>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - pos - 1;
    size_t acc = get_size(std::get<idx>(obj));
    // recur
    return acc + get_tuple_size(obj, int_ < pos - 1 > ());
}

template <class tuple_type>
void serialize_tuple(const tuple_type& obj, StreamType::iterator& res, int_<0>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - 1;
    serializer(std::get<idx>(obj), res);
}

template <class tuple_type, size_t pos>
void serialize_tuple(const tuple_type& obj, StreamType::iterator& res, int_<pos>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - pos - 1;
    serializer(std::get<idx>(obj), res);
    // recur
    serialize_tuple(obj, res, int_ < pos - 1 > ());
}


template <class T>
T deserialize(StreamType::const_iterator& begin, const StreamType::const_iterator& end)
{
    return deserialize_helper<T>::apply(begin, end);
}

template <class T>
T deserialize(const StreamType& res)
{
    StreamType::const_iterator it = res.begin();
    return deserialize<T>(it, res.end());
}

template <class tuple_type>
void deserialize_tuple(tuple_type& obj,
                       StreamType::const_iterator& begin,
                       StreamType::const_iterator end, int_<0>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - 1;
    typedef typename std::tuple_element<idx, tuple_type>::type T;
    std::get<idx>(obj) = std::move(deserialize_helper<T>::apply(begin, end));
}

template <class tuple_type, size_t pos>
void deserialize_tuple(tuple_type& obj,
                       StreamType::const_iterator& begin,
                       StreamType::const_iterator end, int_<pos>)
{
    constexpr size_t idx = std::tuple_size<tuple_type>::value - pos - 1;
    typedef typename std::tuple_element<idx, tuple_type>::type T;
    std::get<idx>(obj) = std::move(deserialize_helper<T>::apply(begin, end));
    // recur
    deserialize_tuple(obj, begin, end, int_ < pos - 1 > ());
}

template <class T>
struct deserialize_helper
{
    static T apply(StreamType::const_iterator& begin,
                   StreamType::const_iterator end)
    {
        assert(begin + sizeof(T) <= end);
        T val;
        uint8_t* ptr = reinterpret_cast<uint8_t*>(&val);
        copy(begin, begin + sizeof(T), ptr);
        begin += sizeof(T);
        return val;
    }
};

template <class T>
struct deserialize_helper<std::vector<T>>
{
    static std::vector<T> apply(StreamType::const_iterator& begin,
                                StreamType::const_iterator end)
    {
        // retrieve the number of elements
        size_t size = deserialize_helper<size_t>::apply(begin, end);
        std::vector<T> vect(size);

        for (size_t i = 0; i < size; ++i)
        {
            vect[i] = std::move(deserialize_helper<T>::apply(begin, end));
        }

        return vect;
    }
};

template <>
struct deserialize_helper<std::string>
{
    static std::string apply(StreamType::const_iterator& begin,
                             StreamType::const_iterator end)
    {
        // retrieve the number of elements
        size_t size = deserialize_helper<size_t>::apply(begin, end);

        if (size == 0u) return std::string();

        std::string str(size, '\0');

        for (size_t i = 0; i < size; ++i)
        {
            str.at(i) = deserialize_helper<uint8_t>::apply(begin, end);
        }

        return str;
    }
};

template <class... T>
struct deserialize_helper<std::tuple<T...>>
{
    static std::tuple<T...> apply(StreamType::const_iterator& begin,
                                  StreamType::const_iterator end)
    {
        //return std::make_tuple(deserialize(begin,begin+sizeof(T),T())...);
        std::tuple<T...> ret;
        deserialize_tuple(ret, begin, end, int_ < sizeof...(T) - 1 > ());
        return ret;
    }
};

template <class... T>
struct serialize_helper<std::tuple<T...>>
{
    static void apply(const std::tuple<T...>& obj, StreamType::iterator& res)
    {
        serialize_tuple(obj, res, int_ < sizeof...(T) - 1 > ());
    }
};

template <>
struct serialize_helper<std::string>
{
    static void apply(const std::string& obj, StreamType::iterator& res)
    {
        // store the number of elements of this vector at the beginning
        serializer(obj.length(), res);

        for (const auto& cur : obj)
        {
            serializer(cur, res);
        }
    }
};

template <class T>
struct serialize_helper<std::vector<T>>
{
    static void apply(const std::vector<T>& obj, StreamType::iterator& res)
    {
        // store the number of elements of this vector at the beginning
        serializer(obj.size(), res);

        for (const auto& cur : obj)
        {
            serializer(cur, res);
        }
    }
};

template <class T>
struct serialize_helper
{
    static void apply(const T& obj, StreamType::iterator& res)
    {
        const uint8_t* ptr = reinterpret_cast<const uint8_t*>(&obj);
        copy(ptr, ptr + sizeof(T), res);
        res += sizeof(T);
    }
};

template <class T>
void serializer(const T& obj, StreamType::iterator& res)
{
    serialize_helper<T>::apply(obj, res);
}

template <class ...T>
struct get_size_helper<std::tuple<T...>>
{
    static size_t value(const std::tuple<T...>& obj)
    {
        return get_tuple_size(obj, int_ < sizeof...(T) - 1 > ());
    }
};

template <class T>
struct get_size_helper
{
    static size_t value(const T& obj)
    {
        return sizeof(T);
    }
};


template <class T>
struct get_size_helper<std::vector<T>>
{
    static size_t value(const std::vector<T>& obj)
    {
        return std::accumulate(obj.begin(), obj.end(), sizeof(size_t),
                               [](const size_t& acc, const T & cur)
                               {
                                   return acc + get_size(cur);
                               });
    }
};

template <>
struct get_size_helper<std::string>
{
    static size_t value(const std::string& obj)
    {
        return sizeof(size_t) + obj.length() * sizeof(uint8_t);
    }
};





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

template void load_serialized_graph(const std::string file_name, vector<vector<int>> & graph);
template void load_serialized_graph(const std::string file_name, vector<int> & graph);
template void load_serialized_graph(const std::string file_name, vector<double> & graph);

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


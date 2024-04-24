#pragma once

#include "JSONParser.h"

#include <Eigen/Eigen>
#include <cmath>
#include <functional>
#include <iostream>
#include <fstream>
#include <list>
#include <tetgen.h>
#include <triangle.h>
#include <vector>
#include <utility>
#include <unordered_map>

using std::vector;
using std::string;
using std::pair;
using std::unordered_map;

constexpr float pi = static_cast<float>(3.1415926535);

// Convenient typedefs
using colCoordMat = Eigen::Matrix2Xf;
using coord = std::pair<float, float>;

// TODO: make this hold the data directly
struct coord3D
{
    float x = 0;
    float y = 0;
    float z = 0;

    coord3D(float x, float y, float z) : x(x), y(y), z(z)
    {
    }

    coord3D operator+(const coord3D &other) const
    {
        return coord3D(x + other.x, y + other.y, z + other.z);
    }
};

struct Hash_Edge_to_Face
{
    unordered_map<int, unordered_map<int, pair<int, int>>> table;
    size_t double_edge_count = 0;

    pair<int, int> &at(const int &a, const int &b)
    {
        if (a > b)
        {
            if (!table[b].contains(a))
            {
                table[b][a] = {-1, -1};
            }
            return table[b][a];
        }
        else if (a < b)
        {
            if (!table[a].contains(b))
            {
                table[a][b] = {-1, -1};
            }
        }
        else
        {
            throw "invalid query";
        }

        return table[a][b];
    }
};

colCoordMat listToMatrix(std::list<coord> source);

colCoordMat vectorToMatrix(std::vector<coord> source);

std::vector<coord> matrixToVector(colCoordMat sourceMatrix);

inline void log2() {}

template <typename First, typename... Rest>
void log2(First first, Rest &&...rest)
{
    std::cout << std::forward<First>(first);
    log2(std::forward<Rest>(rest)...);
}

template <typename... Rest>
void log(Rest &&...rest)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(0);
    std::cout << "";
    log2(std::forward<Rest>(rest)...);
    std::cout << "" << std::endl;
}

// List Mapping Functions and overloads
template <typename T, typename G>
std::list<G> mapList(const std::list<T> &vec, const std::function<G(const T &)> &op)
{
    std::list<G> c;
    for (auto &item : vec)
    {
        c.push_back(op(item));
    }
    return c;
}

template <typename T, typename G>
std::list<G> operator<<(const std::list<T> &vec, const std::function<G(const T &)> &op) { return mapList(vec, op); }

// Vector Mapping Functions and overloads
template <typename T, typename G>
std::vector<G> mapVector(const std::vector<T> &vec, const std::function<G(const T &, size_t)> &op)
{
    std::vector<G> c;
    c.reserve(vec.size());
    for (size_t i = 0; i < vec.size(); i++)
    {
        c.push_back(op(vec[i], i));
    }
    return std::move(c);
}

template <typename T, typename G>
std::vector<G> mapVector(const std::vector<T> &vec, const std::function<G(const T &)> &op)
{
    std::function<G(const T &, size_t)> new_op([op](const T &val, size_t)
                                               { return op(val); });
    return mapVector(vec, new_op);
}

// mapThread and its overloads
template <typename A, typename B, typename C>
std::vector<C> mapThread(const std::vector<A> &vec1, const std::vector<B> &vec2,
                         const std::function<C(const A &, const B &, size_t)> &op)
{
    if (vec1.size() != vec2.size())
        throw "SIZE MISMATCH";

    std::vector<C> new_vector;
    new_vector.reserve(vec1.size());
    for (size_t i = 0; i < vec1.size(); i++)
    {
        new_vector.push_back(op(vec1[i], vec2[i], i));
    }
    return new_vector;
}

template <typename A, typename B, typename C>
std::vector<C> mapThread(const std::vector<A> &vec1, const std::vector<B> &vec2,
                         const std::function<C(const A &, const B &)> &op)
{
    return mapThread(vec1, vec2, std::function([&](const A &a, const B &b, size_t)
                                               { return op(a, b); }));
}

// To Cartesian Space
Eigen::Vector2f getPoint(const Eigen::Vector2f &coord, const Eigen::Vector2f &origin, const Eigen::Vector2f &v1,
                         const Eigen::Vector2f &v2);

// To Hex Space
Eigen::Vector2f getCoords(const Eigen::Vector2f &pt, const Eigen::Vector2f &origin, const Eigen::Vector2f &v1,
                          const Eigen::Vector2f &v2);

Eigen::Vector2i roundPtToCoord(const Eigen::Vector2f &pt, const Eigen::Vector2f &origin, const Eigen::Vector2f &v1,
                               const Eigen::Vector2f &v2);

template <typename T>
std::vector<T> filter(const std::vector<T> &vec, std::function<bool(const T &)> op)
{
    std::vector<T> res;
    for (const T &item : vec)
    {
        if (op(item))
        {
            res.push_back(item);
        }
    }
    return res;
}

template <typename T>
std::list<T> filter(const std::list<T> &vec, std::function<bool(const T &)> op)
{
    std::list<T> res;
    for (const T &item : vec)
    {
        if (op(item))
        {
            res.push_back(item);
        }
    }
    return res;
}

template <typename T>
std::vector<T> concat(std::vector<T> &&a, const std::vector<T> &b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
};

template <typename T>
std::vector<T> concat(const std::vector<T> &a, const std::vector<T> &b)
{
    std::vector<T> bothNeighbors;
    bothNeighbors.reserve(a.size() + b.size());
    bothNeighbors.insert(bothNeighbors.end(), a.begin(), a.end());
    bothNeighbors.insert(bothNeighbors.end(), b.begin(), b.end());
    return bothNeighbors;
};

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
    return concat(a, b);
}

// Runs function over range []
//  template <typename T>
//  std::vector<T> table(size_t length, std::function<T(size_t)> op)
//  {
//      return std::vector(length, 0)
//          << std::function<size_t(int, size_t)>([](int, size_t i) { return i; })
//          << op;
//  }

// Functions for ease of use with triangle and tetgen
inline Eigen::Vector2f getCornerVector(const triangulateio &obj, int corner_index)
{
    const auto x_index = corner_index * 2;
    return Eigen::Vector2d(obj.pointlist[x_index], obj.pointlist[x_index + 1]).cast<float>();
}

inline Eigen::Vector3f getCornerVector(const tetgenio &obj, int corner_index)
{
    const auto x_index = corner_index * 3;
    return Eigen::Vector3d(obj.pointlist[x_index], obj.pointlist[x_index + 1], obj.pointlist[x_index + 2]).cast<float>();
}

inline std::vector<int> getTriangleCornerIndices(const triangulateio &obj, int triangleIndex)
{
    const int *ptr = obj.trianglelist + static_cast<ptrdiff_t>(triangleIndex * obj.numberofcorners);
    return {ptr, ptr + obj.numberofcorners};
}

inline std::vector<int> getTetCornerIndices(const tetgenio &obj, int tetIndex)
{
    const int *ptr = obj.tetrahedronlist + static_cast<ptrdiff_t>(tetIndex * obj.numberofcorners);
    return {ptr, ptr + obj.numberofcorners};
}

inline Eigen::Matrix2Xf getTriangleMatrix(const triangulateio &obj, int triangleIndex)
{
    const auto corners = getTriangleCornerIndices(obj, triangleIndex);

    Eigen::Matrix2Xf ret = Eigen::Matrix2Xf::Zero(2, corners.size());

    for (size_t j = 0; j < corners.size(); j++)
    {
        ret.col(static_cast<Eigen::Index>(j)) = getCornerVector(obj, corners[j]);
    }
    return ret;
}

inline triangulateio triangulateMatrix(const Eigen::Matrix2Xf &mat)
{
    // TODO: Investigate flags further
    std::string flags =
        std::string("z") // Start arrays at 0
        // + std::string("V") //Verbose
        + std::string("Q") // Quiet
        ;

    const auto numPoints = mat.cols();
    auto points = new double[numPoints * 2];
    for (int i = 0; i < mat.cols(); i++)
    {
        points[i * 2] = mat.col(i)(0);
        points[i * 2 + 1] = mat.col(i)(1);
    }
    triangulateio in = {};
    in.pointlist = points;
    in.numberofpoints = static_cast<int>(numPoints);
    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;    // Might be able to use this to associate points with indices in the original matrix
    in.pointattributelist = nullptr; // Might be able to use this to associate points with indices in the original matrix

    triangulateio out = {};
    out.pointlist = nullptr;
    out.trianglelist = nullptr;
    out.pointmarkerlist = nullptr;
    out.pointattributelist = nullptr;
    out.trianglelist = nullptr;

    triangulate(&flags[0], &in, &out, nullptr);

    return out;
}

using json = nlohmann::json;

template <typename T, int G>
std::vector<T> eigenToVec(Eigen::Vector<T, G> mat) { return std::vector<T>(mat.data(), mat.data() + mat.rows()); }

template <typename T>
json toJson(std::vector<T> v)
{
    auto a = json::array();
    for (size_t i = 0; i < v.size(); i++)
    {
        a.push_back(v[i]);
    }
    return a;
}

inline json extractTriangleMathematicaMesh(const triangulateio &obj)
{
    std::vector<Eigen::Matrix<float, 2, 1, 0>> points;
    {
        points.reserve(obj.numberofpoints);
        for (int i = 0; i < obj.numberofpoints; ++i)
        {
            points.push_back(getCornerVector(obj, i));
        }
    }

    std::vector<std::vector<int>> triangles;
    {
        triangles.reserve(obj.numberoftriangles);
        for (int i = 0; i < obj.numberoftriangles; ++i)
        {
            triangles.push_back(getTriangleCornerIndices(obj, i));
        }
    }

    json pointJson = json::array();
    for (auto &pt : points)
    {
        pointJson.push_back(toJson(eigenToVec(pt)));
    }

    json triangleJson = json::array();
    for (auto &set : triangles)
    {
        triangleJson.push_back(toJson(set));
    }
    json ret = json::array({pointJson, triangleJson});
    return ret;
};

inline std::vector<std::vector<int>> tetgenToTetVector(const tetgenio &obj)
{
    std::vector<std::vector<int>> tets;
    tets.reserve(obj.numberoftetrahedra);
    for (int i = 0; i < obj.numberoftetrahedra; i++)
    {
        tets.push_back(getTetCornerIndices(obj, i));
    }
    return tets;
}

inline json extractTetMathematicaMesh(const tetgenio &obj)
{
    std::vector<Eigen::Vector3f> points;
    for (int i = 0; i < obj.numberofpoints; i++)
    {
        points.push_back(getCornerVector(obj, i));
    }

    std::vector<std::vector<int>> tets;
    for (int i = 0; i < obj.numberoftetrahedra; i++)
    {
        tets.push_back(getTetCornerIndices(obj, i));
    }

    auto pointJson = json(points);

    auto tetJson = json(tets);

    json ret = json::array({pointJson, tetJson});
    return ret;
};


// Return a subset of elements at the given indices
template <typename T>
std::vector<T> subset(const std::vector<T> &set, const std::vector<int> &indices)
{
    std::vector<T> new_vector;
    new_vector.reserve(indices.size());
    for (const auto &index : indices)
    {
        new_vector.push_back(set[index]);
    }
    return new_vector;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &input)
{
    unsigned int sum = 0;
    for (auto &layer : input)
    {
        sum += layer.size();
    }
    std::vector<T> result;
    result.reserve(sum);
    unsigned int i = 0;
    for (const auto &layer : input)
    {
        for (const auto &pt : layer)
        {
            result.push_back(pt);
        }
    }
    return result;
}

inline int exportObj(string path, vector<pair<vector<Eigen::Vector3f>, vector<vector<int>>>> data, vector<string> names)
{
    const int nums = data.size();

    for (int i = 0; i < nums; i++)
    {
        if (!std::filesystem::exists(path))
        {
            std::filesystem::create_directories(path);
        }
        std::ofstream target(path + names.at(i) + ".obj", std::ofstream::trunc);

        // export the points
        vector<Eigen::Vector3f> &vertices = data.at(i).first;
        Hash_Edge_to_Face hash;
        for (Eigen::Vector3f vertex : vertices)
        {
            std::stringstream ss;
            ss.setf(std::ios::fixed);
            ss.precision(16);
            ss << "v ";
            ss << vertex(0) << " ";
            ss << vertex(1) << " ";
            ss << vertex(2) << "\n";
            target << ss.rdbuf();
        }

        // export the faces
        vector<vector<int>> &faces = data.at(i).second;
        for (int j = 0; j < faces.size(); j++)
        {
            vector<int> &face = faces.at(j);
            std::stringstream ss;
            ss << "f ";
            for (int k = 0;k < face.size(); k++)
            {
                pair<int, int> edge;
                if (k == face.size() - 1)
                {
                    edge = {face.at(k), face.at(0)};
                }
                else
                {
                    edge = {face.at(k), face.at(k+1)};
                }

                int &target = edge.first < edge.second ? hash.at(edge.first, edge.second).first : hash.at(edge.first, edge.second).second;

                if (target != -1)
                {
                    log("Orientation error detected!");
                }
                else
                {
                    target = j;
                }
                ss << face.at(k) + 1 << " ";
            }
            ss << "\n";
            target << ss.rdbuf();
        }

        target.close();
    }

    return 0;
}
#pragma once

#include "Contour2D.h"
#include "UtilityFunctions.h"

#include <Eigen/Eigen>
#include <ranges>

using std::pair;
using std::tuple;
using std::vector;
using std::unordered_map;

#define LOADING_SIZE 10
#define DEBUG true
#ifndef __nullptr
#define __nullptr nullptr
#endif

// 1) Generate vertices, don't append, create flat array w/ fixed size
// 2) Remove 2d/3d vector to hashmaps (DONE)
// 3) Store indices of the other two vertices not on the edge (DONE)
// 4) Figure out how to reduce number of data structures

template <typename T>
class rigidVector
{
    T *data = __nullptr;
    size_t max = 0;
    size_t s = 0;

public:
    rigidVector()
    {
    }

    void reserve(size_t newMax)
    {
        max = newMax;
        delete[] data;
        data = new T[newMax];
    }

    T *start() { return data + 0; }
    T *end() { return data + s; }

    T &operator[](size_t index)
    {
        return data[index];
    }

    void push_back(const T &i)
    {
        data[s++] = i;
    }

    void push_back(T &&i)
    {
        data[s++] = std::move(i);
        if (s > max)
        {
            throw "Out of Range";
        }
    }
};

inline int orientation(const Eigen::Vector3f &u, const Eigen::Vector3f &v, const Eigen::Vector3f &w)
{
    float val = u.cross(v).dot(w);
    if (val < 0)
        return -1;
    if (val > 0)
        return 1;
    return 0;
}

inline int polyhedronOrientation(const Eigen::Vector3f &a, const Eigen::Vector3f &b, const Eigen::Vector3f &c, const Eigen::Vector3f &d)
{
    const Eigen::Vector3f u = b - a;
    const Eigen::Vector3f v = c - a;
    const Eigen::Vector3f w = d - a;

    Eigen::Matrix3f m;
    m << u, v, w;

    float val = m.determinant();
    if (val < 0)
        return -1;
    if (val > 0)
        return 1;
    return 0;
};

// Slower, but doesn't copy construct. Use when vectors are small.
inline vector<int> complementByReference(const vector<int> &source, const vector<int> &target)
{
    vector<int> accumulator;
    accumulator.reserve(source.size());
    for (const auto &item : source)
    {
        bool wasFound = false;
        for (const auto &t : target)
        {
            if (t == item)
                wasFound = true;
        }
        if (!wasFound)
            accumulator.push_back(item);
    }
//    std::ranges::sort(accumulator);
    return accumulator;
}

// O(nlog(n)+mlog(m) operation, but every usage n<=m<=4, so probably not a problem
inline vector<int> complement(vector<int> source, vector<int> target)
{
    std::ranges::sort(source);
    std::ranges::sort(target);

    vector<int> accumulator;
    accumulator.reserve(source.size());
    int i = 0, j = 0;
    while (i < source.size() && j < target.size())
    {
        auto &element = source[i];
        // If the element is ahead of the target comparison, there are no more matches with it
        if (element > target[j])
        {
            j++;
        }
        // If the element is behind the target comparison, we can safely use it
        else if (element < target[j])
        {
            accumulator.push_back(element);
            i++;
        }
        // If the element is equal, we know it cannot be in the set and we skip it
        else if (element == target[j])
        {
            i++;
        }
    }

    // Push back any remaining points
    for (; i < source.size(); i++)
    {
        accumulator.push_back(source[i]);
    }
    return accumulator;
}

inline pair<vector<int>, vector<int>> orderTets(const pair<int, int> &edge,
                                                const vector<vector<int>> &tets)
{
    pair<vector<int>, vector<int>> ret;
    vector<int> &orderedTets = ret.first;
    orderedTets.reserve(tets.size());

    // Which 2 other corners does each tet have?
    vector<vector<int>> cornersByTet; // This should have two non-edge vertices per tet
    {
        cornersByTet.reserve(tets.size());
        vector<int> edgeEndpoints = {edge.first, edge.second};
        for (const auto &tet : tets)
        {
            cornersByTet.push_back(complementByReference(tet, edgeEndpoints));
        }
    }

    // Stores which tet(s) each corner belongs to
    // Which tet, which index in tet
    unordered_map<int, vector<pair<int, int>>> tets_by_corner;
    {
        // For Each tet
        // O(2n) b/c each tet has exactly two stored corners
        size_t corner_index = 0;
        for (int tetIndex = 0; tetIndex < cornersByTet.size(); tetIndex++) // O(n)
        {
            // For each of the two corners in that tet
            for (corner_index = 0; corner_index < cornersByTet[tetIndex].size(); corner_index++) // O(2)
            {
                // Push back the corners that it touches
                const auto &corner = cornersByTet[tetIndex][corner_index];
                tets_by_corner[corner].reserve(10);
                tets_by_corner[corner].emplace_back(tetIndex, corner_index);
            }
        }
    }

    // Find the two endpoints (if they exist)
    auto &endpoints = ret.second;
    {
        // For each entry in tets_by_corner
        for (const auto &[corner, owningTets] : tets_by_corner)
        {
            if (owningTets.size() == 1) // If the corner only belongs to one tet
            {
                endpoints.push_back(corner); // Then save it as an endpoint

                if (endpoints.size() == 2)
                    break; // We know there will be 0 or 2
            }
        }
    }

    int endPoint;
    int nextCorner;
    // If there are no boundary spaces
    if (endpoints.empty())
    {
        // If there are no unique corners, just start with the first tet and go from there
        orderedTets.push_back(0);
        nextCorner = cornersByTet[0][1]; // start with the first tet's second point
        endPoint = cornersByTet[0][0];   // And end with its first point
    }
    // If there are boundary spaces
    else
    {
        // We know the associated set of tets will have exactly one entry
        const auto &temp = tets_by_corner[endpoints[0]][0]; // The tet associated with the first unique corner
        orderedTets.push_back(temp.first);
        nextCorner = cornersByTet[temp.first][1 - temp.second]; // Start with the other corner in that tet
        endPoint = endpoints[1];                                // End with the other unique corner
    }

    // Run through all the connected tets to order them correctly
    // Stop when we get to the endpoint (Which will be on the first tet we started at or on the other unconnected tet)
    while (nextCorner != endPoint) // O(n)
    {
        const auto &tet_entry = tets_by_corner[nextCorner];

        // If we have already looked at the first tet in the entry, then we know it is the second and vice versa
        // First is tet index, second is the corner's index in that tet (0 or 1 b/c we removed the two we are rotating around)
        auto [tet_index, corner_index] = tet_entry[0].first == orderedTets[orderedTets.size() - 1]
                                             ? tet_entry[1]
                                             : tet_entry[0];

        orderedTets.push_back(tet_index); // Push back the tet index
        nextCorner = cornersByTet[tet_index][static_cast<size_t>(1) - corner_index];
        // Save the other corner from that tet
    }

    return ret;
}

class Hash3d
{
public:
    // vector<vector<vector<int>>> table;
    unordered_map<int, unordered_map<int, unordered_map<int, int>>> table;
    // std::map<pair<size_t, size_t>, vector<int>> table;
    size_t size;

    Hash3d(size_t size) : size(size)
    {
        // table = vector(size, vector(size, vector(size,-1)));
        table.reserve(size);
    }

    int &at(size_t a, size_t b, size_t c)
    {
        // Order arguments
        if (b > a)
        {
            std::swap(a, b);
        }
        if (c > a)
        {
            std::swap(c, a);
        }
        if (c > b)
        {
            std::swap(b, c);
        }

        if (!table[a][b].contains(c))
        {
            table[a][b][c] = -1;
        }
        return table[a][b][c];
    }
};

struct Hash2d
{
    unordered_map<int, unordered_map<int, int>> table;
    // TODO: This is slow
    int &at(const int &a, const int &b)
    {
        if (a > b)
        {
            if (!table[b].contains(a))
            {
                table[b][a] = -1;
            }
            return table[b][a];
        }
        if (!table[a].contains(b))
        {
            table[a][b] = -1;
        }
        return table[a][b];
    }
};

inline tuple<vector<Eigen::Vector3f>, vector<vector<int>>, vector<pair<int, int>>> contourTetMultiDC(const vector<Eigen::Vector3f> &points_by_index,
                                                                                                     const vector<vector<int>> &tets_by_index,
                                                                                                     const vector<vector<float>> &vals_by_point_index)
{
    log("Contouring.");
    size_t number_of_materials = vals_by_point_index[0].size();
    size_t number_of_points = points_by_index.size();
    size_t number_of_tets = tets_by_index.size();

    vector<pair<int, int>> corner_combinations = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

    // get material index at each point
    vector<int> primary_material_by_point_index;
    {
        primary_material_by_point_index.reserve(number_of_points);
        for (size_t point_index = 0; point_index < points_by_index.size(); point_index++)
        {
            primary_material_by_point_index.push_back(getMaxPos(vals_by_point_index[point_index]));
        }
    }

    // create adj table from edges to faces

    Hash2d edge_index_by_endpoint_indices;
    edge_index_by_endpoint_indices.table.reserve(points_by_index.size());
    vector<pair<int, int>> edges_by_index;
    vector<vector<int>> tets_by_edge_index;
    {
        edges_by_index.reserve(number_of_tets * 6);
        tets_by_edge_index.reserve(number_of_tets * 6);
        const int display_fraction = number_of_tets % LOADING_SIZE == 0
                                         ? number_of_tets / LOADING_SIZE
                                         : (number_of_tets / LOADING_SIZE) + 1;

        // For each tet
        // O(n)
        log("Mapping Tets to Edges");
        for (size_t tet_index = 0; tet_index < number_of_tets; tet_index++)
        {
            const vector<int> &tet = tets_by_index[tet_index];

            // For each edge in that tet
            // O(6)
            for (const auto &corner_combination : corner_combinations)
            {
                // Get the edge's endpoints
                const int &firstEndpoint = tet[corner_combination.first];
                const int &secondEndpoint = tet[corner_combination.second];

                // See if we have looked at this edge before
                int &edge_index = edge_index_by_endpoint_indices.at(firstEndpoint, secondEndpoint);

                if (edge_index == -1) // If the edge does not exist
                {
                    // Add it to the various indexes we are building
                    edge_index = static_cast<int>(edges_by_index.size());

                    edges_by_index.emplace_back(firstEndpoint, secondEndpoint);

                    tets_by_edge_index.emplace_back();
                    tets_by_edge_index[edge_index].reserve(10);
                }

                // Track the tet part no matter what
                tets_by_edge_index[edge_index].push_back(tet_index);
            }
            if (tet_index % display_fraction == 0)
                log("  ", static_cast<float>(100 * tet_index / display_fraction) / LOADING_SIZE, "%");
        }
    }

    // create interpolation points_by_index, one per edge with material change
    unordered_map<int, Eigen::Vector3f> edgePoint_by_edge_index;
    {
        edgePoint_by_edge_index.reserve(edges_by_index.size());
        // O(n)
        for (int i = 0; i < edges_by_index.size(); i++)
        {
            pair<int, int> &edge = edges_by_index[i];
            const int &p_material_index = primary_material_by_point_index[edge.first];
            const int &q_material_index = primary_material_by_point_index[edge.second];

            // If there is a change between the p and q materials
            if (p_material_index != q_material_index)
            {
                const Eigen::Vector3f &p = points_by_index[edge.first];
                const Eigen::Vector3f &q = points_by_index[edge.second];
                const vector<float> &pvals = vals_by_point_index[edge.first];
                const vector<float> &qvals = vals_by_point_index[edge.second];

                // Calculate the difference between them
                edgePoint_by_edge_index[i] = interpEdge2Mat<3>(
                    p, q,
                    {pvals[p_material_index], pvals[q_material_index]},
                    {qvals[p_material_index], qvals[q_material_index]});
            }
        }
    }

    // create vertices, one per tet with material change
    vector<Eigen::Vector3f> vertex_by_index;
    vector<int> vertex_index_by_tet_index;
    {
        vector<Eigen::Vector3f> temp;
        vertex_by_index.reserve(number_of_tets * 6);
        vertex_index_by_tet_index.reserve(number_of_tets);
        for (const auto &tet : tets_by_index)
        {
            // If there is no material change within the tet
            if (primary_material_by_point_index[tet[0]] == primary_material_by_point_index[tet[1]] && primary_material_by_point_index[tet[0]] == primary_material_by_point_index[tet[2]] && primary_material_by_point_index[tet[0]] == primary_material_by_point_index[tet[3]])
            {
                vertex_index_by_tet_index.push_back(-1); // There is no central vertex
            }
            else
            {
                temp.reserve(corner_combinations.size());
                for (const auto &corner_combination : corner_combinations)
                {
                    const int &pointIndex = edge_index_by_endpoint_indices.at(
                        tet[corner_combination.first], tet[corner_combination.second]);

                    if (edgePoint_by_edge_index.contains(pointIndex))
                        temp.push_back(edgePoint_by_edge_index[pointIndex]);
                }
                vertex_by_index.push_back(getMassPoint<3>(temp));
                vertex_index_by_tet_index.push_back(vertex_by_index.size() - 1);
                temp.resize(0);
            }
        }
    }

    // Create Segments

    Hash3d bdFaceHash(edges_by_index.size());
    vector<vector<int>> segments_by_index;
    vector<pair<int, int>> segment_materials_by_edge_index;
    {
        log("Creating Contour Segements.");
        segments_by_index.reserve(edges_by_index.size());
        segment_materials_by_edge_index.reserve(edges_by_index.size());

        const vector<pair<int, int>> triangleEdges = {{0, 1}, {1, 2}, {2, 0}};

        const int print_constant = edges_by_index.size() % LOADING_SIZE == 0
                                       ? edges_by_index.size() / LOADING_SIZE
                                       : (edges_by_index.size() / LOADING_SIZE) + 1;

        // O(n), total O(n^2)
        // This function has the largest slowdown
        for (int edgeIndex = 0; edgeIndex < edges_by_index.size(); edgeIndex++)
        {
            const pair<int, int> &edge = edges_by_index[edgeIndex];

            // Material change
            if (primary_material_by_point_index[edge.first] != primary_material_by_point_index[edge.second])
            {
                segment_materials_by_edge_index.emplace_back(primary_material_by_point_index[edge.first], primary_material_by_point_index[edge.second]);

                const vector<int> &tets_around_edge_by_index = tets_by_edge_index[edgeIndex];
                const vector<vector<int>> tets_surrounding_edge_set = subset(tets_by_index, tets_around_edge_by_index);

                // O(n)
                const auto [ordered_tets_in_set, endpoints] = orderTets(edge, tets_surrounding_edge_set);

                // Each tet surrounding the edge
                const vector<int> tets_touching_edge_set = subset(tets_around_edge_by_index, ordered_tets_in_set);

                // Stores the vertices in order that make up the new segment
                vector<int> new_segment_vertices_set_ordered = subset(vertex_index_by_tet_index, tets_touching_edge_set);

                // Just in case cover and grow failed
                if (!endpoints.empty())
                {
                    // if there are boundary faces
                    vector<int> endpoint_vertices = {-1, -1};
                    // O(2)
                    for (int j = 0; j < endpoint_vertices.size(); j++)
                    {
                        vector<int> tri = {edge.first, edge.second, endpoints[j]};
                        std::ranges::sort(tri);

                        int &hashValue = bdFaceHash.at(tri[0], tri[1], tri[2]);
                        endpoint_vertices[j] = hashValue;

                        // If it doesn't exist
                        if (hashValue == -1)
                        {
                            // Add Face Vertex
                            vector<Eigen::Vector3f> massedPoints;
                            massedPoints.reserve(triangleEdges.size());
                            // O(3)
                            for (const auto &triangleEdgeEndpoints : triangleEdges)
                            {
                                int edge_index = edge_index_by_endpoint_indices.at(
                                    tri[triangleEdgeEndpoints.first], tri[triangleEdgeEndpoints.second]);
                                if (edgePoint_by_edge_index.contains(edge_index))
                                    massedPoints.push_back(edgePoint_by_edge_index[edge_index]);
                            }
                            vertex_by_index.push_back(getMassPoint<3>(massedPoints));
                            endpoint_vertices[j] = vertex_by_index.size() - 1;
                            hashValue = endpoint_vertices[j];
                        }
                    }

                    new_segment_vertices_set_ordered = concat(
                        concat({endpoint_vertices[0]}, new_segment_vertices_set_ordered), {endpoint_vertices[1]});
                }

                // Orient the polygon
                const vector<int> &tet1 = tets_by_index[tets_around_edge_by_index[ordered_tets_in_set[0]]];
                const vector<int> &tet2 = tets_by_index[tets_around_edge_by_index[ordered_tets_in_set[1]]];
                int p1 = complementByReference(tet1, tet2)[0];
                int p2 = complementByReference(tet1, {edge.first, edge.second, p1})[0];
                int p3 = edge.first;
                int p4 = edge.second;
                Eigen::Vector3f v1 = points_by_index.at(p1);
                Eigen::Vector3f v2 = points_by_index.at(p2);
                Eigen::Vector3f v3 = points_by_index.at(p3);
                Eigen::Vector3f v4 = points_by_index.at(p4);

                if (polyhedronOrientation(v1, v2, v3, v4) < 0)
                {
                    std::ranges::reverse(new_segment_vertices_set_ordered);
                }

                segments_by_index.push_back(new_segment_vertices_set_ordered);
            }

            if (edgeIndex % print_constant == 0)
                log("  ", static_cast<float>(100 * edgeIndex / print_constant) / LOADING_SIZE, "%");
        }
    }
    return {vertex_by_index, segments_by_index, segment_materials_by_edge_index};
}

inline void tetralizeMatrix(const Eigen::Matrix3Xf &pts, tetgenio &out)
{
    log("Tetrahedralizing Data.");
    tetgenio in;

    in.firstnumber = 0; // Arrays start at 0

    in.numberofpoints = static_cast<int>(pts.cols());
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            in.pointlist[3 * i + j] = static_cast<double>(pts.col(i)(j));
        }
    }
    tetgenbehavior *tmp = new tetgenbehavior();
    tmp->zeroindex = 1;
    tmp->quiet = 1;
    tetrahedralize(tmp, &in, &out);
}

inline Eigen::Vector3f getFaceNorm(const vector<Eigen::Vector3f> &pts)
{
    Eigen::Vector3f sum = {0, 0, 0};
    for (int i = 0; i < pts.size(); i++)
    {
        const auto offset = (i + 1) % pts.size();
        sum += pts[i].cross(pts[offset]);
    }
    return sum;
}

inline pair<vector<Eigen::Vector3f>, vector<vector<int>>> getContourByMat3D(
    const vector<Eigen::Vector3f> &verts,
    const vector<vector<int>> &segs,
    const vector<pair<int, int>> &segmats,
    int mat,
    float shrink)
{
    // Select segments by target material
    vector<int> first_index_match;
    vector<int> second_index_match;
    {
        for (int i = 0; i < segs.size(); i++)
        {
            if (segmats[i].first == mat)
                first_index_match.push_back(i);
            if (segmats[i].second == mat)
                second_index_match.push_back(i);
        }
    }

    vector<vector<int>> reversed_second_matches = subset(segs, second_index_match);
    for (vector<int> &seg : reversed_second_matches)
    {
        std::ranges::reverse(seg);
    }
    vector<vector<int>> matching_material_segments = concat(subset(segs, first_index_match), reversed_second_matches);
    // All these segments start with the target material

    vector<int> oldIndices_by_newIndex; // Indices of the vertices we care about
    {
        vector<bool> vertices_used_by_index(verts.size(), false);
        {
            // O(n)
            for (const auto &seg : matching_material_segments)
            {
                // O(2)
                for (const auto &ind : seg)
                {
                    vertices_used_by_index[ind] = true;
                }
            }
        }
        // O(n)
        for (int i = 0; i < vertices_used_by_index.size(); i++)
        {
            if (vertices_used_by_index[i])
            {
                oldIndices_by_newIndex.push_back(i);
            }
        }
    }

    vector<int> newIndices_by_oldIndex(verts.size(), -1); // Map the old indices to the new indices
    for (int i = 0; i < oldIndices_by_newIndex.size(); i++)
    {
        newIndices_by_oldIndex[oldIndices_by_newIndex[i]] = i;
    }
    vector<Eigen::Vector3f> new_vertices = subset(verts, oldIndices_by_newIndex);

    // New vertices (can be mapped)
    vector<vector<int>> revised_new_segments;
    {
        for (const auto &seg : matching_material_segments)
        {
            revised_new_segments.emplace_back();
            auto &newItem = revised_new_segments[revised_new_segments.size() - 1];
            newItem.reserve(seg.size());

            for (const int &pt : seg)
            {
                newItem.push_back(newIndices_by_oldIndex[pt]);
            }
        }
    }

    // shrink
    vector<Eigen::Vector3f> vertex_normals(new_vertices.size(), {0, 0, 0});
    for (const auto &seg : revised_new_segments)
    {
        vector<Eigen::Vector3f> points = subset(new_vertices, seg);
        Eigen::Vector3f nm = getFaceNorm(points);
        for (auto index : seg)
        {
            vertex_normals[index] += nm;
        }
    }

    for (size_t i = 0; i < vertex_normals.size(); i++)
    {
        const Eigen::Vector3f faceVector = vertex_normals[i];
        const Eigen::Vector3f normalized = faceVector.normalized();
        new_vertices[i] += shrink * normalized;
    }

    return {new_vertices, revised_new_segments};
}

inline vector<pair<vector<Eigen::Vector3f>, vector<vector<int>>>> getContourAllMats3D(
    const vector<Eigen::Vector3f> &verts,
    const vector<vector<int>> &segs,
    const vector<pair<int, int>> &segmats,
    const int &number_of_materials,
    const float &shrink)
{
    vector<pair<vector<Eigen::Vector3f>, vector<vector<int>>>> a;
    for (int i = 0; i < number_of_materials; i++)
    {
        a.push_back(getContourByMat3D(verts, segs, segmats, i, shrink));
    }
    return a;
}

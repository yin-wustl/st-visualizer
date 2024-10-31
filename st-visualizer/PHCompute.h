//
// Created by Yin Li on 9/25/24.
//

#ifndef ST_VISUALIZER_PHCOMPUTE_H
#define ST_VISUALIZER_PHCOMPUTE_H

#include <vector>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stack>
#include <set>

#include <Eigen/Eigen>

#include <phat/compute_persistence_pairs.h>
#include <phat/boundary_matrix.h>
#include <phat/representations/default_representations.h>

#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>
#include <phat/algorithms/swap_twist_reduction.h>
#include <phat/algorithms/exhaustive_compress_reduction.h>
#include <phat/algorithms/lazy_retrospective_reduction.h>

#include <phat/helpers/dualize.h>

using std::vector;
using std::tuple;
using std::unordered_map;
using std::unordered_set;
using std::stack;
using std::set;

struct SimplexHash {
    std::size_t operator()(const vector<int>& simplex) const
    {
        size_t result = 0;
        for (int vertex : simplex)
        {
            // FIXME: bad implementation, frequent collisions
            result = result ^ std::hash<int>{}(vertex);
        }

//        vector<int> sorted_simplex = simplex;
//        std::sort(sorted_simplex.begin(), sorted_simplex.end());
//        for(const int &i : sorted_simplex) {
//            result ^= i + 0x9e3779b9 + (result << 6) + (result >> 2);
//        }

        return result;
    }
};

struct SimplexEqual {
    bool operator()(const std::vector<int>& a, const std::vector<int>& b) const
    {
        if (a.size() != b.size()) return false;

        std::unordered_map<int, int> count;
        for (int num: a) count[num]++;
        for (int num: b)
        {
            if (--count[num] < 0) return false;
        }

        return true;
    }
};



vector<vector<int>> get_boundary(vector<int> &simplex)
{
    vector<vector<int>> result;
    int k = simplex.size() - 1;
    stack<pair<int, vector<int>>> stack;
    stack.push({0, {}});

    while (!stack.empty())
    {
        auto [start, current] = stack.top();
        stack.pop();

        if (current.size() == k)
        {
            result.push_back(current);
            continue;
        }

        for (int i = start; i < simplex.size(); i++)
        {
            vector<int> next = current;
            next.push_back(simplex[i]);
            stack.push({i + 1, next});
        }
    }

    return result;
}

vector<phat::index> get_boundary_idx(vector<int> &simplex, unordered_map<vector<int>, int, SimplexHash, SimplexEqual> &simplex_to_ind)
{
    vector<phat::index> result;
    if (simplex.size() == 1)
    {
        return result;
    }

    vector<vector<int>> boundaries = get_boundary(simplex);
    for (vector<int> &boundary_simplex : boundaries)
    {
        result.push_back(simplex_to_ind[boundary_simplex]);
    }

    return result;
}

float get_alpha(const vector<int> &simplex, const vector<float> &alphas)
{
    float largest = std::numeric_limits<float>::lowest();
    for (int vertex : simplex)
    {
        float alpha = alphas[vertex];
        if (alpha > largest)
        {
            largest = alpha;
        }
    }
    return largest;
}

bool compare_alpha(const tuple<vector<int>, float> &t_1, const tuple<vector<int>, float> &t_2)
{
    float alpha_1 = std::get<1>(t_1);
    float alpha_2 = std::get<1>(t_2);
    int dimension_1 = std::get<0>(t_1).size();
    int dimension_2 = std::get<0>(t_2).size();

    return std::tie(alpha_1, dimension_1) < std::tie(alpha_2, dimension_2);
}

void print_vector(const vector<int>& vec)
{
    vector<int> sorted = vec;
    std::sort(sorted.begin(), sorted.end());

    std::cout << "{ ";
    for (const int& val : sorted) {
        std::cout << val << " ";
    }
    std::cout << "}";
}

// for diagnostic purpose only
void write_tuples(const vector<tuple<vector<int>, float>>& tuples)
{
    string filename = "tmp.txt";
    std::ofstream outfile(filename, std::ios_base::out);

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const tuple<vector<int>, float> tuple : tuples)
    {
        vector<int> sorted = std::get<0>(tuple);
        std::sort(sorted.begin(), sorted.end());
        for (const int& val : sorted) {
            outfile << val << " ";
        }
        outfile << std::fixed << std::setprecision(5) << std::get<1>(tuple) << std::endl;
    }

    outfile.close();
}

void print_tuples(const vector<tuple<vector<int>, float>>& tuples)
{
    for (auto tuple : tuples)
    {
        const std::vector<int>& vec = std::get<0>(tuple);
        float f = std::get<1>(tuple);

        print_vector(vec);
        std::cout << std::fixed << std::setprecision(5) << f << std::endl;
    }
}

void print_filtration(const vector<tuple<vector<int>, float>> &filtration, bool print, bool write)
{
    if (print)
    {
        print_tuples(filtration);
    }
    if (write)
    {
        write_tuples(filtration);
    }
}

template <typename ColumnType>
void boundary_matrix_diagnose(const phat::boundary_matrix<ColumnType> &boundary_matrix, const vector<tuple<vector<int>, float>> &filtration, unordered_map<vector<int>, int, SimplexHash, SimplexEqual> simplex_to_ind)
{
    const int num_simplices = boundary_matrix.get_num_cols();
    for(int i = num_simplices - 1; i >= 0; i--)
    {
        vector<phat::index> col;
        boundary_matrix.get_col(i, col);

        vector<int> simplex = std::get<0>(filtration[i]);
        vector<phat::index> boundaries = get_boundary_idx(simplex, simplex_to_ind);

        std::vector<int> col_int(col.begin(), col.end());
        std::vector<int> boundaries_int(boundaries.begin(), boundaries.end());
        if (!SimplexEqual()(col_int, boundaries_int))
        {
            std::cout << "boundary error! " << std::endl;
            print_vector(col_int);
            print_vector(boundaries_int);
        }
    }
}


void compute_ph(const vector<vector<float>> &materials, const vector<vector<int>> &tets)
{
    int num_materials = materials[0].size();
    int num_points = materials.size();

    const vector<vector<int>> edge_combinations = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    const vector<vector<int>> triangle_combinations = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

    // FIXME: why didn't the simplex hash work?
    unordered_map<vector<int>, bool, SimplexHash, SimplexEqual> simplex_visited;
    vector<vector<int>> edges;
    vector<vector<int>> triangles;

    for (vector<int> tet : tets)
    {
        for (const vector<int> &combination : edge_combinations)
        {
//            vector<int> edge = {tet[combination[0]], tet[combination[1]]};
//            std::sort(edge.begin(), edge.end());
//            edges.insert(edge);

            const vector<int> edge = {tet[combination[0]], tet[combination[1]]};
            if (!simplex_visited[edge])
            {
                edges.push_back(edge);
                simplex_visited[edge] = true;
            }
        }
        for (const vector<int> &combination : triangle_combinations)
        {
//            vector<int> triangle = {tet[combination[0]], tet[combination[1]], tet[combination[2]]};
//            std::sort(triangle.begin(), triangle.end());
//            triangles.insert(triangle);

            const vector<int> triangle = {tet[combination[0]], tet[combination[1]], tet[combination[2]]};
            if (!simplex_visited[triangle])
            {
                triangles.push_back(triangle);
                simplex_visited[triangle] = true;
            }
        }
    }

    std::cout << "there are " << num_points << " points, " << edges.size() << " edges, " << triangles.size() << " triangles, " << tets.size() << " tets in filtratiion" << std::endl;

    // tuple: largest value, second largest value, index of largest value
    vector<tuple<float, float, int>> points_properties(num_points);

    for (int i = 0; i < num_points; i++)
    {
        const vector<float> &values = materials[i];

        float largest = std::numeric_limits<float>::lowest();
        float second_largest = std::numeric_limits<float>::lowest();
        int largest_index = -1;

        for (int j = 0; j < num_materials; j++)
        {
            if (values[j] > largest)
            {
                second_largest = largest;
                largest = values[j];
                largest_index = j;
            } else if (values[j] > second_largest && values[j] != largest)
            {
                second_largest = values[j];
            }
        }

        points_properties[i] = {largest, second_largest, largest_index};
    }

    // compute persistent homology for each material except the last one (no tissue)
    for (int material_idx = 0; material_idx < num_materials - 1; material_idx++)
    {
        vector<tuple<vector<int>, float>> filtration;
        vector<float> alphas(num_points);

        // determine the material specific alpha values for each point
        for (int i = 0; i < num_points; i++)
        {
            const vector<float> &values = materials[i];
            tuple<float, float, int> &point_properties = points_properties[i];

            float alpha;
            if (material_idx == std::get<2>(point_properties))
            {
                alpha = std::get<0>(point_properties) - std::get<1>(point_properties);
            } else
            {
                alpha = values[material_idx] - std::get<0>(point_properties);
            }
            alphas[i] = 1 - alpha;
        }

        // insert points to filtration
        for (int i = 0; i < num_points; i++)
        {
            filtration.push_back({{i}, alphas[i]});
        }

        // insert edges to filtration
        for (vector<int> edge : edges)
        {
            filtration.push_back({edge, get_alpha(edge, alphas)});
        }

        // insert triangles to filtration
        for (vector<int> triangle : triangles)
        {
            filtration.push_back({triangle, get_alpha(triangle, alphas)});
        }

        // insert tets to filtration
        for (vector<int> tet : tets)
        {
            filtration.push_back({tet, get_alpha(tet, alphas)});
        }

        int num_simplices = filtration.size();
        std::sort(filtration.begin(), filtration.end(), compare_alpha);
        print_filtration(filtration, false, true);

        unordered_map<vector<int>, int, SimplexHash, SimplexEqual> simplex_to_ind;
        for (int i = 0; i < filtration.size(); i++)
        {
            vector<int> &simplex = std::get<0>(filtration[i]);
            simplex_to_ind[simplex] = i;
        }

        // initialize boundary matrix
        phat::boundary_matrix<phat::sparse_pivot_column> boundary_matrix;
        boundary_matrix.set_num_cols(num_simplices);

        for (int i = 0; i < num_simplices; i++)
        {
            tuple<vector<int>, float> &tuple = filtration[i];
            vector<int> &simplex = std::get<0>(tuple);
            int dimension = simplex.size() - 1;
            vector<phat::index> boundary_idx = get_boundary_idx(simplex, simplex_to_ind);
            std::sort(boundary_idx.begin(), boundary_idx.end());
            boundary_matrix.set_col(i, boundary_idx);
            boundary_matrix.set_dim(i, dimension);
        }

        boundary_matrix_diagnose(boundary_matrix, filtration, simplex_to_ind);

        phat::persistence_pairs pairs;
        phat::compute_persistence_pairs< phat::standard_reduction >( pairs, boundary_matrix );
        pairs.sort();

        // print the pairs:
        std::cout << std::endl;
        std::cout << "There are " << pairs.get_num_pairs() << " persistence pairs: " << std::endl;
        for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
            std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
    }


}

#endif //ST_VISUALIZER_PHCOMPUTE_H

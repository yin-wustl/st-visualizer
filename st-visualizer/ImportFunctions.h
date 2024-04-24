#pragma once

#include "GrowAndCover.h"
#include "UtilityFunctions.h"

#include <algorithm>
#include <chrono>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

using std::pair;
using std::string;
using std::vector;
using std::function;

struct tsv_return_type
{
    int num_points;
	vector<string> names;
    vector<string> clusterNames;
	vector<Eigen::Matrix3Xf> slices;
	vector<vector<vector<float>>> clusters;
	vector<vector<vector<float>>> values;
};

vector<pair<vector<coord>, vector<coord>>> importAlignments(const string &alignment_file);

tsv_return_type loadTsv(const string &file_name,
                        const vector<string> &slice_names,
                        unsigned int slice_index,
                        unsigned int tissue_index,
                        pair<unsigned, unsigned> row_col_indices,
                        unsigned int cluster_ind,
                        const vector<unsigned> &feature_indices,
                        unsigned int z_distance,
                        vector<pair<vector<coord>, vector<coord>>> source_targets);

function<vector<coord>(vector<coord>)> getTransSVD(const vector<coord> &source, const vector<coord> &target);

// this helper function produces an array of length n, which is all zero except a 1 in the i-th spot.
vector<float> getClusterArray(size_t length, size_t i);

Eigen::Matrix2f getSVDRotation(colCoordMat source_matrix, colCoordMat target_matrix);

colCoordMat translateToZeroCentroid(colCoordMat source_matrix);

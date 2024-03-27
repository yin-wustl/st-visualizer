
#include "ImportFunctions.h"
#include "GrowAndCover.h"
#include "UtilityFunctions.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

using std::pair;
using std::string;
using std::vector;

constexpr int wid_buffer = 2;
constexpr int num_ransac = 20;

// From https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
vector<string> splitString(const string &s, const string &delimiter)
{
    size_t pos_start = 0, pos_end;
    const size_t delimiter_length = delimiter.length();
    vector<string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != string::npos)
    {
        string token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delimiter_length;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

vector<coord> zip(const vector<float> &top, const vector<float> &bottom)
{
    if (top.size() != bottom.size())
    {
        throw std::range_error("Size mismatch");
    }
    vector<coord> res;
    res.reserve(top.size());
    for (size_t i = 0; i < top.size(); i++)
    {
        res.emplace_back(coord(top[i], bottom[i]));
    }
    return res;
}

// Any impossible fields default to inf
vector<float> convertRowToFloat(const vector<string> &input)
{
    vector<float> res;
    res.reserve(input.size());
    for (size_t i = 0; i < input.size(); i++)
    {
        try
        {
            res.emplace_back(std::stof(input[i]));
        }
        catch (...)
        {
            res.emplace_back(std::numeric_limits<float>::infinity());
        }
    }
    return res;
}

// Import data from alignment json file
// TODO: handle invalid alignment file input
// Split into cells excluding the first row, assuming the alignment file has more than 1 lines
vector<pair<vector<coord>, vector<coord>>> importAlignments(const string &alignment_file)
{
    std::ifstream aFile(alignment_file);
    vector<vector<string>> csvCells;
    vector<string> csvFeatures;
    if (aFile.is_open())
    {
        string line;
        std::getline(aFile, line);
        csvFeatures = splitString(line, ",");
        while (std::getline(aFile, line))
        {
            vector<string> row = splitString(line, ",");
            row.erase(row.begin(), row.begin() + 3);
            csvCells.push_back(row);
        }
        aFile.close();
    }
    else
    {
        throw "ALIGNMENT FILE NOT FOUND";
    }

    // Transform cells into float
    const vector<vector<float>> transformCells = mapVector(csvCells, std::function(convertRowToFloat));

    // Now we have n columns and m*4 rows of ints
    // This part transforms that into an m by 2 by n structure of coords

    // Transpose pairs
    vector<vector<coord>> coordinateSet;
    for (size_t i = 0; i < transformCells.size(); i += 2)
    {
        coordinateSet.push_back(zip(transformCells[i], transformCells[i + 1]));
    }

    vector<pair<vector<coord>, vector<coord>>> finalSet;
    // line up in pairs
    for (size_t i = 0; i < coordinateSet.size(); i += 2)
    {
        finalSet.emplace_back(coordinateSet[i], coordinateSet[i + 1]);
    }

    return finalSet;
}

tsv_return_type loadTsv(const string &file_name,
                        const vector<string> &slice_names,
                        unsigned int slice_index,
                        unsigned int tissue_index,
                        pair<unsigned, unsigned> row_col_indices,
                        unsigned int cluster_ind,
                        const vector<unsigned> &feature_indices,
                        unsigned int z_distance,
                        vector<pair<vector<coord>, vector<coord>>> source_targets)
{
    // Import raw file data
    std::ifstream aFile(file_name);
    log("Loading TSV.");
    vector<vector<string>> rawData;
    if (aFile.is_open())
    {
        string line;
        while (std::getline(aFile, line))
        {
            rawData.push_back(splitString(line, "\t"));
        }
        aFile.close();
    }
    else
    {
        throw "TSV FILE NOT FOUND";
    }

    log("Loading names.");
    vector<string> names = mapVector(feature_indices, std::function([&rawData](const unsigned &index, size_t)
                                                                    { return rawData.front()[index]; }));
    names.emplace_back("No Tissue");

    rawData.erase(rawData.begin());
    vector<vector<string>> &tab = rawData;

    // Get the max cluster number
    unsigned int max = 0;
    for (const auto &row : tab)
    {
        {
            try
            {
                int n = std::stoi(row[cluster_ind]);
                if (n > max)
                    max = n;
            }
            catch (...)
            {
            }
        }
    }

    // Add 1 for no tissue
    unsigned int newClusters = max + 1;
    size_t newFeatures = feature_indices.size();

    // Generate cluster name array
    vector<int> clusterIndices(newClusters, 0);
    std::iota(clusterIndices.begin() + 1, clusterIndices.end(), 1);
    vector<string> clusterNames = mapVector(clusterIndices, std::function([](const int &index)
                                                                          { return std::to_string(index); }));
    clusterNames.emplace_back("No Tissue");

    // Extract the relevant records
    // Records should hold slices of data rows that match the right name of the slice and is a tissue sample
    vector<vector<vector<string>>> sliced_records = mapVector(slice_names, std::function(
                                                                               [&tab, &slice_index, &tissue_index](const string &name, size_t)
                                                                               {
                                                                                   // Filter out the rows that are not tissues
                                                                                   return filter(tab, std::function(
                                                                                                          [&tab, &slice_index, &tissue_index, &name](const vector<string> &row)
                                                                                                          {
                                                                                                              return row[slice_index] == name && row[tissue_index] == "1";
                                                                                                          }));
                                                                               }));

    // This is the point where parallel vectors are created

    // Pull out the corresponding xy data from the records, adjusted to be a unified coordinate system
    pair<unsigned int, unsigned int> xy_indices(row_col_indices.second, row_col_indices.first);
    vector<Eigen::Matrix2Xf> slices = mapVector(sliced_records, std::function(
                                                                    [&xy_indices, &source_targets](const vector<vector<string>> &record, size_t i)
                                                                    {
                                                                        const vector<coord> raw_slice_coordinates_vector = mapVector(
                                                                            record, std::function(
                                                                                        [&xy_indices](const vector<string> &row, size_t)
                                                                                        {
                                                                                            return pair(
                                                                                                std::stof(row[xy_indices.first]),
                                                                                                std::stof(row[xy_indices.second]));
                                                                                        }));
                                                                        // If it's the first slice, no adjustment necessary
                                                                        if (i == 0)
                                                                        {
                                                                            return vectorToMatrix(raw_slice_coordinates_vector);
                                                                        }

                                                                        // Base the remaining slices coordinate adjustment off of the previous one
                                                                        // Basically the other slices get transformed so that they match the coordinates of the first
                                                                        // All the transforms are independent.
                                                                        const std::function transform = getTransSVD(
                                                                            source_targets[i - 1].first, source_targets[i - 1].second);
                                                                        return vectorToMatrix(transform(raw_slice_coordinates_vector));
                                                                    }));

    // Convert the clusters into an array of 1/0 based on the cluster index (one hot indexing)
    // Clusters are represented as vectors with all values zero, except a single 1 in the ith place where i is the cluster it belongs to
    vector<vector<vector<float>>> original_clusters = mapVector(sliced_records, std::function(
                                                                                    [&](const vector<vector<string>> &record, size_t)
                                                                                    {
                                                                                        return mapVector(record, std::function(
                                                                                                                     [&](const vector<string> &row, size_t)
                                                                                                                     {
                                                                                                                         // If the data point doesn't contain tissue, give it a "none" value, otherwise use the cluster value it is a part of
                                                                                                                         // FIXME: why do we need +1 here?
                                                                                                                         return getClusterArray(newClusters + 1, row[tissue_index] == "0" ? newClusters : std::stoi(row[cluster_ind]));
                                                                                                                     }));
                                                                                    }));

    log("Parsing Values.");
    // Values are represented as arrays too, but the values are not just 1 or 0, but are all floats (except for the last index)
    vector<vector<vector<float>>> values = mapVector(sliced_records, std::function([&](const vector<vector<string>> &record, size_t)
                                                                                   { return mapVector(record, std::function([&](const vector<string> &row, size_t)
                                                                                                                            {
            if(row[tissue_index] == "0")
            {
                return getClusterArray(newFeatures + 1, newFeatures);
            }

            vector<float> a = mapVector(feature_indices, std::function([&](const unsigned& index, size_t)
            {
                try {
                    return std::stof(row[index]);
                } catch (std::exception &e) {
                    return 0.0f; 
                } 
            }));
            a.emplace_back(0);
            return a; })); }));

    log("Growing Slices.");
    // Add buffer to each slice and grow and cover neighboring slices
    vector<Eigen::Matrix2Xf> new_slice_data = mapVector(slices, std::function([&](const Eigen::Matrix2Xf &, size_t i)
                                                                              {
        log("  ", i + 1, "/", slices.size(), " slices");
        // If it's the first slice
        if(i == 0)
        {
            return growAndCover(slices[i], slices[i + 1], wid_buffer, num_ransac);
        }

        // If it's the last slice
        else if(i == slices.size() - 1)
        {
            return growAndCover(slices[i], slices[i - 1], wid_buffer, num_ransac);
        }

        // All other slices, first combine points from the previous and next slices together
        Eigen::Matrix2Xf top_and_bottom_slice(2, slices[i + 1].cols() + slices[i - 1].cols());
        top_and_bottom_slice << slices[i + 1], slices[i - 1];
        return growAndCover(slices[i], top_and_bottom_slice, wid_buffer, num_ransac); }));

    vector<Eigen::Matrix3Xf> slices3d = mapThread(
        new_slice_data, slices, std::function([z_distance](const Eigen::Matrix2Xf &new_slice, const Eigen::Matrix2Xf &old_slice, size_t i)
                                              {
                Eigen::Matrix2Xf layer2d(
                    2, new_slice.cols() + old_slice.cols());
                layer2d << old_slice, new_slice;
                //The old and new points on that layer (in that order)
                Eigen::Matrix3Xf layer3d(3, layer2d.cols());
                for(int j = 0; j < layer2d.cols(); j++)
                {
                    layer3d.col(j)(0) = layer2d.col(j)(0);
                    layer3d.col(j)(1) = layer2d.col(j)(1);
                    layer3d.col(j)(2) = static_cast<float>(z_distance * i);
                }

                return layer3d; }));

    // Associating data with the new points on their respective slices
    vector<vector<vector<float>>> grown_clusters = mapThread(
        new_slice_data, original_clusters, std::function([&](const Eigen::Matrix2Xf &new_coordinates, const vector<vector<float>> &old_clusters)
                                                         {
                //Associate null data with the areas that we grew to and append them to the list of locations
                const vector new_clusters(new_coordinates.cols(), getClusterArray(newClusters + 1, newClusters));
                return concat(old_clusters, new_clusters); }));

    vector<vector<vector<float>>> grown_values = mapThread(
        new_slice_data, values, std::function([&](const Eigen::Matrix2Xf &new_coordinates, const vector<vector<float>> &old_clusters)
                                              {
                //Same process as the previous step, but this time for vals instead of clusters
                const vector tempVector(
                    new_coordinates.cols(),
                    getClusterArray(newFeatures + 1, newFeatures));
                return concat(old_clusters, tempVector); }));

    log("Adding bounding slices.");

    // Add empty top and bottom slices
    {
        auto &topSlice = slices3d[slices3d.size() - 1];
        Eigen::Matrix3Xf top(3, topSlice.cols());
        for (int i = 0; i < topSlice.cols(); i++)
        {
            top.col(i) = topSlice.col(i) + Eigen::Vector3f({0, 0, static_cast<float>(z_distance)});
        }
        slices3d.push_back(top);

        auto &bottomSlice = slices3d[0];
        Eigen::Matrix3Xf bottom(3, bottomSlice.cols());
        for (int i = 0; i < bottomSlice.cols(); i++)
        {
            bottom.col(i) = bottomSlice.col(i) - Eigen::Vector3f({0, 0, static_cast<float>(z_distance)});
        }
        slices3d.insert(slices3d.begin(), bottom);
        // TODO: Inefficient O(n) call. Consider changing in all three locations.
    }
    {
        auto &topSlice = grown_clusters[grown_clusters.size() - 1];
        vector top(topSlice.size(), getClusterArray(newClusters + 1, newClusters));
        grown_clusters.push_back(top);

        auto &bottomSlice = grown_clusters[0];
        vector bottom(bottomSlice.size(), getClusterArray(newClusters + 1, newClusters));
        grown_clusters.insert(grown_clusters.begin(), bottom);
    }

    {
        auto &topSlice = grown_values[grown_values.size() - 1];
        vector top(topSlice.size(), getClusterArray(newFeatures + 1, newFeatures));
        grown_values.push_back(top);

        auto &bottomSlice = grown_values[0];
        vector bottom(bottomSlice.size(), getClusterArray(newFeatures + 1, newFeatures));
        grown_values.insert(grown_values.begin(), bottom);
    }

    // TODO: std::move for speed reasons
    tsv_return_type ret;
    ret.names = names;
    ret.clusterNames = clusterNames;
    ret.slices = slices3d;
    ret.values = grown_values;
    ret.clusters = grown_clusters;
    log("TSV Import Complete.");
    return ret;
}

Eigen::Vector2f getCentroid(colCoordMat sourceMatrix)
{
    return sourceMatrix.rowwise().mean();
}

colCoordMat translateToZeroCentroid(colCoordMat source_matrix)
{
    // Get the average of each row
    const auto centroid = getCentroid(source_matrix);

    // Subtract the centroid from each column;
    return source_matrix.colwise() - centroid;
}

Eigen::Matrix2f getSVDRotation(colCoordMat source_matrix, colCoordMat target_matrix)
{
    // Row 0 is x, row 1 is y

    //(* getting the centroid *)
    const colCoordMat zeroSource = translateToZeroCentroid(source_matrix);
    colCoordMat zeroTarget = translateToZeroCentroid(target_matrix);

    const Eigen::Matrix2f mat = zeroSource * zeroTarget.transpose();

    //(* SVD decomposition *)
    const Eigen::JacobiSVD<colCoordMat> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //(* obtaining the rotation *)
    Eigen::Matrix2f r = (svd.matrixU() * svd.matrixV().transpose()).transpose(); // This is definitely a rotation matrix

    return r;
}

std::function<vector<coord>(vector<coord>)> getTransSVD(const vector<coord> &source,
                                                        const vector<coord> &target)
{
    const auto sourceMatrix = vectorToMatrix(source);
    const auto targetMatrix = vectorToMatrix(target);

    // Convert to matrices
    const auto r = getSVDRotation(sourceMatrix, targetMatrix);
    const Eigen::Rotation2D<float> rotation(r);

    const Eigen::Vector2f sourceCentroid = getCentroid(sourceMatrix);
    const Eigen::Vector2f targetCentroid = getCentroid(targetMatrix);
    const Eigen::Vector2f targetCentroidTransform = targetCentroid - sourceCentroid;

    const Eigen::Translation2f netTranslation(targetCentroidTransform);
    const Eigen::Translation2f toZero(-1 * sourceCentroid);
    const Eigen::Translation2f fromZero(sourceCentroid);
    Eigen::Transform<float, 2, Eigen::Affine> finalTransform = netTranslation * fromZero * rotation * toZero;
    // Translate after rotate

    //(* transform *) Creating the function
    return {
        [finalTransform](vector<coord> points)
        {
            return matrixToVector(finalTransform * vectorToMatrix(std::move(points)));
        }};
}

vector<float> getClusterArray(size_t length, size_t i)
{
    vector<float> ret(length, 0);
    if (i < length)
    {
        ret[i] = 1;
    }
    return ret;
}

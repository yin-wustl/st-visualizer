#pragma once

#include "UtilityFunctions.h"

#include <map>
#include <vector>
#include <chrono>

using std::pair;
using std::tuple;
using std::vector;



// orientation (not sure if this should be i or f)
// Returns 1 if left, 0 if none, -1 if right (right hand rule)
int orientation(const Eigen::Vector2f &a, const Eigen::Vector2f &b);

int getMaxPos(const vector<float> &material_values);

// getMassPoint
// Might need work based on input type
template<unsigned int N>
Eigen::Matrix<float, N, 1> getMassPoint(const vector<Eigen::Matrix<float, N, 1>> &a)
{
    Eigen::Matrix<float, N, 1> temp = Eigen::Matrix<float, N, 1>::Zero(N);
    for (const Eigen::Matrix<float, N, 1> &v: a)
    {
        temp += v;
    }
    return temp / a.size();
}

// InterpEdge2Mat
// Does a linear interpolation of two materials to produce the midpoint location where the "cross over"
template<unsigned int N>
Eigen::Matrix<float, N, 1> interpEdge2Mat(const Eigen::Matrix<float, N, 1> &p,
                                          const Eigen::Matrix<float, N, 1> &q,
                                          const pair<float, float> &p_vals,
                                          const pair<float, float> &q_vals)
{
    const float top = std::abs(p_vals.first - p_vals.second);
    const float bottom = std::abs(p_vals.first - p_vals.second) + std::abs(q_vals.first - q_vals.second);
    const float t = top / bottom;
    return p * (1 - t) + q * t;
}

// contourTriMultiDC
struct contourTriMultiDCStruct
{
    vector<Eigen::Vector2f> verts;
    vector<pair<int, int>> segs;
    vector<pair<int, int>> segMats;
    vector<Eigen::Vector2f> fillVerts;
    vector<vector<int>> fillTris;
    vector<int> fillMats;
};

contourTriMultiDCStruct contourTriMultiDC(const Eigen::Matrix2Xf &pointIndexToPoint,
                                          const vector<vector<int>> &triangleIndexToCornerIndices,
                                          const vector<vector<float>> &pointIndexToMaterialValues);

inline Eigen::Vector2f perp(const Eigen::Vector2f &a)
{ return {-1 * a[1], a[0]}; }

// getContourByMat2D returns new vertices, new segments
inline pair<vector<Eigen::Vector2f>, vector<pair<int, int>>> getContourByMat2D(
        const vector<Eigen::Vector2f> &verts,
        const vector<pair<int, int>> &segs,
        const vector<pair<int, int>> &segmats,
        const int &mat,
        const float &shrink)
{
    // select segments by mat
    vector<pair<int, int>> newSegments;
    {
        newSegments.reserve(segmats.size() * 2);
        for (size_t i = 0; i < segmats.size(); i++)
        {
            if (segmats[i].first == mat)
            {
                const auto &segment = segs[i];
                newSegments.emplace_back(segment.first, segment.second);
            }
        }

        for (size_t i = 0; i < segmats.size(); i++)
        {
            if (segmats[i].second == mat)
            {
                const auto &segment = segs[i];
                newSegments.emplace_back(segment.second, segment.first);
            }
        }
    }

    /*prune unused vertices*/
    vector<bool> vertsUsed(verts.size(), false);
    for (auto &seg: newSegments)
    {
        vertsUsed[seg.first] = true;
        vertsUsed[seg.second] = true;
    }

    vector<int> nvertInds;
    for (int i = 0; i < vertsUsed.size(); i++)
    {
        if (vertsUsed[i])
        {
            nvertInds.push_back(i);
        }
    }

    vector<size_t> vertNewInds(verts.size(), 0);
    for (size_t i = 0; i < nvertInds.size(); i++)
    {
        vertNewInds[nvertInds[i]] = i;
    }

    vector<Eigen::Vector2f> nverts = subset(verts, nvertInds);
    vector<pair<int, int>> adjusted_nsegs;
    {
        adjusted_nsegs.reserve(newSegments.size());
        for (auto &seg: newSegments)
        {
            adjusted_nsegs.emplace_back(vertNewInds[seg.first], vertNewInds[seg.second]);
        }
    }

    /*shrink*/
    vector<pair<float, float>> vertnorms(nverts.size(), {0, 0});
    for (auto &seg: adjusted_nsegs)
    {
        Eigen::Vector2f nm = perp(nverts[seg.second] - nverts[seg.first]);
        vertnorms[seg.first].first += nm(0);
        vertnorms[seg.first].second += nm(1);

        vertnorms[seg.second].first += nm(0);
        vertnorms[seg.second].second += nm(1);
    }

    for (int i = 0; i < nverts.size(); i++)
    {
        Eigen::Vector2f &pt = nverts[i];
        pt = pt + shrink * Eigen::Vector2f(vertnorms[i].first, vertnorms[i].second).normalized();
    }

    return {nverts, adjusted_nsegs};
}

// getContourAllMats2D each index is new vertices, new segments
inline vector<pair<vector<Eigen::Matrix<float, 2, 1, 0>>, vector<pair<int, int>>>>
getContourAllMats2D(const vector<Eigen::Vector2f> &verts, const vector<pair<int, int>> &segs,
                    const vector<pair<int, int>> &segmats, const int &nmat, const float &shrink)
{
    vector<pair<vector<Eigen::Matrix<float, 2, 1, 0>>, vector<pair<int, int>>>> ret;
    {
        ret.reserve(nmat);
        for (int i = 0; i < nmat; ++i)
        {
            ret.push_back(getContourByMat2D(verts, segs, segmats, i, shrink));
        }
    }
    return ret;
}

pair<
        vector<
                pair<
                        vector<Eigen::Vector3f>,
                        vector<pair<int, int>>>>,
        tuple<
                vector<Eigen::Vector3f>,
                vector<vector<int>>,
                vector<int>
        >>
getSectionContours(const Eigen::Matrix3Xf &pts, const vector<vector<float>> &vals, float shrink);

pair<vector<vector<pair<vector<Eigen::Vector3f>, vector<pair<int, int>>>>>,
        vector<tuple<vector<Eigen::Vector3f>, vector<vector<int>>, vector<int>>>>
getSectionContoursAll(vector<Eigen::Matrix3Xf> sections,
                      vector<vector<vector<float>>> vals,
                      float shrink);

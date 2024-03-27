
#include "Contour2D.h"

using std::vector;
using std::pair;

int orientation(const Eigen::Vector2f &a, const Eigen::Vector2f &b)
{
    Eigen::Vector3f intermediateA;
    intermediateA(0) = a(0);
    intermediateA(1) = a(1);
    intermediateA(2) = 0;

    Eigen::Vector3f intermediateB;
    intermediateB(0) = b(0);
    intermediateB(1) = b(1);
    intermediateB(2) = 0;
    const float result = intermediateA.cross(intermediateB).eval()(2);
    return result == 0.0f ? 0 : (result < 0.0f ? -1 : 1);
}

int getMaxPos(const vector<float> &material_values)
{
    int max_index = 0;
    for (size_t i = 0; i < material_values.size(); i++)
    {
        if (material_values[max_index] < material_values[i])
        {
            max_index = static_cast<int>(i);
        }
    }

    return max_index;
}

contourTriMultiDCStruct contourTriMultiDC(const Eigen::Matrix2Xf &pointIndexToPoint,
                                          const vector<vector<int>> &triangleIndexToCornerIndices,
                                          const vector<vector<float>> &pointIndexToMaterialValues)
{
    // Step 1: Set up a structure to define geometry
    const size_t numberOfTriangles = triangleIndexToCornerIndices.size();

    // This is a list of the combinations of edges on a triangle by corner
    const vector<pair<int, int>> triangle_edges = {{0, 1}, {1, 2}, {2, 0}};

    // Primary Material index at each point
    const vector<int> primaryMaterialIndexByPointIndex = mapVector(pointIndexToMaterialValues, std::function(getMaxPos));

    const size_t number_of_points = pointIndexToPoint.cols();
    vector<vector<int>> endpointIndicesToEdgeIndex(number_of_points, vector<int>(number_of_points, -1));
    
    vector<pair<int, int>> edgeIndexToEndpointIndices;
    edgeIndexToEndpointIndices.reserve(triangle_edges.size());
    // Stores which face index an edge belongs to
    vector<vector<int>> edgeIndexToFaceIndices;
    edgeIndexToFaceIndices.reserve(triangle_edges.size());

    // Make connections between edges/faces, edges/endpoints, and endpoints/edges
    for (int triangleSide = 0; triangleSide < 3; triangleSide++)
    {
        for (int faceIndex = 0; faceIndex < numberOfTriangles; faceIndex++)
        {
            int firstEndpointIndex = triangleIndexToCornerIndices[faceIndex][triangle_edges[triangleSide].first];
            int secondEndpointIndex = triangleIndexToCornerIndices[faceIndex][triangle_edges[triangleSide].second];

            if (endpointIndicesToEdgeIndex[firstEndpointIndex][secondEndpointIndex] == -1) // If the edge doesn't already exist
            {
                const int num_of_edges = edgeIndexToEndpointIndices.size();

                const pair<int, int> &cornerPair = triangle_edges[triangleSide];
                edgeIndexToEndpointIndices.push_back({triangleIndexToCornerIndices[faceIndex][cornerPair.first],
                                                      triangleIndexToCornerIndices[faceIndex][cornerPair.second]});
                edgeIndexToFaceIndices.push_back({static_cast<int>(faceIndex)});
                endpointIndicesToEdgeIndex[firstEndpointIndex][secondEndpointIndex] = num_of_edges;
                endpointIndicesToEdgeIndex[secondEndpointIndex][firstEndpointIndex] = num_of_edges;
            }
            else
            {
                int existingEdgeIndex = endpointIndicesToEdgeIndex[firstEndpointIndex][secondEndpointIndex];
                edgeIndexToFaceIndices[existingEdgeIndex].push_back(faceIndex); // Connect it to its other face
            }
        }
    }

    /*create interpolation points if they exist, one per edge with material change*/

    vector<pair<Eigen::Vector2f, bool>> edgeIndexToMidPoints(edgeIndexToFaceIndices.size(), {{0, 0}, false}); // Second value is whether it has been set
    for (int edgeIndex = 0; edgeIndex < edgeIndexToEndpointIndices.size(); edgeIndex++)
    {
        const pair<int, int> &endpointIndices = edgeIndexToEndpointIndices[edgeIndex];

        // If the two materials at the adjacent points are different
        if (primaryMaterialIndexByPointIndex[endpointIndices.first] != primaryMaterialIndexByPointIndex[endpointIndices.second])
        {
            const int &endpt0Index = endpointIndices.first;
            const int &endpt1Index = endpointIndices.second;

            const int &endpt0PrimaryValueIndex = primaryMaterialIndexByPointIndex[endpt0Index];
            const int &endpt1PrimaryValueIndex = primaryMaterialIndexByPointIndex[endpt1Index];

            const vector<float> &primaryValues0 = pointIndexToMaterialValues[endpt0Index];
            const vector<float> &primaryValues1 = pointIndexToMaterialValues[endpt1Index];

            edgeIndexToMidPoints[edgeIndex] = {interpEdge2Mat<2>(
                                                   pointIndexToPoint.col(endpointIndices.first),
                                                   pointIndexToPoint.col(endpointIndices.second),
                                                   {primaryValues0[endpt0PrimaryValueIndex], primaryValues0[endpt1PrimaryValueIndex]},
                                                   {primaryValues1[endpt0PrimaryValueIndex], primaryValues1[endpt1PrimaryValueIndex]}),
                                               true};
        }
    }

    vector<int> edgeIndexToFacePointIndex(edgeIndexToEndpointIndices.size(), -1);
    /*create vertices in faces, one per triangle with material change*/
    // TODO: Increase speed
    // 1) Make slice tests
    // 2) Make tests on real data
    // 3) Time Benchmark for speed
    vector<Eigen::Vector2f> facePointByIndex;
    facePointByIndex.reserve(numberOfTriangles * 2);
    vector triangleIndexToFacePointIndex(numberOfTriangles, -1);
    for (size_t triangleIndex = 0; triangleIndex < triangleIndexToCornerIndices.size(); triangleIndex++)
    {
        vector<int> cornerIndices = triangleIndexToCornerIndices[triangleIndex];
        // If there are material changes in the triangle
        if (primaryMaterialIndexByPointIndex[cornerIndices[0]] != primaryMaterialIndexByPointIndex[cornerIndices[1]] || primaryMaterialIndexByPointIndex[cornerIndices[1]] != primaryMaterialIndexByPointIndex[cornerIndices[2]])
        {
            // Generate the center point of the existing midpoint
            vector<Eigen::Vector2f> triangleMidpoints;
            triangleMidpoints.reserve(triangle_edges.size());
            for (auto edge : triangle_edges)
            {
                int endpointIndices = endpointIndicesToEdgeIndex[cornerIndices[edge.first]][cornerIndices[edge.second]];
                if (edgeIndexToMidPoints[endpointIndices].second)
                    triangleMidpoints.push_back(edgeIndexToMidPoints[endpointIndices].first);
            }

            Eigen::Vector2f centerPoint = getMassPoint<2>(triangleMidpoints);
            facePointByIndex.push_back(centerPoint);

            // Create a face vertex at that generated center point
            triangleIndexToFacePointIndex[triangleIndex] = static_cast<int>(facePointByIndex.size()) - 1;
        }
    }

    /*create center segments*/
    vector<pair<int, int>> centerSegmentIndexToEndpointIndices;
    centerSegmentIndexToEndpointIndices.reserve(facePointByIndex.size() * 2);
    vector<pair<int, int>> centerSegmentToEndpointPrimaryMaterialIndices;
    centerSegmentToEndpointPrimaryMaterialIndices.reserve(facePointByIndex.size() * 2);

    for (int edge_index = 0; edge_index < edgeIndexToEndpointIndices.size(); edge_index++)
    {
        // If there is a change in material between the edges of the segments
        if (primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[edge_index].first] != primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[edge_index].second])
        {
            // The index is going to be the size of the list of segments before the push

            const pair<int, int> &endpointIndices = edgeIndexToEndpointIndices[edge_index];
            centerSegmentToEndpointPrimaryMaterialIndices.emplace_back(
                primaryMaterialIndexByPointIndex[endpointIndices.first],
                primaryMaterialIndexByPointIndex[endpointIndices.second]);

            pair<int, int> newSegmentEndpointIndices;
            if (edgeIndexToFaceIndices[edge_index].size() == 1)
            {
                // Boundary edge, connect edge point and triangle point
                facePointByIndex.push_back(edgeIndexToMidPoints[edge_index].first);
                const int facePointIndex = facePointByIndex.size() - 1;
                edgeIndexToFacePointIndex[edge_index] = facePointIndex;
                newSegmentEndpointIndices = {triangleIndexToFacePointIndex[edgeIndexToFaceIndices[edge_index][0]], facePointIndex};
            }
            else
            {
                newSegmentEndpointIndices = {triangleIndexToFacePointIndex[edgeIndexToFaceIndices[edge_index][0]], triangleIndexToFacePointIndex[edgeIndexToFaceIndices[edge_index][1]]};
            }

            // Ensure consistent orientation among endpoints
            if (orientation(pointIndexToPoint.col(edgeIndexToEndpointIndices[edge_index].first) - pointIndexToPoint.col(edgeIndexToEndpointIndices[edge_index].second), facePointByIndex[newSegmentEndpointIndices.first] - pointIndexToPoint.col(edgeIndexToEndpointIndices[edge_index].second)) < 0)
            {
                newSegmentEndpointIndices = {newSegmentEndpointIndices.second, newSegmentEndpointIndices.first};
            }

            centerSegmentIndexToEndpointIndices.push_back(std::move(newSegmentEndpointIndices));
        }
    }

    /*Create Fill Triangles -> Solid fill triangles for displaying areas of material type, not the contour*/
    vector<Eigen::Vector2f> resultingPointsByIndex;
    resultingPointsByIndex.reserve(facePointByIndex.size() + pointIndexToPoint.cols()); // We know how big this will be
    for (const auto &centerPoint : facePointByIndex)
    {
        resultingPointsByIndex.push_back(centerPoint);
    }
    // Join both sets of points into all the existing points
    for (const auto &pt : pointIndexToPoint.colwise())
    {
        resultingPointsByIndex.emplace_back(pt);
    }
    vector<vector<int>> resultingTriangleIndexToResultingCornerIndices;
    resultingTriangleIndexToResultingCornerIndices.reserve(numberOfTriangles * 6);
    vector<int> fillMats;
    fillMats.reserve(numberOfTriangles * 6);

    /*first type of triangles : dual to mesh edges with a material change*/
    for (int currentEdgeIndex = 0; currentEdgeIndex < edgeIndexToEndpointIndices.size(); currentEdgeIndex++)
    {
        // If the materials on either end of an edge don't match, there will be a triangle
        if (primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[currentEdgeIndex].first] != primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[currentEdgeIndex].second])
        {
            pair<int, int> segment_endpoints;
            if (edgeIndexToFaceIndices[currentEdgeIndex].size() == 1) // if boundary edge
            {
                if (edgeIndexToFacePointIndex[currentEdgeIndex] == -1)
                    throw "Logic error, midpoint index not found in face points";
                /*boundary edge : connect edge point and center point*/
                segment_endpoints = {triangleIndexToFacePointIndex[edgeIndexToFaceIndices[currentEdgeIndex][0]], edgeIndexToFacePointIndex[currentEdgeIndex]};
            }
            else
            {
                /*interior edge : connect two center points*/
                // We know both of these will exist b/c there is always a center point if the triangle has any material changes
                segment_endpoints = {
                    triangleIndexToFacePointIndex[edgeIndexToFaceIndices[currentEdgeIndex][0]],
                    triangleIndexToFacePointIndex[edgeIndexToFaceIndices[currentEdgeIndex][1]]};
            }

            // Insert two triangles&materials using the two endpoints of the currentEdge as the third corner

            // Insert the "top" triangle
            resultingTriangleIndexToResultingCornerIndices.push_back(vector({
                segment_endpoints.first,
                segment_endpoints.second,
                static_cast<int>(facePointByIndex.size()) + edgeIndexToEndpointIndices[currentEdgeIndex].first // B/c the edge indices start after the face point indices
            }));
            fillMats.push_back(primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[currentEdgeIndex].first]);

            // Insert the "bottom" triangle
            resultingTriangleIndexToResultingCornerIndices.push_back(vector({segment_endpoints.first,
                                                                                segment_endpoints.second,
                                                                                static_cast<int>(facePointByIndex.size()) + edgeIndexToEndpointIndices[currentEdgeIndex].second}));
            fillMats.push_back(primaryMaterialIndexByPointIndex[edgeIndexToEndpointIndices[currentEdgeIndex].second]);
        }
    }

    /* second type of triangles: original mesh triangle, if there is no material change,
     or a third of the triangle, if there is some edge with no material change */
    const int startingSize = facePointByIndex.size();
    for (int currentTriangleIndex = 0; currentTriangleIndex < triangleIndexToCornerIndices.size(); currentTriangleIndex++)
    {
        // If There is no change at all
        if (primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][1]] == primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][2]] && primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][2]] == primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][0]])
        {
            const vector<int> &corner_indices = triangleIndexToCornerIndices[currentTriangleIndex];
            vector<int> resulting_corner_indices;
            {
                resulting_corner_indices.reserve(corner_indices.size());
                for (int corner_index : corner_indices)
                {
                    resulting_corner_indices.push_back(corner_index + startingSize);
                }
            }
            resultingTriangleIndexToResultingCornerIndices.push_back(resulting_corner_indices);
            fillMats.push_back(primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][0]]);
        }
        else
        {
            for (int j = 0; j < 3; j++)
            {
                // If one of the edges has no change
                const pair<int, int> &triangle_side = triangle_edges[j];
                const vector<int> &cornerIndices = triangleIndexToCornerIndices[currentTriangleIndex];
                if (primaryMaterialIndexByPointIndex[cornerIndices[triangle_side.first]] ==
                    primaryMaterialIndexByPointIndex[cornerIndices[triangle_side.second]])
                {
                    resultingTriangleIndexToResultingCornerIndices.push_back({cornerIndices[triangle_side.first] + static_cast<int>(facePointByIndex.size()),
                                                                              cornerIndices[triangle_side.second] + static_cast<int>(facePointByIndex.size()),
                                                                              triangleIndexToFacePointIndex[currentTriangleIndex]});
                    fillMats.push_back(primaryMaterialIndexByPointIndex[triangleIndexToCornerIndices[currentTriangleIndex][triangle_side.second]]);
                }
            }
        }
    }

    return {std::move(facePointByIndex),
            std::move(centerSegmentIndexToEndpointIndices),
            std::move(centerSegmentToEndpointPrimaryMaterialIndices),
            std::move(resultingPointsByIndex),
            std::move(resultingTriangleIndexToResultingCornerIndices),
            std::move(fillMats)};
}

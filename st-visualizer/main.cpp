
#include "Contour3D.h"
#include "ImportFunctions.h"
#include "JSONParser.h"
#include "Stats.h"
#include "UtilityFunctions.h"
#include "Timing.h"
#include "PHExport.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace nlohmann;

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::stringstream;
using std::vector;

unsigned long contour_2d;
unsigned long contour_3d;
unsigned long stats;
unsigned long export_io;

bool ph_toggle;
string ph_points_path;
string ph_tets_path;
int wid_buffer;
int num_ransac;

// Mode 0: ./st-visualizer 0 <config.json file path>
// Mode 1: ./st-visualizer 1 <config.json file content>
int main(int argc, char *argv[])
{
    json config;
    if (strcmp(argv[1], "0") == 0)
    {
        ifstream file(argv[2]);
        if (file.is_open())
        {
            string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
            file.close();
            config = json::parse(content);
        }
        else
        {
            cerr << "Error opening the file" << endl;
            return -1;
        }
    }
    else if (strcmp(argv[1], "1") == 0)
    {
        config = json::parse(argv[2]);
    }
    else
    {
        std::cout << "Usage: " << argv[0] << " <mode> <argument>" << std::endl;
        return 1;
    }

    string alignmentFile = config.at("alignmentFile").get<string>();
    string target = config.at("target").get<string>();
    float shrink = config.at("shrink").get<float>();
    vector<string> sliceNames;
    for (const auto &name : config.at("sliceNames"))
    {
        if (name.is_string())
        {
            sliceNames.push_back(name.get<std::string>());
        }
    }
    vector<unsigned> featureCols;
    for (const auto &name : config.at("featureCols"))
    {
        featureCols.push_back(name.get<unsigned>());
    }

    ph_toggle = config.at("PHExport").get<bool>();
    ph_points_path = config.at("PHPoints").get<string>();
    ph_tets_path = config.at("PHTets").get<string>();
    wid_buffer = config.at("GrowWidth").get<int>();
    num_ransac = config.at("NumRansac").get<int>();

    const vector<pair<vector<coord>, vector<coord>>> alignmentValues = importAlignments(alignmentFile);

    const tsv_return_type results = loadTsv(
        config.at("fileName").get<std::string>(),
        sliceNames,
        config.at("sliceIndex").get<int>(),
        config.at("tissueIndex").get<int>(),
        std::pair<unsigned, unsigned>(config.at("rowIndex").get<int>(), config.at("colIndex").get<int>()),
        config.at("clusterIndex").get<int>(),
        featureCols,
        config.at("zDistance").get<int>(),
        alignmentValues);

    std::chrono::steady_clock::time_point start_contour_2d = std::chrono::high_resolution_clock::now();
    auto [ctrs2dVals, tris2dVals] = getSectionContoursAll(results.slices, results.values, shrink);
    auto [ctrs2dclusters, tris2dclusters] = getSectionContoursAll(results.slices, results.clusters, shrink);
    std::chrono::steady_clock::time_point end_contour_2d = std::chrono::high_resolution_clock::now();
    contour_2d = duration_cast<std::chrono::microseconds>(end_contour_2d - start_contour_2d).count();

    std::chrono::steady_clock::time_point start_contour_3d = std::chrono::high_resolution_clock::now();
    auto allpts = concatMatrixes(results.slices);
    auto ctrs3dVals = getVolumeContours(allpts, flatten<std::vector<float>>(results.values), shrink, true);
    auto ctrs3dClusters = getVolumeContours(allpts, flatten<std::vector<float>>(results.clusters), shrink, false);
    auto ptClusIndex = mapVector(results.clusters, std::function(
                                                       [](const std::vector<std::vector<float>> &layer, size_t)
                                                       {
                                                           return mapVector(layer, std::function(getMaxPos));
                                                       }));
    auto ptValIndex = mapVector(results.values, std::function([](const std::vector<std::vector<float>> &layer)
                                                              { return mapVector(layer, std::function(getMaxPos)); }));
    auto slices = mapVector(results.slices, std::function([](const Eigen::Matrix3Xf &layer)
                                                          {
                                                              std::vector<Eigen::Vector3f> temp;
                                                              temp.reserve(layer.cols());
                                                              for (const auto &pt: layer.colwise())
                                                              {
                                                                  temp.emplace_back(pt);
                                                              }
                                                              return temp; }));

    auto convertCtrs = [](
                           std::vector<std::vector<std::pair<
                               std::vector<Eigen::Matrix<float, 3, 1, 0>>, std::vector<std::pair<int, int>>>>> &ctrs2dVals)
    {
        json ctrs2dValsJson = json::array();
        for (auto &ctrSlice : ctrs2dVals)
        {
            json ctrJson = json::array();
            for (auto &ctr : ctrSlice)
            {
                json temp = json::array();
                temp.push_back(ctr.first);
                for (auto &i : ctr.second)
                {
                    // i.first++;
                    // i.second++;
                }
                temp.push_back(ctr.second);
                ctrJson.push_back(temp);
            }

            ctrs2dValsJson.push_back(ctrJson);
        }

        return ctrs2dValsJson;
    };

    auto convertTris = [](std::vector<std::tuple<
                              std::vector<Eigen::Matrix<float, 3, 1, 0>>,
                              std::vector<std::vector<int>>,
                              std::vector<int>>> &tris2dVals)
    {
        json tris2dValsJson = json::array();
        for (auto &tris : tris2dVals)
        {
            json a = json::array();
            
            json b = json::array();
            for (auto &elem : std::get<0>(tris))
            {
                b.push_back(std::vector(elem.data(), elem.data() + elem.rows()));
            }
            a.push_back(b);

            auto &triangles = std::get<1>(tris);
            a.push_back(triangles);

            auto &materials = std::get<2>(tris);
            a.push_back(materials);

            tris2dValsJson.push_back(a);
        }

        return tris2dValsJson;
    };

    auto convert3D = [](
                         std::vector<std::pair<std::vector<Eigen::Vector3f>, std::vector<std::vector<int>>>> &ctrs3d)
    {
        json ctrs3dJson = json::array();
        for (auto &ctr : ctrs3d)
        {
            json a = json::array();

            a.push_back(ctr.first);
            auto &segs = ctr.second;
            for (auto &a : segs)
            {
                for (auto &b : a)
                {
                    // b++;
                }
            }
            a.push_back(segs);

            ctrs3dJson.push_back(a);
        }

        return ctrs3dJson;
    };
    std::chrono::steady_clock::time_point end_contour_3d = std::chrono::high_resolution_clock::now();
    contour_3d = duration_cast<std::chrono::microseconds>(end_contour_3d - start_contour_3d).count();

    std::chrono::steady_clock::time_point start_stats = std::chrono::high_resolution_clock::now();
    vector<float> surface_area_features = computeSurfaceArea(ctrs3dVals);
    vector<float> surface_area_clusters = computeSurfaceArea(ctrs3dClusters);
    vector<float> volume_features = computeVolume(ctrs3dVals);
    vector<float> volume_clusters = computeVolume(ctrs3dClusters);
    auto [componentsVals, handlesVals] = connectedComponent(ctrs3dVals);
    auto [componentsClusters, handlesClusters] = connectedComponent(ctrs3dClusters);
    std::chrono::steady_clock::time_point end_stats = std::chrono::high_resolution_clock::now();
    stats = duration_cast<std::chrono::microseconds>(end_stats - start_stats).count();

    std::chrono::steady_clock::time_point start_export_io = std::chrono::high_resolution_clock::now();
    json ret = json::object();
    ret["nat"] = results.values[0][0].size(); // nMat,
    ret["shrink"] = shrink;                   // shrink,
    ret["clusters"] = results.clusters;       // clusters,
    ret["slices"] = slices;                   // slices,
    ret["sliceNames"] = sliceNames;           // sliceNames
    ret["values"] = results.values;
    ret["featureNames"] = results.names;                 // featureNames,
    ret["featureCols"] = featureCols;                    // featureCols,
    ret["ptValIndex"] = ptValIndex;                      // ptValIndex,
    ret["ctrs2Dvals"] = convertCtrs(ctrs2dVals);         // ctrs2Dvals,
    ret["tris2Dvals"] = convertTris(tris2dVals);         // tris2Dvals
    ret["ctrs3Dvals"] = convert3D(ctrs3dVals);           // ctrs3Dvals,
    ret["ptClusIndex"] = ptClusIndex;                    // ptClusIndex,
    ret["nClusters"] = results.clusters[0][0].size();    // nClusters,
    ret["ctrs2Dclusters"] = convertCtrs(ctrs2dclusters); // ctrs2Dclusters,
    ret["tris2Dclusters"] = convertTris(tris2dclusters); // tris2Dclusters,
    ret["ctrs3Dclusters"] = convert3D(ctrs3dClusters);   // ctrs3Dclusters,
                                                         //    ret["ctrsVolumeClusters"] = getVolumes(ctrs3dClusters);

    if (config.at("objExport").get<bool>())
    {
        log("Exporting obj files.");
        exportObj(config.at("featureObj").get<string>(), ctrs3dVals, results.names);
        exportObj(config.at("clusterObj").get<string>(), ctrs3dClusters, results.clusterNames);
    }

    ret["ctrsSurfaceAreaVals"] = surface_area_features;
    ret["ctrsSurfaceAreaClusters"] = surface_area_clusters;
    ret["ctrsVolumeVals"] = volume_features;
    ret["ctrsVolumeClusters"] = volume_clusters;
    ret["componentsVals"] = componentsVals;
    ret["componentsClusters"] = componentsClusters;
    ret["handlesVals"] = handlesVals;
    ret["handlesClusters"] = handlesClusters;

    log("Calculations complete.");

    if (config.at("resultExport").get<bool>())
    {
        log("Complete. Writing to file.");
        std::ofstream f(target);
        f << ret;
    }
    std::chrono::steady_clock::time_point end_export_io = std::chrono::high_resolution_clock::now();
    export_io = duration_cast<std::chrono::microseconds>(end_export_io - start_export_io).count();

    if (config.at("timingExport").get<bool>())
    {
        string timing_path = config.at("timing").get<string>();
        export_timing(results.num_points, timing_path);
    }

    log("Exiting.");
    return 0;
}

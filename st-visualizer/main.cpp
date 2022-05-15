#include <string>
#include <fstream>
#include <iostream>
#include "Contour2D.h"
#include "JSONParser.h"
#include "UtilityFunctions.h"
using namespace nlohmann;
//https://wias-berlin.de/software/tetgen/


int main(int argc, char** argv)
{
    using json = json;
    int i = 0;

    auto jsonToMatrix = [](const json& source)
    {
        //This should only be used for testing. Not logic safe.
        Eigen::Matrix2Xf a(2, source.size());

        for(int i = 0; i < source.size(); i++)
        {
            a(0, i) = source[i][0];
            a(1, i) = source[i][1];
        }
        return a;
    };

    auto jsonToVector = [](const json& source)
    {
        std::vector<std::vector<float>> ret;
        for(auto row : source)
        {
            std::vector<float> temp;
            for(auto item : row)
            {
                temp.push_back(item);
            }
            ret.push_back(temp);
        }
        return ret;
    };

    auto jsonToTris = [](const json& source)
    {
        std::vector<std::vector<int>> ret;
        for(auto row : source)
        {
            std::vector<int> temp;
            for(auto item : row)
            {
                temp.push_back(static_cast<int>(item) - 1); //Subtract 1 b/c mathematica indices start a 1
            }
            ret.push_back(temp);
        }
        return ret;
    };
    // read a JSON file
    std::ifstream
        file("C:\\Users\\Aiden McIlraith\\Documents\\GitHub\\st-visualizer\\UnitTest\\singleContourTest.json");
    auto a = file.is_open();
    json j2 = json::parse(file);
    json ret2 = json::array();
    for(auto j : j2)
    {
        auto pts = jsonToMatrix(j[0]);
        auto vals = jsonToVector(j[1]);
        auto delaunay = triangulateMatrix(pts);
        auto tris = table(static_cast<size_t>(delaunay.numberoftriangles), std::function([delaunay](size_t i)
        {
            return getTriangleCornerIndices(delaunay, i);
        }));

        auto first = contourTriMultiDC(pts, tris, vals);

        json ret = json::array();

        std::vector<Eigen::Vector2f> temp;
        for(auto& vert : first.fillVerts)
        {
            temp.push_back(vert);
        }
        ret.push_back(temp);
        ret.push_back(first.fillTris);

        ret.push_back(first.fillMats);
        ret2.push_back(ret);
    }
    std::ofstream f(
        "C:\\Users\\Aiden McIlraith\\Documents\\GitHub\\st-visualizer\\UnitTest\\singleContourResultsFaces.json");
    f << ret2;


    // const auto alignmentValues = importAlignments("C:/Users/Aiden McIlraith/Documents/GitHub/st-visualizer/NMK_F_transformation_pt_coord.csv");
    // const auto results = loadTsv(
    // 	"C:/Users/Aiden McIlraith/Documents/GitHub/st-visualizer/NMK_20201201_cell_type_coord_allspots.tsv",
    // 	std::vector<std::string>({ "NMK_F_U1","NMK_F_U2","NMK_F_U3","NMK_F_U4" }),
    // 	1,
    // 	2,
    // 	std::pair<unsigned, unsigned>(3, 
    // 	5,
    // 	std::vector<unsigned>({ 6,7,8,9, 10, 11, 12, 13, 14, 15 }),
    // 	60,
    // 	alignmentValues
    // );

    // getSectionContoursAll(results.slices, results.values, 1);
}

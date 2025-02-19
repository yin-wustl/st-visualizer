//
// Created by Yin Li on 5/18/24.
//

#ifndef ST_VISUALIZER_PHEXPORT_H
#define ST_VISUALIZER_PHEXPORT_H

#include <Eigen/Eigen>
#include <vector>
#include <iostream>
#include <fstream>

using std::vector;
using std::string;
using std::stringstream;

extern bool ph_toggle;
extern string ph_points_path;
extern string ph_tets_path;

inline void export_ph(const vector<Eigen::Vector3f> &points, const vector<vector<float>> &materials, const vector<vector<int>> &tets)
{
    if (!ph_toggle)
    {
        return;
    }

    std::ofstream file_points(ph_points_path, std::ios_base::out);
    std::ofstream file_tets(ph_tets_path, std::ios_base::out);
    if (file_points.is_open() && file_tets.is_open())
    {
        for (int i = 0; i < points.size(); i++)
        {
            const Eigen::Vector3f &coord = points.at(i);
            stringstream ss;
            ss.setf(std::ios::fixed);
            ss.precision(16);
            ss << coord(0) << "," << coord(1) << "," << coord(2);

            for (float material : materials.at(i))
            {
                ss << "," << material;
            }
            ss << std::endl;
            file_points << ss.rdbuf();
        }

        for (int i = 0; i < tets.size(); i++)
        {
            const vector<int> &tet = tets.at(i);
            file_tets << tet.at(0) << "," << tet.at(1) << "," << tet.at(2) << "," << tet.at(3) << std::endl;
        }

        file_points.close();
        file_tets.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file");
    }
}

#endif //ST_VISUALIZER_PHEXPORT_H

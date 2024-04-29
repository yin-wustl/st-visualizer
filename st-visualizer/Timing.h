//
// Created by Yin Li on 4/22/24.
//

#include "Contour2D.h"

#include <iostream>
#include <fstream>
#include <vector>

#ifndef ST_VISUALIZER_TIMING_H
#define ST_VISUALIZER_TIMING_H

using std::vector;

extern unsigned long import_io;
extern unsigned long preprocessing;
extern unsigned long cover_and_grow;
extern unsigned long contour_2d;
extern unsigned long contour_3d;
extern unsigned long stats;
extern unsigned long export_io;
extern vector<unsigned long> contour_tetgen;
extern vector<unsigned long> contour_triangle;

inline void export_timing(int num_points) {
    std::ofstream file("timing.csv", std::ios_base::app);
    if (file.is_open())
    {
        file << num_points << ",";
        file << import_io << ",";
        file << preprocessing << ",";
        file << cover_and_grow << ",";
        file << contour_2d << ",";
        file << contour_3d << ",";
        file << stats << ",";
        file << export_io << ",";
        file << contour_tetgen[0] << "," << contour_tetgen[1] << ",";
        file << contour_triangle[0] << "," << contour_triangle[1] << "\n";
//        file << std::endl;
        file.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file");
    }
}

#endif //ST_VISUALIZER_TIMING_H

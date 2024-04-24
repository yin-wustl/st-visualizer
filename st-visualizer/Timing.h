//
// Created by Yin Li on 4/22/24.
//

#include <iostream>
#include <fstream>

#ifndef ST_VISUALIZER_TIMING_H
#define ST_VISUALIZER_TIMING_H

extern unsigned long import_io;
extern unsigned long preprocessing;
extern unsigned long cover_and_grow;
extern unsigned long contour_2d;
extern unsigned long contour_3d;
extern unsigned long stats;
extern unsigned long export_io;
extern unsigned long contour_tetgen;

inline void export_timing(int num_points) {
//    std::cout << num_points << std::endl;
//    std::cout << import_io << std::endl;
//    std::cout << preprocessing << std::endl;
//    std::cout << cover_and_grow << std::endl;
//    std::cout << contour_2d << std::endl;
//    std::cout << contour_3d << std::endl;
//    std::cout << stats << std::endl;
//    std::cout << export_io << std::endl;

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
        file << contour_tetgen << "\n";
        file.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file");
    }
}

#endif //ST_VISUALIZER_TIMING_H

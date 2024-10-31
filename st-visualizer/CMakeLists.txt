cmake_minimum_required(VERSION 3.25)
project(st-visualizer)

set(CMAKE_CXX_STANDARD 23)

include_directories(.)
include_directories(eigen-3.4.0)
include_directories(tetgen-1.6.0)
include_directories(triangle-1.6)
include_directories(phat-1.7/include)
#include_directories(phat-1.6/include)

add_executable(st-visualizer
        tetgen-1.6.0/predicates.cxx
        tetgen-1.6.0/tetgen.cxx
        tetgen-1.6.0/tetgen.h
        triangle-1.6/triangle.c
        triangle-1.6/triangle.h
        triangle-1.6/tricall.c
        Contour2D.cpp
        Contour2D.h
        Contour3D.cpp
        Contour3D.h
        GrowAndCover.cpp
        GrowAndCover.h
        ImportFunctions.cpp
        ImportFunctions.h
        JSONParser.h
        main.cpp
        Stats.h
        UtilityFunctions.cpp
        UtilityFunctions.h
        Timing.h
        PHExport.h
        PHCompute.h)

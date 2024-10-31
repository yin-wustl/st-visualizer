# FIXME: obsolete, do not use makefile to compile the project

CC = clang
CXX = clang++
BIN = ./bin/
CCFLAGS = -O3 -g -std=c99 -Wfatal-errors -w -lm -arch x86_64 -arch arm64
CXXFLAGS = -O3 -g -std=c++2a -Wfatal-errors -w -lm -arch x86_64 -arch arm64
INCLUDES = -Ieigen-3.4.0 -Itetgen-1.6.0 -Itriangle-1.6
DEPENDENCIES = $(BIN)tetgen.o $(BIN)predicates.o $(BIN)tricall.o $(BIN)triangle.o
OBJECTS = $(wildcard $(BIN)*.o)

ifeq ($(OS), Windows_NT)
	EXECUTABLE=st-visualizer.exe
	CCFLAGS = -Og -g -std=c99 -Wfatal-errors -w 
	CXXFLAGS = -Og -g -std=c++2a -Wfatal-errors -w 
else
	EXECUTABLE=st-visualizer
endif

release: $(wildcard ./*.cpp) dependencies
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(BIN)$(EXECUTABLE) $(wildcard ./*.cpp) $(DEPENDENCIES)

test: main.cpp dependencies objects
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(BIN)$(EXECUTABLE) $(OBJECTS)

dependencies: $(BIN)tetgen.o $(BIN)tricall.o

objects: $(BIN)import.o $(BIN)utility.o $(BIN)contour-2d.o $(BIN)contour-3d.o $(BIN)grow-and-cover.o $(BIN)main.o

$(BIN)import.o: ImportFunctions.cpp ImportFunctions.h $(BIN)utility.o $(BIN)grow-and-cover.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)import.o ImportFunctions.cpp $(BIN)utility.o $(BIN)grow-and-cover.o

$(BIN)utility.o: UtilityFunctions.cpp UtilityFunctions.h JSONParser.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)utility.o UtilityFunctions.cpp

$(BIN)contour-2d.o: Contour2D.cpp Contour2D.h $(BIN)utility.o 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)contour-2d.o Contour2D.cpp $(BIN)utility.o 

$(BIN)contour-3d.o: Contour3D.cpp Contour3D.h 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)contour-3d.o Contour3D.cpp $(BIN)utility.o $(BIN)contour-2d.o

$(BIN)grow-and-cover.o: GrowAndCover.cpp GrowAndCover.h $(BIN)utility.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)grow-and-cover.o GrowAndCover.cpp $(BIN)utility.o

$(BIN)stats.o: Stats.cpp Stats.h $(BIN)utility.o $(BIN)import.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)stats.o Stats.cpp $(BIN)utility.o $(BIN)import.o

$(BIN)main.o: $(BIN)contour-3d.o $(BIN)import.o $(BIN)utility.o 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BIN)main.o main.cpp $(BIN)contour-3d.o $(BIN)import.o $(BIN)utility.o 


$(BIN)tetgen.o: tetgen-1.6.0/tetgen.cxx tetgen-1.6.0/tetgen.h $(BIN)predicates.o
	$(CXX) $(CXXFLAGS) -c -o $(BIN)tetgen.o tetgen-1.6.0/tetgen.cxx $(BIN)predicates.o

$(BIN)predicates.o: tetgen-1.6.0/predicates.cxx
	$(CXX) $(CXXFLAGS) -c -o $(BIN)predicates.o tetgen-1.6.0/predicates.cxx

$(BIN)tricall.o: triangle-1.6/tricall.c $(BIN)triangle.o
	$(CC) $(CCFLAGS) -c -o $(BIN)tricall.o triangle-1.6/tricall.c $(BIN)triangle.o

$(BIN)triangle.o: triangle-1.6/triangle.c triangle-1.6/triangle.h 
	$(CC) $(CCFLAGS) -c -o $(BIN)triangle.o triangle-1.6/triangle.c

.PHONY: clean
clean:
	rm -rf $(BIN)*

.PHONY: clean-main
clean-main:
	rm -rf $(BIN)st-visualizer
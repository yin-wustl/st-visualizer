#!/bin/zsh

echo $(uname -a) >> time.txt

# {time ./bin/st-visualizer 0; } 2>> time.txt

{time ./bin/st-visualizer 1 '{"fileName":"../picture/AlignmentImages/NMKimage/NMK_20201201_cell_type_coord_allspots.tsv","alignmentFile": "../data/NMK_F_transformation_pt_coord.csv","target": "bin/nmk-test-output.json","shrink":0,"sliceNames":["NMK_F_U1","NMK_F_U2","NMK_F_U3","NMK_F_U4"],"featureCols":[6,7,8,9,10,11,12,13,14,15],"sliceIndex":1,"tissueIndex":2,"rowIndex":3,"colIndex":4,"clusterIndex":5,"zDistance":100}'; } 2>> time.txt

{time ./bin/st-visualizer 1 '{"fileName":"../picture/AlignmentImages/CRChi_res/CRC_112C1_cell_type_coord_allspots.tsv","alignmentFile": "../data/NMK_F_transformation_pt_coord.csv","target": "bin/nmk-test-output.json","shrink":0,"sliceNames":["CRC_HT112C1_1","CRC_HT112C1_2"],"featureCols":[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],"sliceIndex":1,"tissueIndex":2,"rowIndex":3,"colIndex":4,"clusterIndex":7,"zDistance":100}'; } 2>> time.txt

echo -en '\n' >> time.txt
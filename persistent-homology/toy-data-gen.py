import numpy as np
import pandas as pd
import math
import json
import os

# grid properties
name = "toy-data-1"
point_distance = 2
v_1 = np.array([point_distance, 0])
v_2 = np.array([point_distance * np.cos(np.pi / 3), point_distance * np.sin(np.pi / 3)])

slices = [
    [
        ((0, 0), [0.3, 0.7]), ((1, 0), [0.3, 0.7]), ((2, 0), [0.3, 0.7]), ((3, 0), [0.3, 0.7]),
        ((0, 1), [0.2, 0.8]), ((1, 1), [0.2, 0.8]), ((2, 1), [0.2, 0.8]), ((0, 2), [0.2, 0.8]), ((1, 2), [0.2, 0.8]), ((0, 3), [0.2, 0.8]),
        ((1, -1), [0.2, 0.8]), ((2, -1), [0.2, 0.8]), ((3, -1), [0.2, 0.8]), ((2, -2), [0.2, 0.8]), ((3, -2), [0.2, 0.8]), ((3, -3), [0.2, 0.8])
    ], 
    [
        ((0, 0), [0.8, 0.2]), ((1, 0), [0.8, 0.2]), ((2, 0), [0.8, 0.2]), ((3, 0), [0.8, 0.2]),
        ((0, 1), [0.2, 0.8]), ((1, 1), [0.2, 0.8]), ((2, 1), [0.2, 0.8]), ((0, 2), [0.2, 0.8]), ((1, 2), [0.2, 0.8]), ((0, 3), [0.2, 0.8]),
        ((1, -1), [0.2, 0.8]), ((2, -1), [0.2, 0.8]), ((3, -1), [0.2, 0.8]), ((2, -2), [0.2, 0.8]), ((3, -2), [0.2, 0.8]), ((3, -3), [0.2, 0.8])
    ],
    [
        ((0, 0), [0.3, 0.7]), ((1, 0), [0.3, 0.7]), ((2, 0), [0.3, 0.7]), ((3, 0), [0.3, 0.7]),
        ((0, 1), [0.2, 0.8]), ((1, 1), [0.2, 0.8]), ((2, 1), [0.2, 0.8]), ((0, 2), [0.2, 0.8]), ((1, 2), [0.2, 0.8]), ((0, 3), [0.2, 0.8]),
        ((1, -1), [0.2, 0.8]), ((2, -1), [0.2, 0.8]), ((3, -1), [0.2, 0.8]), ((2, -2), [0.2, 0.8]), ((3, -2), [0.2, 0.8]), ((3, -3), [0.2, 0.8])
    ], 
]

num_slices = len(slices)
num_features = len(slices[0][0][1])

# convert grid coordinates to cartesian coordinates
def grid_to_cartesian(a, b):
    v = a * v_1 + b * v_2
    return v[0], v[1]


headers = ["slice_id", "is_tissue", "cluster", "row", "column"] + [f"feature-{i}" for i in range(num_features)]
output = pd.DataFrame(columns=headers)

for i, slice in enumerate(slices):
    data = []
    slice_id = "slice_{}".format(i)
    for j, point in enumerate(slice):
        x, y = grid_to_cartesian(*point[0])
        is_tissue = 1
        features = point[1]
        cluster = 0
        data.append([slice_id, is_tissue, cluster, x, y] + features)
    table = pd.DataFrame(data, columns=headers)
    output = pd.concat([output, table])

output.to_csv(f"persistent-homology/assets/{name}/data.tsv", index=False, sep="\t")

# generate alignment file
with open(f"persistent-homology/assets/{name}/alignment.csv", "w") as file:
    file.write(" , , , , , ,\n")  # Writing a single line
    for i in range(num_slices - 1):
        file.writelines(
            [
                ", , ,0,100,-100,-100,100\n",
                ", , ,0,100,-100,100,-100\n",
                ", , ,0,100,-100,-100,100\n",
                ", , ,0,100,-100,100,-100\n",
            ]
        )
file.close()

# generate metadata file
current_path = os.getcwd()
output = {
    "fileName": f"{current_path}/persistent-homology/assets/{name}/data.tsv",
    "alignmentFile": f"{current_path}/persistent-homology/assets/{name}/alignment.csv",
    "target": f"{current_path}/persistent-homology/assets/{name}/output.json",
    "shrink": 0,
    "GrowWidth": 2,
    "NumRansac": 2,
    "sliceNames": [f"slice_{i}" for i in range(num_slices)],
    "featureCols": [5, 6],
    "sliceIndex": 0,
    "tissueIndex": 1,
    "rowIndex": 3,
    "colIndex": 4,
    "clusterIndex": 2,
    "zDistance": point_distance,
    "resultExport": False,
    "objExport": False,
    "featureObj": f"{current_path}/persistent-homology/assets/{name}/features/",
    "clusterObj": f"{current_path}/persistent-homology/assets/{name}/clusters/",
    "timingExport": False,
    "timing": f"{current_path}/persistent-homology/assets/{name}/timing.csv",
    "PHExport": True,
    "PHPoints": f"{current_path}/persistent-homology/assets/{name}/points.csv",
    "PHTets": f"{current_path}/persistent-homology/assets/{name}/tets.csv",
}

with open(f"persistent-homology/assets/{name}/config.json", "w") as file:
    file.write(json.dumps(output))
file.close()

# run st-visualizer
os.system(
    f"{current_path}/st-visualizer/cmake-build-release/st-visualizer 0 {current_path}/persistent-homology/assets/{name}/config.json"
)

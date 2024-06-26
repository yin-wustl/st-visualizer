import numpy as np
import pandas as pd
import math
import json
import os

# torus properties
R = 80
r = 20

# grid properties
point_distance = 4
slice_diameter = math.ceil(1.5 * (R + r) / point_distance)
num_slices_one_side = 4
num_slices = 2 * num_slices_one_side - 1
v_1 = np.array([point_distance, 0])
v_2 = np.array([point_distance * np.cos(np.pi / 3), point_distance * np.sin(np.pi / 3)])

# distance of each slice from the center slice
slice_distance = [i * (r / (num_slices_one_side + 1)) for i in range(-num_slices_one_side, num_slices_one_side + 1)]

# width of intersection of torus on each slice
widths = [2 * np.sqrt(r**2 - d**2) for d in slice_distance]

headers = ["slice_id", "is_tissue", "cluster", "row", "column", "feature", "feature_void"]
output = pd.DataFrame(columns=headers)


# convert grid coordinates to cartesian coordinates
def grid_to_cartesian(a, b):
    v = a * v_1 + b * v_2
    return v[0], v[1]


# TODO: verify if this is correct
# convert cartesian coordinates to grid coordinates
def cartesian_to_grid(x, y):
    return np.dot([x, y], np.linalg.inv(v_1)), np.dot([x, y], np.linalg.inv(v_2))


# generate data
for i in range(num_slices):
    data = []
    slice_id = "slice_{}".format(i)
    for a in range(-slice_diameter, slice_diameter):
        for b in range(-slice_diameter, slice_diameter):
            x, y = grid_to_cartesian(a, b)
            distance = math.sqrt(x**2 + y**2)
            is_tissue = 1
            if R - widths[i] <= distance <= R + widths[i]:
                features = [math.fabs(r - math.sqrt((distance - R) ** 2 + slice_distance[i] ** 2)) / r, 0]
            else:
                features = [0, 1]
            cluster = 0
            data.append([slice_id, is_tissue, cluster, x, y] + features)
    table = pd.DataFrame(data, columns=headers)
    output = pd.concat([output, table])

min_feature = output["feature"].min()
max_feature = output["feature"].max()
print(f"Min feature: {min_feature}, Max feature: {max_feature}")

output.to_csv(f"persistent-homology/assets/data_torus.tsv", index=False, sep="\t")

# generate alignment file
with open("persistent-homology/assets/alignment.csv", "w") as file:
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
    "fileName": f"{current_path}/persistent-homology/assets/data_torus.tsv",
    "alignmentFile": f"{current_path}/persistent-homology/assets/alignment.csv",
    "target": f"{current_path}/persistent-homology/assets/output_torus.json",
    "shrink": 0,
    "sliceNames": [f"slice_{i}" for i in range(num_slices)],
    "featureCols": [5, 6],
    "sliceIndex": 0,
    "tissueIndex": 1,
    "rowIndex": 3,
    "colIndex": 4,
    "clusterIndex": 2,
    "zDistance": r / (num_slices_one_side + 1),
    "resultExport": False,
    "objExport": False,
    "featureObj": "{current_path}/persistent-homology/assets/features/",
    "clusterObj": "{current_path}/persistent-homology/assets/clusters/",
}

with open("persistent-homology/assets/config.json", "w") as file:
    file.write(json.dumps(output))
file.close()

# run st-visualizer
# /Users/yin/Code/masters-project/st-visualizer/st-visualizer/cmake-build-release/st-visualizer 0 /Users/yin/Code/masters-project/st-visualizer/persistent-homology/assets/config.json
# os.system(f"{current_path}/st-visualizer/cmake-build-release/st-visualizer 0 {current_path}/persistent-homology/assets/config.json")

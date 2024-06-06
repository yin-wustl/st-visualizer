import os
import math
from config_gen import config_gen
from data_gen import data_gen

start = 2000
end = 5000
step = 5

repeat = 5
num_slices = 10

exec_path = "/Users/yin/Code/masters-project/st-visualizer/st-visualizer/cmake-build-release/st-visualizer"
data_path = "/Volumes/tmp/"  # end in "/"
alignment_path = data_path + "alignments.csv"
size_set = set()

for size in range(start, end + 1, step):
    slice_diameter = math.ceil(math.sqrt(size))
    num_points = num_slices * (2 * slice_diameter) ** 2
    if num_points in size_set:
        continue
    size_set.add(num_points)
    print(num_points)

    this_data_path = data_path + str(num_points) + "/"
    data_gen(
        num_slices=num_slices,
        slice_diameter=slice_diameter,
        point_distance=10,
        num_points=num_points,
        num_clusters=20,
        num_features=10,
        prob_tissue=0.8,
        rotation_degree_range=2,
        path=this_data_path,
    )
    config_gen(str(num_points), this_data_path, alignment_path)
    for _ in range(repeat):
        command = f"{exec_path} 0 {this_data_path}config_{num_points}.json"
        print(command)
        os.system(f"{exec_path} 0 {this_data_path}config_{num_points}.json")

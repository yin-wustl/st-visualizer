import numpy as np
import pandas as pd

# Settings
num_slices = 4
slice_diameter = 50
point_distance = 10
num_clusters = 20
num_features = 10
prob_tissue = 0.8
headers = [
    'slice_id',
    'is_tissue',
    'cluster',
    'row_grid',
    'column_grid',
] + ['feature_{}'.format(i) for i in range(num_features)]


# Generate data
index = 0
data = []
for i in range(num_slices):
    slice_id = 'slice_{}'.format(i)
    for x in range(-slice_diameter, slice_diameter):
        for y in range(-slice_diameter, slice_diameter):
            if np.random.rand() > prob_tissue:
                is_tissue = 0
            else:
                is_tissue = 1
            cluster = np.random.randint(num_clusters)
            row_grid = x
            column_grid = y
            features = [np.random.rand() for _ in range(num_features)]
            features = [i/np.sum(features) for i in features]
            # table.loc[index] = [slice_id, is_tissue, cluster, row_grid, column_grid] + features
            data.append([slice_id, is_tissue, cluster, row_grid, column_grid] + features)
            index += 1

table = pd.DataFrame(data, columns=headers)

# Convert grid coordinates to cartesian coordinates
origin = np.array([0, 0])
v_1 = np.array([point_distance, 0])
v_2 = np.array([point_distance * np.cos(np.pi/3), point_distance * np.sin(np.pi/3)])

# convert coords to numpy array
grid_coords = table[['column_grid', 'row_grid']].to_numpy()
base = np.array([v_1, v_2])
cartesian_coords = np.dot(grid_coords, base)

# append the cartesian coordinates to the table
table['column'] = cartesian_coords[:, 0]
table['row'] = cartesian_coords[:, 1]

table.to_csv('data.tsv', index=False, sep="\t")


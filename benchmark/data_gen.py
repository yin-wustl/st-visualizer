import numpy as np
import pandas as pd
import os


def data_gen(
    num_slices: int,
    slice_diameter: int,
    point_distance: int,
    num_points: int,
    num_clusters: int,
    num_features: int,
    prob_tissue: float,
    rotation_degree_range: float,
    path: str,
):
    headers = [
        "slice_id",
        "is_tissue",
        "cluster",
        "row_grid",
        "column_grid",
    ] + ["feature_{}".format(i) for i in range(num_features)]

    output = pd.DataFrame(columns=headers)

    # Generate data
    for i in range(num_slices):
        data = []
        slice_id = "slice_{}".format(i)
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
                features = [i / np.sum(features) for i in features]
                data.append(
                    [slice_id, is_tissue, cluster, row_grid, column_grid] + features
                )
        table = pd.DataFrame(data, columns=headers)

        # Convert grid coordinates to cartesian coordinates
        v_1 = np.array([point_distance, 0])
        v_2 = np.array(
            [point_distance * np.cos(np.pi / 3), point_distance * np.sin(np.pi / 3)]
        )

        # convert coords to numpy array
        grid_coords = table[["column_grid", "row_grid"]].to_numpy()
        base = np.array([v_1, v_2])
        cartesian_coords = np.dot(grid_coords, base)
        ones = np.ones((cartesian_coords.shape[0], 1))
        cartesian_coords = np.concatenate((cartesian_coords, ones), axis=1)

        # rotate the slices around a random pivot
        if rotation_degree_range:
            rotation = np.random.uniform(-rotation_degree_range, rotation_degree_range)
            pivot = np.array(
                [
                    np.random.randint(-slice_diameter, slice_diameter),
                    np.random.randint(-slice_diameter, slice_diameter),
                ]
            )
            translation_matrix = np.array(
                [[1, 0, -pivot[0]], [0, 1, -pivot[1]], [0, 0, 1]]
            )
            rotattion_matrix = np.array(
                [
                    [np.cos(rotation), -np.sin(rotation), 0],
                    [np.sin(rotation), np.cos(rotation), 0],
                    [0, 0, 1],
                ]
            )
            translation_back_matrix = np.array(
                [[1, 0, pivot[0]], [0, 1, pivot[1]], [0, 0, 1]]
            )
            transformation_matrix = np.dot(
                translation_back_matrix, np.dot(rotattion_matrix, translation_matrix)
            )
            cartesian_coords = np.dot(cartesian_coords, transformation_matrix)

        # append the cartesian coordinates to the table
        table["column"] = cartesian_coords[:, 0]
        table["row"] = cartesian_coords[:, 1]

        output = pd.concat([output, table])

    os.makedirs(os.path.dirname(path), exist_ok=True)
    output.to_csv(f"{path}data_{num_points}.tsv", index=False, sep="\t")

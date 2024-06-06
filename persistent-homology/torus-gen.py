import random
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay

# torus parameters
R = 1
r = 0.5
num_points = 1000


width = 2 * (R + r)
height = 2 * r


# test if a point is inside the torus
def torus(x, y, z):
    return (R - np.sqrt(x**2 + y**2)) ** 2 + z**2 < r**2


# generate random points in the bounding box of the torus, with a random value appended
points = np.array([p for p in np.random.rand(num_points, 3) * [width, height, width] - [R + r, r, R + r] if torus(*p)])

# delauany tetrahedralization
tets_data = pd.DataFrame(Delaunay(points).simplices)

# generate random values for each point
values = pd.Series(np.random.rand(len(points)))
points_data = pd.concat([pd.DataFrame(points), values], axis=1)


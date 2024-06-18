import dionysus as d
import numpy as np
from scipy.spatial import Delaunay
from itertools import combinations
from scipy.spatial.distance import euclidean
import miniball
import re
import sys
import mrcfile
from PIL import Image
import scipy.sparse as sp

def add_perturbation(points, magnitude=0.001):
    perturbation = np.random.uniform(-magnitude, magnitude, size=np.array(points).shape)
    return points + perturbation

if len(sys.argv) < 6:
    print("Usage: python main.py input.txt output.txt diff_threshold min_alpha_diff output_verts.txt")
    sys.exit(1)

file_path = sys.argv[1]

data = []

def in_range(low, high, point, data):
    x,y,z = point[0],point[1],point[2]

    if low <= data[x][y][z] <= high:
        return True
    return False

# Open the image
image = Image.open(file_path)

# Convert the image to grayscale
image = image.convert('L')

# Convert the image to a numpy array
data = np.array(image)

print("Data shape:", data.shape)

shape = data.shape

all_points = []
all_edges = []
all_tris = []

points_alpha = []
edges_alpha = []
tris_alpha = []

lower = float(sys.argv[3])
upper = float(sys.argv[4])

points_index = [[-1 for _ in range(data.shape[1])] for _ in range(data.shape[0])]

#max_data = np.max(data)
max_data = 255

#verts
for i in range(data.shape[0]):
    for j in range(data.shape[1]):

        v = np.array([i, j])

        points_index[v[0]][v[1]] = len(all_points)
        all_points.append(v)
        #points_alpha.append(int(data[v[0],v[1]]))
        points_alpha.append(int(max_data) - int(data[v[0],v[1]]))

print("Vertices created")

edge_offset = [[1,0],[0,1]]

for i in range(data.shape[0]):
    for j in range(data.shape[1]):
            
        curr_v = np.array([i, j])

        for offset in edge_offset:
            other_v = curr_v + np.array(offset)

            if not np.any(other_v >= shape): 
                edge = [points_index[curr_v[0]][curr_v[1]],points_index[other_v[0]][other_v[1]]]
                all_edges.append(edge)
                edges_alpha.append(max(points_alpha[edge[0]],points_alpha[edge[1]]))

print("Edges created")

square_verts = [[0,0],[1,0],[0,1],[1,1]]

for i in range(data.shape[0]):
    for j in range(data.shape[1]):

        curr_v = np.array([i, j])

        v1, v2, v3, v4 = curr_v+np.array(square_verts[0]), curr_v+np.array(square_verts[1]), curr_v+np.array(square_verts[2]), curr_v+np.array(square_verts[3])

        if not np.any(v4 >= shape): 

            square = [points_index[v1[0]][v1[1]],points_index[v2[0]][v2[1]],points_index[v3[0]][v3[1]],points_index[v4[0]][v4[1]]]
            #diagonal of square
            diagonal = [square[i] for i in [0, 3]]
            #the two triangle faces
            face1 = [square[i] for i in [0, 1, 3]]
            face2 = [square[i] for i in [0, 2, 3]]

            alpha = max(points_alpha[square[0]],points_alpha[square[1]],points_alpha[square[2]],points_alpha[square[3]])
            
            all_edges.append(diagonal)
            edges_alpha.append(alpha)

            all_tris.append(face1)
            all_tris.append(face2)
            tris_alpha.append(alpha)
            tris_alpha.append(alpha)
            
print("Faces and square diagonal created")
print(len(all_points))
print(len(all_edges))
print(len(all_tris))

f = d.Filtration()

for i in range(len(all_points)):
    f.append(d.Simplex([i],points_alpha[i]))
for i in range(len(all_edges)):
    f.append(d.Simplex(all_edges[i],edges_alpha[i]))
for i in range(len(all_tris)):
    f.append(d.Simplex(all_tris[i],tris_alpha[i]))

#filtration complete
f.sort()

print("Starting obtaining PH")

m = d.homology_persistence(f,method = 'column')


print("Starting writing to file")

with open(sys.argv[2], 'w') as file:
    
    for i,c in enumerate(m):
        numbers_after_asterisk = re.findall(r'\*(\d+)', str(c))
        # Convert the matched numbers to integers
        numbers_after_asterisk = [int(num)+1 for num in numbers_after_asterisk]
        result_line = 'Column '+' '.join(map(str, numbers_after_asterisk))+ '\n'
        file.write(result_line)

    for s in f:
        matches = re.findall(r'<(.*?)>', str(s))[0]
        verts = matches.split(',')
        result_line = 'Simplex '+' '.join(map(str, verts))+ ' ' + str(s.data) + '\n'
        file.write(result_line)

    count = 0

    for i in range(len(m)):
        if m.pair(i) < i: continue  # skip negative simplices
        
        if m.pair(i) != m.unpaired:
            if(f[m.pair(i)].data - f[i].data) >= float(sys.argv[5]) and f[i].dimension() == 1:
                result_line = 'Interval '+ str(i+1)+ ' ' + str(m.pair(i)+1) + ' '+str(f[i].data) +' ' +str(f[m.pair(i)].data) + '\n'
                file.write(result_line)
                count+= 1

    print(count)

with open(sys.argv[6], 'w') as f:

    for m in range(len(all_points)):
        point = all_points[m]
        value = points_alpha[m]
        i, j = point[0], point[1]
        f.write(f"{i} {j} {value}\n")

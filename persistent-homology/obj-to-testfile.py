
# load the obj file
from asyncio import sleep


vertices = []

with open('persistent-homology/assets/torus.obj', 'r') as f:
    for line in f:
        line = line.split()
        if line[0] == 'v':
            vertices.append((float(line[1]), float(line[2]), float(line[3])))
        elif line[0] == 'f':
            face = [int(i) for i in line[1:]]
            triangles = []
            for i in range(len(face) - 2):
                triangles.append((face[0], face[i+1], face[i+2]))
                pass

sleep(1)
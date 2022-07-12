from objects import Vec3f, number_on_vector_mult
from primitives import Triangle

def read(filename):
    verts = []
    faces = []
    scale = 1
    pos = ''
    
    with open(filename, "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i].strip().split()
            if i == 0:
                pos = Vec3f(float(line[0]), float(line[1]), float(line[2]))
                scale = float(line[3])
            else:
                if line[0] == "v":
                    v = Vec3f(float(line[1]), float(line[2]), float(line[3]))
                    v = number_on_vector_mult(scale, v) + pos
                    verts.append(v)

                if line[0] == "f":
                    faces.append([int(line[1]),int(line[2]),int(line[3])])

    return verts, faces
                    
def load_model(filename, material):            
    v, f = read(filename)
    objects = []
    
    for i in range(len(f)):
        id_1 = f[i][0] - 1
        id_2 = f[i][1] - 1
        id_3 = f[i][2] - 1
        objects.append(Triangle(v[id_1], v[id_2], v[id_3], material))

    return objects

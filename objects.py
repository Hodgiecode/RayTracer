from abc import ABC, abstractmethod
import math

class Vec3f:
    def __init__(self, x = 1, y = 1, z = 1):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Vec3f(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vec3f(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        if not isinstance(other, int) and not isinstance(other, float):
            return self.x * other.x + self.y * other.y + self.z * other.z
        else:
            return Vec3f(self.x * other, self.y * other, self.z * other)
            

    def length(self):
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalize(self):
        leng = self.length()
        if leng > 0:
           self.x = self.x / leng
           self.y = self.y / leng
           self.z = self.z / leng


class Light:
    def __init__(self, position, intensity):
        self.position = position
        self.intensity = intensity


class Object(ABC):
    def __init__(self, material):
        self.material = material
        self.normal = ''
        self.intersection_point = ''
        
    def get_material(self):
        return self.material

    def get_normal(self):
        return self.normal

    def get_intersection_point(self):
        return self.intersection_point
    
        
    @abstractmethod
    def check_intersect(self, ray_src, ray_dir):
        pass

    @abstractmethod
    def isIntersect(self, ray_src, ray_dir):
        pass

def number_on_vector_mult(num, vect):
    return Vec3f(num * vect.x, num * vect.y, num * vect.z)

def crossProduct(vect_A, vect_B):
    p_x = vect_A.y * vect_B.z - vect_A.z * vect_B.y
    p_y = vect_A.z * vect_B.x - vect_A.x * vect_B.z
    p_z = vect_A.x * vect_B.y - vect_A.y * vect_B.x
    return Vec3f(p_x, p_y, p_z)

def transpose_vector(vector):
    result = []
    for i in range(len(vector)):
        result.append([vector[i]])

    return result
    
def matrix_mult(m1, m2):
    result = []
    for i in range(len(m1)):
       tmp = []
       for j in range(len(m2[0])):
           tmp.append(0)

       result.append(tmp)

    for i in range(len(m1)):
       for j in range(len(m2[0])):
           for k in range(len(m2)):
               result[i][j] += m1[i][k] * m2[k][j]

    return result


def rotation(angle, vect):
    vector = [vect.x, vect.y, vect.z]
    
    angle = angle * math.pi / 180
    m_x = [[1,0,0],
           [0,math.cos(angle),-math.sin(angle)],
           [0,math.sin(angle), math.cos(angle)]
           ]
    
    a = matrix_mult(m_x, transpose_vector(vector))
    b = Vec3f(a[0][0], a[1][0], a[2][0])
    return b

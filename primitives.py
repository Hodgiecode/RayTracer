from objects import Object, Vec3f, crossProduct, number_on_vector_mult
from solver import Solver
import math
from numpy.polynomial import polynomial as P
import numpy as np

### Sphere

class Sphere(Object):
    def __init__(self, center, radius, material):
        self.center = center
        self.radius = radius
        self.material = material
        self.tau = -1
        self.eps = 0.001

        
    def check_intersect(self, ray_src, ray_dir):
        L = self.center - ray_src
        tca = L * ray_dir
        d2 = L * L - tca * tca

        if d2 > self.radius ** 2:
            return False

        thc = math.sqrt(self.radius ** 2 - d2)

        t0 = tca - thc
        t1 = tca + thc

        if t0 > self.eps:
            self.tau = t0
            return True

        if t1 > self.eps:
            self.tau = t1
            return True

        return False

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal = self.intersection_point - self.center
            self.normal.normalize()
            return True

        return False

### Triangle
class Triangle(Object):
    def __init__(self, A, B, C, material):
        self.A = A
        self.B = B
        self.C = C
        self.material = material
        self.tau = -1
        self.eps = 0.00001

    def check_intersect(self, ray_src, ray_dir):
        e1 = self.B - self.A
        e2 = self.C - self.A
        p = crossProduct(ray_dir, e2)

        det = e1 * p
        if det < self.eps:
            return False

        tvec = ray_src - self.A
        u = tvec * p

        if u < 0 or u > det:
            return False

        qvec = crossProduct(tvec, e1)
        v = ray_dir * qvec

        if v < 0 or u + v > det:
            return False
            
        self.tau = e2 * qvec * (1.0 / det)
        
        if self.tau < self.eps:
            return False
            
        return True

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal = crossProduct(self.B - self.A, self.C - self.A)
            self.normal.normalize()
            return True
        
        return False

#Cube
class Cube(Object):
    def __init__(self, center, side, material):
        self.center = center
        self.side = side
        self.material = material
        self.tau = -1
        self.eps = 0.001

    def check_intersect(self, ray_src, ray_dir):
        normals = [ Vec3f(1.0, 0.0, 0.0), Vec3f(-1.0, 0.0, 0.0),
                    Vec3f(0.0, 1.0, 0.0), Vec3f(0.0, -1.0, 0.0),
                    Vec3f(0.0, 0.0, 1.0), Vec3f(0.0, 0.0, -1.0)]

        t_min = - 1.0
        int_normal = ''
        o_to_c = self.center - ray_src
        
        for i in range(6):
            n = normals[i]
            n_dot_d = n * ray_dir
            if n_dot_d >= 0.0:
                continue

            t = (n * o_to_c + self.side) / n_dot_d

            if t > t_min and t_min >= 0.0:
                continue

            diff = ray_src + ray_dir * t - self.center
            if max(math.fabs(diff.x), math.fabs(diff.y), math.fabs(diff.z)) < self.side + self.eps:
                t_min = t
                int_normal = n

        if t_min < 0.0:
            return False

        if int_normal == '':
            return False

        self.tau = t_min
        self.normal = int_normal
        return True
                

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal.normalize()
            return True
        
        return False

## Cone
class Cone(Object):
    def __init__(self, center, height, radius, material):
        self.center = center
        self.height = height
        self.radius = radius
        self.material = material
        self.tau = -1
        self.eps = 1e-9

    def check_intersect(self, ray_src, ray_dir):
        A = ray_src.x - self.center.x
        B = ray_src.z - self.center.z
        D = (-1) * self.height + ray_src.y - self.center.y

        rh = float(self.radius)/self.height
        tan = rh * rh

        a = (ray_dir.x * ray_dir.x) + (ray_dir.z * ray_dir.z) - (tan * (ray_dir.y * ray_dir.y))
        b = (2 * A * ray_dir.x) + (2 * B * ray_dir.z) - (2 * tan * D * ray_dir.y)
        c = (A * A) + (B * B) - (tan * (self.height - ray_src.y + self.center.y) *
                                 (self.height - ray_src.y + self.center.y))

        delta = b * b - 4 * a * c
        if math.fabs(delta) < self.eps:
            return False

        if delta < self.eps:
            return False

        t1 = (-b - math.sqrt(delta))/(2*a)
        t2 = (-b + math.sqrt(delta))/(2*a)

        if (t1 > self.eps):
            self.tau = t1
            r = ray_src.y + self.tau * ray_dir.y
            if ((r > self.center.y) and (r < self.center.y + self.height)):
                return True
        else:
            if (t2 > self.eps):
                self.tau = t2
                r = ray_src.y + self.tau * ray_dir.y
                if ((r > self.center.y) and (r < self.center.y + self.height)):
                    return True
           
        return False

    def isIntersect(self, ray_src, ray_dir):
        top = Vec3f(self.center.x, self.center.y + self.height, self.center.z)
        
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            r_ = math.sqrt((self.intersection_point.x - self.center.x) * (self.intersection_point.x - self.center.x) + (self.intersection_point.z - self.center.z) * (self.intersection_point.z - self.center.z));
            self.normal = Vec3f(self.intersection_point.x - self.center.x, r_*(float(self.radius)/self.height), self.intersection_point.z - self.center.z)
            self.normal.normalize()
            return True
        
        return False


## Cylinder
class Cylinder(Object):
    def __init__(self, center, height, radius, material):
        self.center = center
        self.height = height   
        self.radius = radius
        
        self.normal = ''
        self.material = material
        self.tau = -1
        self.eps = 1e-9

    
    def check_intersect(self, ray_src, ray_dir):
        a = (ray_dir.x * ray_dir.x) + (ray_dir.z * ray_dir.z)
        b = 2 * (ray_dir.x * (ray_src.x - self.center.x) + ray_dir.z * (ray_src.z - self.center.z))
        c = (ray_src.x - self.center.x) * (ray_src.x - self.center.x) + (ray_src.z - self.center.z) * (ray_src.z - self.center.z) - (self.radius * self.radius)

        delta = b * b - 4 * a * c

        if math.fabs(delta) < self.eps:
            return False

        if delta < self.eps:
            return False

        t1 = (-b - math.sqrt(delta))/(2*a)
        t2 = (-b + math.sqrt(delta))/(2*a)

        if (t1 > self.eps):
            self.tau = t1
            r = ray_src.y + self.tau * ray_dir.y
            if ((r > self.center.y) and (r < self.center.y + self.height)):
                return True
        else:
            if (t2 > self.eps):
                self.tau = t2
                r = ray_src.y + self.tau * ray_dir.y
                if ((r > self.center.y) and (r < self.center.y + self.height)):
                    return True
           
        return False


    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal = Vec3f(self.intersection_point.x - self.center.x, 0, self.intersection_point.z - self.center.z) 
            self.normal.normalize()
            return True

        return False

#Эллипсоид
class Ellipsoid(Object):
    def __init__(self, center, a, b, c, material):
        self.center = center
        self.a = a
        self.b = b
        self.c = c
        self.material = material
        self.normal = None
        self.tau = -1
        self.eps = 1e-6

    def check_intersect(self, ray_src, ray_dir):
        axis = Vec3f(self.a, self.b, self.c)
        axis_dir = Vec3f(float(ray_dir.x) / self.a, float(ray_dir.y) / self.b, float(ray_dir.z) / self.c)
        axis_src = Vec3f(float(ray_src.x) / self.a, float(ray_src.y) / self.b, float(ray_src.z) / self.c)
        axis_center = Vec3f(float(self.center.x) / self.a, float(self.center.y) / self.b, float(self.center.z) / self.c)

        temp = math.sqrt(axis_dir.x * axis_dir.x + axis_dir.y * axis_dir.y + axis_dir.z * axis_dir.z)
        A =  temp * temp
        B = 2 * (axis_dir * (axis_src - axis_center))

        tmp = axis_src - axis_center;
        temp = math.sqrt(tmp.x * tmp.x + tmp.y * tmp.y + tmp.z * tmp.z)
        C = temp * temp - 1
        disc = (B * B) - (4 * A * C)
        
        if (disc < 0):
           return False

        t1 = (-B - math.sqrt(disc))/(2 * A);
        t2 = (-B + math.sqrt(disc))/(2 * A);
    
        if t1 > self.eps:
            self.tau = t1
            return True

        if t2 > self.eps:
            self.tau = t2
            return True
    
        return False

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal = self.intersection_point - self.center
            self.normal = Vec3f(2 * self.normal.x / (self.a * self.a),
                                2 * self.normal.y / (self.b * self.b),
                                2 * self.normal.z / (self.c * self.c))
            self.normal.normalize()
            return True

        return False
    

#Парабалоид
class Paraboloid(Object):
    def __init__(self, center, size, max_y, material):
        self.center = center
        self.max_y = max_y
        self.size = size
        
        self.material = material
        self.normal = None
        self.tau = -1
        self.eps = 1e-6

    def check_intersect(self, ray_src, ray_dir):
        dist = ray_src - self.center
        A =   ray_dir.x * ray_dir.x + ray_dir.z * ray_dir.z
        B = 2 * ((dist.x * ray_dir.x) + (dist.z * ray_dir.z)) - ray_dir.y
        C =  (dist.x * dist.x) + (dist.z * dist.z) - (dist.y + self.size)

        disc = B * B - (4 * A * C)
        if (disc < 0):
           return False

        t1 = (-B - math.sqrt(disc))/(2 * A);
        t2 = (-B + math.sqrt(disc))/(2 * A);
    
        if t1 > self.eps:
            self.tau = t1
            return True

        if t2 > self.eps:
            self.tau = t2
            return True
    
        return False

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau

            if (self.intersection_point.y < self.max_y):
                r = self.intersection_point - self.center;

                self.normal = Vec3f(2 * r.x, -1, 2 * r.z)
                self.normal.normalize();

                return True
    
        return False
    
#Гиперболоид
class Hyperboloid(Object):
    def __init__(self, center, size, max_y, min_y, type, material):
        self.center = center
        self.max_y = max_y
        self.min_y = min_y
        self.size = size
        self.type = type
        self.material = material
        self.normal = None
        self.tau = -1
        self.eps = 1e-6

    def check_intersect(self, ray_src, ray_dir):
        dist = ray_src - self.center
        A = (ray_dir.x * ray_dir.x) + (ray_dir.z * ray_dir.z) - (ray_dir.y * ray_dir.y)
        B = 2.0 * ((dist.x * ray_dir.x) + (dist.z * ray_dir.z) - (dist.y * ray_dir.y))
        C = 0.0

        if (self.type == 0):
            C = (dist.x * dist.x) + (dist.z * dist.z) - ((dist.y * dist.y) + self.size)
        else:
            C = (dist.x * dist.x) + (dist.z * dist.z) - ((dist.y * dist.y) - self.size)
            
        disc = B * B - (4 * A * C)
        if (disc < 0):
           return False

        t1 = (-B - math.sqrt(disc))/(2 * A);
        t2 = (-B + math.sqrt(disc))/(2 * A);
    
        if t1 > self.eps:
            self.tau = t1
            return True

        if t2 > self.eps:
            self.tau = t2
            return True
    
        return False

    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau

            if (self.intersection_point.y > self.min_y and self.intersection_point.y < self.max_y):
                r = self.intersection_point - self.center;

                self.normal = Vec3f(2 * r.x, -2 * r.y, 2 * r.z)
                self.normal.normalize();

                return True
    
        return False
        
#Torus
class Torus(Object):   
    def __init__(self, center, r1, r2, axis, material):
        self.center = center
        self.r = r1
        self.R = r2
        self.axis = axis
        self.material = material
        
        self.normal = None
        self.tau = -1
        self.eps = 1e-6

    def check_intersect(self, ray_src, ray_dir):
        self.axis.normalize()
        
        e = ray_src - self.center
        d = ray_dir

        d.normalize()

        iR = self.R * self.R - self.r * self.r
        oR = 4 * self.R * self.R
        L = (e.x * e.x + e.y * e.y + e.z * e.z) + iR
        e_dot_d = e * d

        c4 = (d.x * d.x + d.y * d.y + d.z * d.z) * (d.x * d.x + d.y * d.y + d.z * d.z)
        c3 = 4 * (d.x * d.x + d.y * d.y + d.z * d.z) * e_dot_d
        c2 = 2 * (2 * ((e_dot_d) * (e_dot_d)) + (d.x * d.x + d.y * d.y + d.z * d.z) * L) - oR * ((d.x * d.x) + (d.y * d.y))
        c1 = 4 * e_dot_d * L - 2 * oR * (e.x * d.x + e.y * d.y)
        c0 = (L * L) - oR * ((e.x * e.x) + (e.y * e.y))
		
        solver = Solver()
        num_roots = solver.solve_p4(c3/c4, c2/c4, c1/c4, c0/c4)
        
        if (num_roots == 0):
          return False

        rts = solver.x
        
        intersect = False
        min_t = 9999999
        
        for i in range(num_roots):
            t = rts[i]
            if (t > 1e-5 and t < min_t):
                min_t = t
                intersect = True
        
        if intersect == False:
            return False
        
        self.tau = min_t
        return True
    
    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            p = self.intersection_point - self.center;

            cos_theta = (self.axis * p)/(math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z));
            sin_theta = math.sqrt(1 - cos_theta * cos_theta);

            p.normalize();

            self.normal = p - number_on_vector_mult((1.0 / sin_theta), number_on_vector_mult(self.R , (p - number_on_vector_mult(cos_theta, self.axis))))
            self.normal = number_on_vector_mult((1.0 / self.r), self.normal)
            self.normal.normalize()
            return True
			
        return False

        
    
           

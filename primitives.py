from objects import Object, Vec3f, crossProduct, number_on_vector_mult
import math

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
        self.eps = 0.001

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

        qvec = crossProduct(tvec,e1)
        v = ray_dir * qvec

        if v<0 or u + v > det:
            return False

        self.tau = e2 * qvec * (1.0 / det)
        return self.tau > self.eps

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
        self.eps = 0.001

    def check_intersect(self, ray_src, ray_dir):
        top = Vec3f(self.center.x, self.center.y + self.height, self.center.z)
        rh = float(self.radius)/self.height
        
        a = (ray_dir.x * ray_dir.x) + (ray_dir.z * ray_dir.z) - rh * rh * (ray_dir.y * ray_dir.y)

        b = 2 *((ray_src.x - self.center.x) * ray_dir.x + (ray_src.z - self.center.z) * ray_dir.z -
                   rh * rh *((-1)*self.height + ray_src.y - self.center.y) * ray_dir.y)

        c = -rh * rh * (self.height - ray_src.y + self.center.y) * (self.height - ray_src.y + self.center.y) + (ray_src.x - self.center.x) * (ray_src.x - self.center.x) + (ray_src.z - self.center.z) * (ray_src.z - self.center.z)

        delta = b * b - (4 * a * c);

        if (math.fabs(delta) < 1e-9):
            return False

        if delta < 0:
            return False
        
        t1 = (-b - math.sqrt(delta))/(2*a)
        t2 = (-b + math.sqrt(delta))/(2*a)

        pt = ray_src + number_on_vector_mult(t1, ray_dir)
        if (pt.y >= self.center.y and pt.y <= self.center.y + self.height):
            self.tau = t1
            return True
    
        pt = ray_src + number_on_vector_mult(t2, ray_dir)
        if (pt.y >= self.center.y and pt.y <= self.center.y + self.height):
            self.tau = t2
            return True

        return False

    def isIntersect(self, ray_src, ray_dir):
        top = Vec3f(self.center.x, self.center.y + self.height, self.center.z)
        
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            r_ = math.sqrt((self.intersection_point.x - self.center.x) * (self.intersection_point.x - self.center.x) + (self.intersection_point.z - self.center.z) * (self.intersection_point.z - self.center.z));
            self.normal = Vec3f(self.intersection_point.x - self.center.x, r_*(self.radius/self.height), self.intersection_point.z - self.center.z)
            self.normal.normalize()
            return True
        
        return False


## Cylinder
class Cylinder(Object):
    def __init__(self, center, height, radius, material):
        0


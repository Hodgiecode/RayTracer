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

        if (t1 > t2):
            self.tau = t2
        else:
            self.tau = t1

        if self.tau < 0:
            return False

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
    def __init__(self, center_top, center_bottom, radius, material):
        self.center_top = center_top
        self.center_bottom = center_bottom
        self.center = Vec3f(0.5 * (self.center_top.x - self.center_bottom.x),
                            0.5 * (self.center_top.y - self.center_bottom.y),
                            0.5 * (self.center_top.z - self.center_bottom.z))

        self.height = (self.center_top.y - self.center_bottom.y)
            
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

        if (t1 > t2):
            self.tau = t2
        else:
            self.tau = t1

        #if self.tau < 0:
            #return False

        r = ray_src.y + self.tau * ray_dir.y
        
        if ((r >= self.center.y) and (r <= self.center.y + self.height)):
            return True

        return False
        
        '''
        h = self.center_top - self.center_bottom
        a_bottom = ray_src - self.center_bottom
        a_bottom_h = crossProduct(a_bottom, h)
        ray_dir_h = crossProduct(ray_dir, h)

        h2 = h * h
        a = ray_dir_h * ray_dir_h
        b = number_on_vector_mult(2, ray_dir_h) * a_bottom_h
        c = a_bottom_h * a_bottom_h - (self.radius * self.radius - h2)
        d = b * b - 4 * a * c

        t1 = 0
        t2 = 0
        ip = ''

        if (d > 0):
            t1 = (-b - math.sqrt(d))/(2*a)
            t2 = (-b + math.sqrt(d))/(2*a)

            if (t1 > t2):
                self.tau = t2
            else:
                self.tau = t1

            if (self.tau > 0):
                ip = ray_src + number_on_vector_mult(self.tau, ray_dir)
                h_ = self.center_top - self.center_bottom
                v = h_ * (ip - self.center_bottom) > 0 and h_ * (ip - self.center_top) < 0
                if not v:
                    self.tau = -1
            else:
                return False

        else:
            return False
            

        v = ip - self.center_bottom
        tmp = number_on_vector_mult(h, (number_on_vector_mult(h, v)))
        tmp.x = tmp.x / (h*h)
        tmp.y = tmp.y / (h*h)
        tmp.z = tmp.z / (h*h)

        normal_ = v - tmp;

        ncap = self.center_top - self.center_bottom;
        tcap = self.cap_hit(self.center_bottom, ncap, ray_src, ray_dir, self.radius);

        if ((tcap >= 0 and tcap < self.tau) or self.tau < 0):
            self.tau = tcap
            normal_ = ncap
    

        ncap = (self.center_bottom - self.center_top)
        tcap = self.cap_hit(self.center_top, ncap, ray_src, ray_dir, self.radius)

        if ((tcap >= 0 and tcap < self.tau) or self.tau < 0):
            self.tau = tcap
            self.normal = ncap

        if (self.tau < 0):
            return False
        '''
        return True


    def isIntersect(self, ray_src, ray_dir):
        if self.check_intersect(ray_src, ray_dir):
            self.intersection_point = ray_src +  ray_dir * self.tau
            self.normal = Vec3f(self.intersection_point.x - self.center.x,
                                0,
                                self.intersection_point.z - self.center.z)
                                
            self.normal.normalize()
            return True

        return False


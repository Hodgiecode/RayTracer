from primitives import Sphere, Cube, Cone, Cylinder, Ellipsoid, Paraboloid, Hyperboloid
from objects import Vec3f, Light, number_on_vector_mult
from materials import Red, Glass, Ivory
import math
from PIL import Image
from loads import load_model

def reflect(i, n):
    return i - n * 2.0 * (i * n)

def refract(i, n, eta_t, eta_i = 1.0):
    cosi = (-1) * max(-1.0, min(1.0, i * n))

    if cosi < 0:
        return refract(i, number_on_vector_mult((-1),n), eta_i, eta_t)

    eta = eta_i / eta_t
    k = 1 - eta * eta * (1 - cosi * cosi)

    if k < 0:
        return Vec3f(1, 0, 0)
    else:
        return i * eta + n * (eta * cosi - math.sqrt(k))

def cast_ray(ray_src, ray_dir, objects, lights, depth = 0):
    if depth > 4:
        return Vec3f(0.2, 0.7, 0.8)

    intersect_flag = False
    m = ''
    n = ''
    pt = ''
    
    for i in range(len(objects)):
        if objects[i].isIntersect(ray_src, ray_dir):
            m = objects[i].get_material()
            n = objects[i].get_normal()
            pt = objects[i].get_intersection_point()
            intersect_flag = True

    if intersect_flag == False:
        return Vec3f(0.2, 0.7, 0.8)

    reflect_dir = reflect(ray_dir, n)
    reflect_dir.normalize()

    refract_dir = refract(ray_dir, n, m.refractive_index)
    refract_dir.normalize()

    reflect_color = cast_ray(pt, reflect_dir, objects, lights, depth + 1)
    refract_color = cast_ray(pt, refract_dir, objects, lights, depth + 1)

    diffuse_light_intensity = 0
    specular_light_intensity = 0

    for i in range(len(lights)):
        light_dir = lights[i].position - pt
        light_dir.normalize()

        shadow_material = ''
        shadow_n = ''
        shadow_pt = ''
        shadow_intersect = False
        
        for q in range(len(objects)):
            if objects[q].isIntersect(pt, light_dir):
                shadow_mat = objects[q].get_material()
                shadow_pt = objects[q].get_intersection_point()
                shadow_n = objects[q].get_normal()
                shadow_intersect = True
                break

        if shadow_intersect == True:
            t1 = shadow_pt - pt
            t2 = lights[i].position - pt

            if math.sqrt(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z) < math.sqrt(t2.x * t2.x + t2.y * t2.y + t2.z * t2.z):
                continue

        diffuse_light_intensity  = diffuse_light_intensity + lights[i].intensity * max(0.0, math.fabs(light_dir * n))
        
        p1 = number_on_vector_mult(-1,light_dir)
        p2 = reflect(p1, n)
        p3 = number_on_vector_mult(-1, p2) * ray_dir

        specular_light_intensity = specular_light_intensity + pow(max(0.0, p3), m.specular) * lights[i].intensity;

    return m.color * diffuse_light_intensity * m.albedo_p1 + Vec3f(1.0, 1.0, 1.0) * specular_light_intensity * m.albedo_p2 + reflect_color * m.albedo_p3 + refract_color * m.albedo_p4
        

    

class Scene:
    def __init__(self, W, H):
        self.objects = []
        self.light_vector = []

        light_1 = Light(Vec3f(-20, 20, 20), 1.5)
        light_2 = Light(Vec3f(30, 50, -25), 1.8)
        light_3 = Light(Vec3f(30, 20, 30), 1.7)
        
        self.light_vector = [light_1, light_2, light_3]

        mat_1 = Ivory()
        mat_2 = Red()
        mat_3 = Glass()
        
        #self.objects = load_model("bunny.txt", mat)
        #s1 = Sphere(Vec3f(-3, 0, -16), 5, mat_3)
        #s2 = Sphere(Vec3f(0, 0, -16), 5, mat_2)
        #s3 = Sphere(Vec3f(3, 3, -16), 5, mat_3)

        #self.objects.append(Cone(Vec3f(0, -5.0, -16), 10, 5, mat_1));
        #self.objects.append(s1)
        #self.objects.append(s2)
        #self.objects.append(s3)

        #self.objects.append(Cylinder(Vec3f(0, -5, -16), 10, 3, mat_1))
        #self.objects.append(Cone(Vec3f(0, -5.0, -16), 10, 5, mat_2));
        #self.objects.append(Cube(Vec3f(-5, 0, -16), 3, mat_1));
        #self.objects.append(Cube(Vec3f(5, 0, -16), 3, mat_2));

        #self.objects.append(Paraboloid(Vec3f(0, 0, -16), 3, 5, mat_3))

        self.objects.append(Ellipsoid(Vec3f(0, 0, -16), 8, 5, 3, mat_3))

        #self.objects.append(Hyperboloid(Vec3f(-5,0,-16), 3, 5, -5, 0, mat_1))
        #self.objects.append(Hyperboloid(Vec3f(5,0,-16), 3, 5, -5, 1, mat_1))
        

    def get_pixel(self, x, y, w, h):
        fov = 3.14 / 3.0;
        dir_x = (x + 0.5) - w / 2.0;
        dir_y = -(y + 0.5) + h / 2.0;
        dir_z =  -h /(2.0 * math.tan(fov / 2.0));

        ray_dir = Vec3f(dir_x, dir_y, dir_z);
        ray_dir.normalize();

        c = cast_ray(Vec3f(0,0,0), ray_dir, self.objects, self.light_vector)

        mx = max(c.x, max(c.y, c.z))
        
        if mx > 1:
            c = c * (1.0 / mx)

        c.x = int(255 * max(0.0, min(1.0, c.x)))
        c.y = int(255 * max(0.0, min(1.0, c.y)))
        c.z = int(255 * max(0.0, min(1.0, c.z)))

        return c
        
        
        
###
W = 200
H = 100
img = Image.new( mode = "RGB", size = (W, H) )
data = []
scene = Scene(W, H)
for y in range(H):
      for x in range(W):
          col = scene.get_pixel(x, y, W, H)
          img.putpixel((x,y), (col.x, col.y, col.z))

print('finished')
#img.show()
img.save("my.png")




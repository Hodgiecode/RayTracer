from objects import Vec3f

class Material:
    def __init__(self, refractive_index = 1, albedo_p1 = 1, albedo_p2 = 0, albedo_p3 = 0, albedo_p4 = 0, specular = 0, color = Vec3f(0.3, 0.1,0.1)):
        self.refractive_index = refractive_index
        self.albedo_p1 = albedo_p1
        self.albedo_p2 = albedo_p2
        self.albedo_p3 = albedo_p3
        self.albedo_p4 = albedo_p4
        self.specular = specular
        self.color = color

class Red(Material):
    def __init__(self):
        self.refractive_index = 1.0
        self.albedo_p1 = 1.4
        self.albedo_p2 = 0.3
        self.albedo_p3 = 0.0
        self.albedo_p4 = 0.0
        self.color = Vec3f(0.3,0.1,0.1)
        self.specular = 10.0


class Glass(Material):
    def __init__(self):
        self.refractive_index = 1.5
        self.albedo_p1 = 0.0
        self.albedo_p2 = 0.9
        self.albedo_p3 = 0.1
        self.albedo_p4 = 0.8
        self.color = Vec3f(0.6,0.7,0.8)
        self.specular = 125.0

class Ivory(Material):
    def __init__(self):
        self.refractive_index = 1.0
        self.albedo_p1 = 0.6
        self.albedo_p2 = 0.3
        self.albedo_p3 = 0.1
        self.albedo_p4 = 0.0
        self.color = Vec3f(0.4, 0.4, 0.3)
        self.specular = 50.0

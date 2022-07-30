import math

class Solver:
    def __init__(self):
        self.two_pi = 6.28318530717958648
        self.eps = 1e-14
        self.x = []

    def solve_p3(self, a, b, c):
        a2 = a * a
        q = (a2 - 3 * b)/ 9.0
        r = (a * (2 * a2 - 9 * b) + 27 * c)/ 54.0
        r2 = r * r
        q3 = q * q * q
        A = 0
        B = 0

        if r2 < q3:
            t = r / math.sqrt(q3)
            if t < -1:
                t = -1
            if t > 1:
                t = 1

            t = math.acos(t)
            a = a/3.0
            q = -2 * math.sqrt(q)
            self.x.append(q * math.cos(t/3.0) - a)
            self.x.append(q * math.cos((t + self.two_pi) / 3.0) - a)
            self.x.append(q * math.cos((t - self.two_pi) / 3.0) - a)
            return 3
        else:
            A = - ((math.fabs(r) + math.sqrt(r2 - q3)) ** (1.0/3.0))
            if r < 0:
                A = - A

            if A == 0:
                B = 0
            else:
                B = q/A

            a = a/3.0
            self.x.append((A + B) - a)
            self.x.append(-0.5 * (A + B) - a)
            self.x.append(0.5 * math.sqrt(3.0) * (A - B))
            
            if math.fabs(self.x[2]) < self.eps:
                self.x[2] = self.x[1]
                return 2
            return 1

    def csqrt(self,x,y):
        r = math.sqrt(x * x + y * y)
        a = 0
        b = 0
        
        if y == 0:
            r = math.sqrt(r)
            if x >= 0:
                a = r
                b = 0
            else:
                a = 0
                b = r
        else:
            a = math.sqrt(0.5 * (x + r))
            b = 0.5 * y/a

        return a,b

    def solve_p4_bi(b, d):
        D = b * b - 4 * d
        if D >= 0:
            sD = math.sqrt(D)
            x1 = (-b + sD) / 2.0
            x2 = (-b - sD) / 2.0

            if x2 >= 0:
                sx1 = math.sqrt(x1)
                sx2 = math.sqrt(x2)

                self.x.append(-sx1)
                self.x.append(sx1)
                self.x.append(-sx2)
                self.x.append(sx2)
                return 1

            if x1 < 0:
                sx1 = math.sqrt(-x1)
                sx2 = math.sqrt(-x2)
                self.x.append(-sx1)
                self.x.append(sx1)
                self.x.append(0)
                self.x.append(sx2)
                return 2
        else:
            sD2 = 0.5 * math.sqrt(-D)
            x0,x1 = self.csqrt(-0.5 * b, sD2)
            x2,x3 = self.csqrt(-0.5 * b, -sD2)
            self.x.append(x0)
            self.x.append(x1)
            self.x.append(x2)
            self.x.append(x3)
            return 0  

    def solve_p4_de(self, b, c, d):
        if math.fabs(c) < self.eps * (math.fabs(b) + math.fabs(d)):
            return self.solve_p4_bi(b,d)

        res3 = self.solve_p3(2*b, b*b-4*d, -c*c)
        if res3 > 1:
            if self.x[0] > self.x[1]:
               self.x[0], self.x[1] = self.x[1], self.x[0]
               
            if self.x[2] < self.x[1]:
                self.x[1], self.x[2] = self.x[2], self.x[1]
                
                if self.x[0] > self.x[1]:
                    self.x[0], self.x[1] = self.x[1], self.x[0]

            if self.x[0] > 0:
                sz1 = math.sqrt(self.x[0])
                sz2 = math.sqrt(self.x[1])
                sz3 = math.sqrt(self.x[2])

                if c > 0:
                    self.x = []
                    self.x.append((-sz1 - sz2 - sz3)/2.0)
                    self.x.append((-sz1 + sz2 + sz3)/2.0)
                    self.x.append(( sz1 - sz2 + sz3)/2.0)
                    self.x.append(( sz1 + sz2 - sz3)/2.0)
                    return 4

                self.x = []
                self.x.append((-sz1 - sz2 + sz3)/2.0)
                self.x.append((-sz1 + sz2 - sz3)/2.0)
                self.x.append(( sz1 - sz2 - sz3)/2.0)
                self.x.append(( sz1 + sz2 + sz3)/2.0)
                return 4

            sz1 = math.sqrt(-self.x[0])
            sz2 = math.sqrt(-self.x[1])
            sz3 = math.sqrt(self.x[2])

            if c > 0:
                self.x = []
                self.x.append(-sz3/2.0)
                self.x.append((-sz1 - sz2 )/2.0)
                self.x.append(sz3/2.0)
                self.x.append(( -sz1 - sz2)/2.0)
                return 0

            self.x = []
            self.x.append(sz3/2.0)
            self.x.append((-sz1 + sz2 )/2.0)
            self.x.append(-sz3/2.0)
            self.x.append(( sz1 + sz2)/2.0)
            return 0

        sz1 = math.sqrt(self.x[0])
        szr, sz11 = self.csqrt(self.x[1], self.x[2])

        if c > 0:
            self.x = []
            self.x.append(-sz1/2.0 - szr)
            self.x.append(-sz1/2.0 + szr)
            self.x.append(sz1/2.0)
            self.x.append(sz1)
            return 2

        self.x = []
        self.x.append(sz1/2.0 - szr)
        self.x.append(sz1/2.0 + szr)
        self.x.append(-sz1/2.0)
        self.x.append(sz1)
        return 2

    def n4_step(self, x,a, b, c, d):
        fxs = ((4*x+3*a)*x+2*b)*x+c
        if fxs == 0:
            return 1e99

        fx = (((x+a)*x+b)*x+c)*x+d
        return x - fx/fxs
                

    def solve_p4(self, a, b, c, d):
        d1 = d + 0.25 * a * ( 0.25 * b * a - 3./64 * a * a * a - c)
        c1 = c + 0.5 * a * (0.25 * a * a - b)
        b1 = b - 0.375 * a * a

        res = self.solve_p4_de(b1,c1,d1)
        
        if res == 4:
            self.x[0] = self.x[0] - a/4
            self.x[1] = self.x[1] - a/4
            self.x[2] = self.x[2] - a/4
            self.x[3] = self.x[3] - a/4

        else:
            if res == 2:
                self.x[0] = self.x[0] - a/4
                self.x[1] = self.x[1] - a/4
                self.x[2] = self.x[2] - a/4
            else:
                self.x[0] = self.x[0] - a/4
                self.x[2] = self.x[2] - a/4

        if res > 0:
            self.x[0] = self.n4_step(self.x[0],a,b,c,d)
            self.x[1] = self.n4_step(self.x[1],a,b,c,d)

        if res > 2:
            self.x[2] = self.n4_step(self.x[2],a,b,c,d)
            self.x[3] = self.n4_step(self.x[3],a,b,c,d)

        return res

        
            
            

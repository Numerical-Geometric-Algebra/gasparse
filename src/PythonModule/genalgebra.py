#!/usr/bin/env python

import gasparse
import jinja2
class algebra:
    def __init__(self,ga,name):
        self.name = name
        self.product_names = ["geometric","inner","outer"]
        self.signs = []
        self.bitmaps = []
        self.zerosigns = []
        for names in self.product_names:
            bitmap,sign = ga.cayley(names+"inverted")
            self.signs.append(sign)
            self.bitmaps.append(bitmap)
            self.zerosigns.append(self.compute_zero_sign(sign))
        self.metric = ga.metric()
        self.metric_size = len(self.metric)
        self.size = len(self.bitmaps[0])
        self.nproducts = len(self.product_names)
        self.grades,self.position,self.gradesize = ga.grademap()
        self.ngrades = len(self.gradesize)
        self.compute_gradesbitmap()
        self.compute_reverse()

    def GRADE(self,bitmap):
        return bin(bitmap).count('1')

    def compute_reverse(self):
        self.reverse = [1]*self.ngrades
        for i in range(self.ngrades):
            if(i&2):
                self.reverse[i] = -1

    def compute_zero_sign(self,sign):
        zerosign = [False]*len(sign)
        for j in range(len(sign)):
            for i in range(len(sign)):
                if(sign[i][j] != 0):
                    zerosign[j] = True
                    break;
        return zerosign

    def compute_gradesbitmap(self):
        lst = [0]*self.ngrades
        self.gradesbitmaplen = [0]*self.ngrades
        for i in range(self.ngrades):
            lst_sub = [0]*self.gradesize[i]
            lst[i] = lst_sub;
            self.gradesbitmaplen[i] = self.gradesize[i]

        for i in range(len(self.position)):
            lst[self.GRADE(i)][self.position[i]] = i;

        self.gradesbitmap = lst



vga = gasparse.GA(3) # 3D VGA
cga = gasparse.GA(4,1) # 3D CGA
dga = gasparse.GA(0,0,1) # Dual Numbers

ga0 = algebra(vga,"3DVGA")
ga1 = algebra(cga,"3DCGA")


typemin = 2

# algebras = [ga0,ga1]
algebras = [ga0]
nalgebras = len(algebras)

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader,trim_blocks=True)


multivectorgend = templateEnv.get_template("multivector_gen.c.src")
multivectorgend_h = templateEnv.get_template("multivector_gen.h.src")

types = ['dense','blades']
Types = ['Dense','Blades']

textmultivectorc = multivectorgend.render(algebras=algebras,nalgebras=nalgebras,types=types,Types=Types,typemin=typemin)
textmultivectorh = multivectorgend_h.render(algebras=algebras,nalgebras=nalgebras,types=types,Types=Types,typemin=typemin)

# to save the results
with open("multivector_gen.c", "w") as fh:
    fh.write(textmultivectorc)

with open("multivector_gen.h", "w") as fh:
    fh.write(textmultivectorh)

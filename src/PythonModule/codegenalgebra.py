#!/usr/bin/env python

import gasparse
import jinja2
class algebra:
    def __init__(self,bitmap,sign):
        self.bitmap = bitmap
        self.sign = sign
        self.size = len(bitmap)
        self.compute_zero_sign()
    def grade(self,grade):
        return bin(grade).count('1')

    def set_grademap(self,position,gradesize):
        self.position = position
        self.gradesize = gradesize
        self.ngrades = len(gradesize)


    def compute_zero_sign(self):
        self.zerosign = [False]*len(self.sign)
        for j in range(len(self.sign)):
            for i in range(len(self.sign)):
                if(self.sign[i][j] != 0):
                    self.zerosign[j] = True
                    break;

    def compute_gradesbitmap(self):
        lst = [0]*self.ngrades
        self.gradesbitmaplen = [0]*self.ngrades
        for i in range(self.ngrades):
            lst_sub = [0]*self.gradesize[i]
            lst[i] = lst_sub;
            self.gradesbitmaplen[i] = self.gradesize[i]

        for i in range(len(self.position)):
            lst[self.grade(i)][self.position[i]] = i;

        self.gradesbitmap = lst


vga = gasparse.GA(3) # 3D VGA
bitmap,sign = vga.cayley("innerinverted")
grade,position,gradesize = vga.grademap()

ga0 = algebra(bitmap,sign)
ga0.set_grademap(position,gradesize)
ga0.compute_gradesbitmap()

algebras = [ga0]
nalgebras = len(algebras)

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader,trim_blocks=True)


multivectorgend = templateEnv.get_template("codegend_multivector.c.src")

textmultivectorc = multivectorgend.render(algebras=algebras,nalgebras=nalgebras)

# to save the results
with open("codegend_multivector.c", "w") as fh:
    fh.write(textmultivectorc)

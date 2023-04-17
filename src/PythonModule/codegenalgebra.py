#!/usr/bin/env python

import gasparse
import jinja2
class algebra:
    def __init__(self,bitmap,sign):
        self.bitmap = bitmap
        self.sign = sign
        self.size = len(bitmap)
    def grade(self,grade):
        return bin(grade).count('1')

    def set_grade_size(self,grade_size,max_grade):
        self.grade_size = grade_size
        self.max_grade = max_grade


vga = gasparse.GA(3) # 3D VGA
bitmap,sign = vga.cayley("inverted")

algebras = [algebra(bitmap,sign)]
nalgebras = len(algebras)

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader,trim_blocks=True)


multivectorgend = templateEnv.get_template("codegend_multivector.c.src")

textmultivectorc = multivectorgend.render(algebras=algebras,nalgebras=nalgebras)

# to save the results
with open("codegend_multivector.c", "w") as fh:
    fh.write(textmultivectorc)

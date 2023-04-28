#!/usr/bin/env python

import gasparse
from gasparse import multivector


ga = gasparse.GA(3,compute_mode="generic")
# dtypes = ['blades','sparse','dense']
dtypes = ['blades','sparse','dense']

alist = []
blist = []

for dtype in dtypes:
    ga.default(dtype)
    blades = ga.blades()
    locals().update(blades)
    a = e1+e2+e3
    b = e3+e123
    c = e12
    alist.append(a)
    blist.append(b)
    print("a&b=",a&b)
    print("a^b=",a^b)
    print("a*b=",a*b)
    print("a|b=",a|b)
    print("a+b=",a+b)
    print("a-b=",a-b)
    print("a.dual()=",a.dual())
    print("a.undual()=",a.undual())
    print("b([0,1])=",b([0,1]))
    print("~b=",~b)
    print("add",multivector.add(e1,e2,e3,e123))
    print("product",multivector.geometric_product(e1,a,b,c))
    print("product",multivector.geometric_product(a,b,c))
    print("scalar add",1.234+a)
    print("multiply by scalar",0.123*a)
    print()

# print(alist[0]&blist[1])
# print(alist[1]&blist[0])
# ga.default("sparse",False)
# blades = ga.blades()
# locals().update(blades)
# cs = e1+e2+e3
# ds = e1+e2+e13


'''
for dtype in dtypes:
    print(dtype  + ':')
    ga.default(dtype,True)

    blades = ga.blades(grades=1)
    locals().update(blades)
    I = e1*e2*e3

    c = e1+e2+e3 + (e1*e2) + (e2*e3*e1)
    d = (3*e1) + (0.123*e2) + (e1*e3)

    print()
    print(c&d)
    print(c.dual()|d)
    print((c.dual()|d) - (c&d))
'''
'''
print('mixed:')
e1 = ga.multivector([1.0],['e1'],dtype='blades')
e2 = ga.multivector([1.0],['e2'],dtype='dense')
e3 = ga.multivector([1.0],['e3'],dtype='sparse')
# e4 = ga.multivector([1.0],['e4'],dtype='blades')
print(e1)

a = e1 + e2
b = (e1*e2) + e3
print(a*(e1+e3))
print(a*(e1-e3))
print(~a*(e1-e3))
print(a|b)
print(a^b)
print(1.23343*a)
print()
'''

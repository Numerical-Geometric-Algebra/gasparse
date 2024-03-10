#!/usr/bin/env python3

import gasparse
from gasparse import multivector

ga = gasparse.GA(3,print_type_mv=1,compute_mode="generated")

alist = []
blist = []
dtypes = ['blades','dense']
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


print(alist[0]&blist[1])
print(alist[1]&blist[0])

# ga.default("blades",True)
# blades = ga.blades()
'''
dtypes = ['blades','dense']
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

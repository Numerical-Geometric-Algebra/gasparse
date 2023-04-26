#!/usr/bin/env python

import gasparse
from gasparse import multivector


ga = gasparse.GA(3,print_type_mv=0,large=True)
# dtypes = ['blades','sparse','dense']

ga.default("dense",False)
blades = ga.blades(grades=1)
locals().update(blades)
print(e1&e2);

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

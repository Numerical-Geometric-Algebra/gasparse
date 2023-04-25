#!/usr/bin/env python

import gasparse
from gasparse import multivector

ga = gasparse.GA(3,print_type_mv=0,large=True)
dtypes = ['blades','sparse','dense']
# dtypes = ['blades']

for dtype in dtypes:
    print(dtype  + ':')
    e1 = ga.multivector([1.0],['e1'],dtype=dtype)
    e2 = ga.multivector([1.0],['e2'],dtype=dtype)
    e3 = ga.multivector([1.0],['e3'],dtype=dtype)
    I = e1*e2*e3

    c = e1+e2+e3 + (e1*e2) + (e2*e3*e1)

    print(c.dual() - c*~I)
    print()
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

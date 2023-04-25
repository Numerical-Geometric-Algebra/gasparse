#!/usr/bin/env python3

import gasparsegen as gasparse
from gasparsegen import multivector


ga = gasparse.GA(3,print_type_mv=0)

dtypes = ['blades','dense']

for dtype in dtypes:
    print(dtype  + ':')
    e1 = ga.multivector([1.0],['e1'],dtype=dtype)
    e2 = ga.multivector([1.0],['e2'],dtype=dtype)
    e3 = ga.multivector([1.0],['e3'],dtype=dtype)
    print(e1)
    I = e1*e2*e3

    c = e1+e2+e3 + (e1*e2) + (e2*e3*e1)

    print(c.dual() - c*~I)
    print()


    a = e1 + e2
    b = (e1*e2) + e3
    print(a*(e1+e3))
    print(a*(e1-e3))
    print(~a*(e1-e3))
    print(a|b)
    print(a^b)
    print(1.23343*a)
    print()

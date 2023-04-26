#!/usr/bin/env python3

import gasparsegen as gasparse
from gasparsegen import multivector


ga = gasparse.GA(3,print_type_mv=0)

# ga.default("blades",True)
# blades = ga.blades()

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

#!/usr/bin/env python

import gasparse
from gasparse import multivector


ga = gasparse.GA(p=3,print_type_mv=1)
# a = ga.multivector([2,1,6,8,9,9,4],['','e1','e12','e123','e1234','e12345','e123456'])
# b = ga.multivector([7.654,54.233,76.888],['e4','e34','e134'])

e1 = ga.multivector([1.0],['e1']) # e1 = ga.blade('e1')
e2 = ga.multivector([1.0],['e2'])
e3 = ga.multivector([1.0],['e3'])

# e3 = ga.multivector([1.0],['e3'])
# e4 = ga.multivector([1.0],['e4'])
# e5 = ga.multivector([1.0],['e5'])
# e6 = ga.multivector([1.0],['e6'])
# a = ga.multivector([1.0,1.0],['e1','e2'],dtype='blades')
b = ga.multivector([1.0,1.0],['e1','e123'],dtype='blades')
e12 = multivector.geometric_product(e1,e2)

print(e12)
print(b)

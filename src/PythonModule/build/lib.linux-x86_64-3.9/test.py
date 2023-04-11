#!/usr/bin/env python

import gasparse
from gasparse import multivector

ga = gasparse.GA(p=3,r=1,print_type_mv=1)

e1 = ga.multivector([1.0],['e1'],dtype='blades')
e2 = ga.multivector([1.0],['e2'],dtype='blades')
e3 = ga.multivector([1.0],['e3'],dtype='blades')
e4 = ga.multivector([1.0],['e4'],dtype='blades')

b = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='blades')
print(b)
# b_sparse = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='sparse')
# b_dense = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='dense')

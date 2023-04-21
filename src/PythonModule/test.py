#!/usr/bin/env python

import gasparse
from gasparse import multivector

ga = gasparse.GA(3,print_type_mv=1)

e1 = ga.multivector([1.0],['e1'],dtype='sparse')
e2 = ga.multivector([1.0],['e2'],dtype='dense')
e3 = ga.multivector([1.0],['e3'],dtype='blades')
# e4 = ga.multivector([1.0],['e4'],dtype='sparse')
print(e1)

a = e1 + e2
b = (e1*e2) + e3
print(a*(e1+e3))
print(a*(e1-e3))
print(~a*(e1-e3))
print(a|b)
print(a^b)
# bitmapinv,signinv = ga.cayley("inverted")
# bitmap,sign = ga.cayley()
# bitmapouter,signouter = ga.cayley("outer")

# b = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'], dtype='sparse')
# c = multivector.geometric_product(e1,b)
# b_sparse = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='sparse')
# b_dense = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='sparse')

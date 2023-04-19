#!/usr/bin/env python

import gasparse
from gasparse import multivector

ga = gasparse.GA(4,1)

e1 = ga.multivector([1.0],['e1'],dtype='sparse')
e2 = ga.multivector([1.0],['e2'],dtype='sparse')
e3 = ga.multivector([1.0],['e3'],dtype='sparse')
# e4 = ga.multivector([1.0],['e4'],dtype='sparse')
print(e1)
# bitmapinv,signinv = ga.cayley("inverted")
# bitmap,sign = ga.cayley()
# bitmapouter,signouter = ga.cayley("outer")

# b = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'], dtype='sparse')
# c = multivector.geometric_product(e1,b)
# b_sparse = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='sparse')
# b_dense = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='sparse')

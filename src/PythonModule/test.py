#!/usr/bin/env python

import gasparse
from gasparse import multivector

# ga = gasparse.GA(3,print_type_mv=1)

ga = gasparse.GA(3,print_type_mv=1)
e1 = ga.multivector([1.0],['e1'],dtype='blades')
print(e1)

# ga = gasparse.GA(3)


# dtypes = ['blades','sparse','dense']
'''
dtypes = ['blades','dense']

for dtype in dtypes:
    print(dtype  + ':')
    e1 = ga.multivector([1.0],['e1'],dtype=dtype)
    e2 = ga.multivector([1.0],['e2'],dtype=dtype)
    e3 = ga.multivector([1.0],['e3'],dtype=dtype)
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
# bitmapinv,signinv = ga.cayley("inverted")
# bitmap,sign = ga.cayley()
# bitmapouter,signouter = ga.cayley("outer")

# b = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'], dtype='blades')
# c = multivector.geometric_product(e1,b)
# b_sparse = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='blades')
# b_dense = ga.multivector([1.0,1.0,1.0],['e1','e34','e123'],dtype='blades')
'''

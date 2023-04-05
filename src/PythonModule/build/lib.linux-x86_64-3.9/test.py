#!/usr/bin/env python

import gasparse

ga = gasparse.GA(p=3)
ga.add_basis(q=1)
a = ga.multivector([1.233,6.454,8.989],['e1','e12','e123'])
b = ga.multivector([7.654,54.233,76.888],['e4','e34','e134'])
c = a*b

#!/usr/bin/env python3

import gasparse
from gasparse import multivector

ga = gasparse.GA(4)
ga.default('sparse')
blades = ga.blades()
locals().update(blades)
a = ga.multivector([1,1],["e12","e13"])
b = ga.multivector([1,1,1],["e12","e34","e13"])
print(multivector.exp(a))
print(multivector.exp(b))

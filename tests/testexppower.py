#!/usr/bin/env python3

import gasparse
from gasparse import multivector

import math
ga = gasparse.GA(8,compute_mode='generic')
ga.default("sparse",1e-32)
blades = ga.blades()
locals().update(blades)

b = e12+e13+e34+e15+e36+e25+e16+e27+e67+e18+e38+e58
d = multivector.exp(b)

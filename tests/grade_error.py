import gasparse
import numpy as np
n = 10
ga = gasparse.GA(3)
arr = np.random.rand(n,ga.size(1,2)) # innermost dimension must be the the size of grades 1 and 2
x = ga.multivector(arr.tolist(),grades=[12]) # only accepts lists as input


import timeit
import gasparse
import numpy as np

# Stack smashing happened when the string of the printing of the multivectors got too big.
# Rule: as of now do not print big multivector arrays.
# Numpy Compresses the output when considering a very big array!!!
# Should Follow numpy's advice
ga = gasparse.GA(3)

size = 1000000
mvsize = ga.size()
ones = 0.5*np.ones([size,mvsize])
x_array = np.random.rand(size,mvsize) - ones
x =  ga.multivector(x_array.tolist())

y = x*x # this will not give error
print(x) # but this will give segmentation fault
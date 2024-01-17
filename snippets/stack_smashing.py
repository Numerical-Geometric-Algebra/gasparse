from geo_algebra import *
import timeit

# Stack smashing happened when the string of the printing of the multivectors got too big.
# Rule: as of now do not print big multivector arrays.
# Numpy Compresses the output when considering a very big array!!!
# Should Follow numpy's advice
x = rdn_kvector_array(vga,[1,2],100000)
y = x*x
#print(x)
import gasparse
import numpy as np
from gasparse import mvarray as mv

ga_large = gasparse.GA(3,compute_mode='large')
ga = gasparse.GA(3,compute_mode='generic')

arr1 = np.random.rand(8).tolist()
arr2 = np.random.rand(8).tolist()
arr3 = np.random.rand(8).tolist()

x_large = ga_large.mvarray(arr1,dtype='sparse')
y_large = ga_large.mvarray(arr2,dtype='sparse')
z_large = ga_large.mvarray(arr3,dtype='sparse')

x = ga.mvarray(arr1,dtype='sparse')
y = ga.mvarray(arr2,dtype='sparse')
z = ga.mvarray(arr3,dtype='sparse')

print(x*y - x_large*y_large)
# print(x_large*y_large)
print()

print((x^y) - (x_large^y_large))
# print(x_large^y_large)
print()

print((x|y) - (x_large|y_large))
# print(x_large|y_large)
print()

print((x&y) - (x_large&y_large))
# print(x_large&y_large)
print()

print((x+y) - (x_large+y_large))

print(((x+y) - (x_large+y_large)).type())

locals().update(ga_large.basis())

print("checking auto discard mechanism:")
w = ga_large.mvarray([1,1,3e-33,4],grades=[0,2],dtype='sparse')

# w = 1e-31 + 1*e12 + 3*e13 + 4*e23
print(w)
print(w+e1)
print()

print("geometric product")
print(z*x*y)
print(z_large*x_large*y_large)
print(mv.tprod(z_large,x_large,y_large,ptype='geometric'))
print(mv.tprod(z,x,y,ptype='geometric'))

print("outer product")
print(z^x^y)
print(z_large^x_large^y_large)
print(mv.tprod(z_large,x_large,y_large,ptype='outer'))
print(mv.tprod(z,x,y,ptype='outer'))

print("inner product")
print((z|x)|y)
print((z_large|x_large)|y_large)
print(mv.tprod(z_large,x_large,y_large,ptype='inner'))
print(mv.tprod(z,x,y,ptype='inner'))

print("regressive product")
print((z&x)&y)
print((z_large&x_large)&y_large)
print(mv.tprod(z,x,y,ptype='regressive'))

# This will give error because the ternary regressive product is not yet implemented for large multivectors
print(mv.tprod(z_large,x_large,y_large,ptype='regressive'))



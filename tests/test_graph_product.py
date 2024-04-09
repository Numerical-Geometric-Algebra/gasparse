import gasparse
import numpy as np

ga_large = gasparse.GA(3,compute_mode='large')
ga = gasparse.GA(3,compute_mode='generic')

arr1 = np.random.rand(8).tolist()
arr2 = np.random.rand(8).tolist()

x_large = ga_large.mvarray(arr1,dtype='sparse')
y_large = ga_large.mvarray(arr2,dtype='sparse')

x = ga.mvarray(arr1,dtype='sparse')
y = ga.mvarray(arr2,dtype='sparse')

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

z = ga_large.mvarray([1,1,3e-33,4],grades=[0,2],dtype='sparse')

# z = 1e-31 + 1*e12 + 3*e13 + 4*e23
print(z)
print(z+e1)
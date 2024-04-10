import gasparse
import numpy as np
from gasparse import mvarray as mv

ga = gasparse.GA(3,compute_mode='large')
arr = np.random.rand(3,8).tolist()
x = ga.mvarray(arr,dtype='sparse')

print("product:")
print(x[0])
print(x[0]*x[1])
print(x[0]*x[1]*x[2])
print(x.prod())
print("sum:")
print(x[0]+x[1]+x[2])
print(x.sum())
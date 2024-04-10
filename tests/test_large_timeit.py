import gasparse
import timeit
import numpy as np

gal = gasparse.GA(12,compute_mode='large')
ga = gasparse.GA(12,compute_mode='generic')

arr1 = np.random.rand(10,ga.size(1)).tolist()
arr2 = np.random.rand(10,ga.size(1)).tolist()

xl1 = gal.mvarray(arr1,grades=1,dtype='sparse')
x1 = ga.mvarray(arr1,grades=1,dtype='sparse')

xl2 = gal.mvarray(arr2,grades=1,dtype='sparse')
x2 = ga.mvarray(arr2,grades=1,dtype='sparse')

time_generic = timeit.timeit(lambda: x1*x2, number=5)/5
time_large = timeit.timeit(lambda: xl1*xl2, number=5)/5

print("generic is ", time_generic/time_large, " times slower then large",sep='')

time_generic = timeit.timeit(lambda: x1+x2, number=5)/5
time_large = timeit.timeit(lambda: xl1+xl2, number=5)/5

print("generic is ", time_generic/time_large, " times slower then large",sep='')

time_generic = timeit.timeit(lambda: x1.prod(), number=5)/5
time_large = timeit.timeit(lambda: xl1.prod(), number=5)/5

print("generic is ", time_generic/time_large, " times slower then large",sep='')

time_generic = timeit.timeit(lambda: x1.sum(), number=5)/5
time_large = timeit.timeit(lambda: xl1.sum(), number=5)/5

print("generic is ", time_generic/time_large, " times slower then large",sep='')
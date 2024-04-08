import gasparse
import numpy as np

ga = gasparse.GA(3)
arr = np.random.rand(3,2,ga.size(1,2))
x = ga.mvarray(arr.tolist(),grades=[1,2])
print(x)

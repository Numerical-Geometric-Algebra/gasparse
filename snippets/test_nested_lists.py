import gasparse

ga = gasparse.GA(3)
x = ga.multivector_array([[[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]]],grades=1)
#x = ga.multivector_array([1,2,3],grades=1)
print(x)
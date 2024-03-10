import gasparse

ga = gasparse.GA(4,1)
x = ga.multivector([-1,1,2,3,4,5],grades=[0,1])
print(type(x(0)))
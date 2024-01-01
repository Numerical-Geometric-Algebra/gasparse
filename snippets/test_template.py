import gasparse
from gasparse import multivector

ga = gasparse.GA(3,compute_mode="generic", print_type_mv=1)

basis = ga.basis()
locals().update(basis)

x = ga.multivector([1,23,4,4,8],[1.0,e1,e2,e12,e13])
y = ga.multivector([1,23,4,4],['e','e1','e2','e12'])
z = ga.multivector([1,2,3,4,5,6,7,8])
u = ga.multivector([-1,1,2,3],grades=[0,2])
v = ga.multivector([1,2,3],grades=2)

print(x.list())
print(x.list(bitmap=True))
print(x.list(1))
print(x.list([1,2]))
print(x.list(1,bitmap=True))

y = x.cast("dense")
#values,basis = x.list(1)



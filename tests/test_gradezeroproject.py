import gasparse
from gasparse import multivector as mv
''' Grade projection to zero creates a new type of multivector: A scalar type multivector. Operations with this elements are optimized.
'''


ga = gasparsegen.GA(3)
x = ga.multivector([[2,2,3,4],[5,6,7,8]],grades=[0,1])
y = ga.multivector([[2,2,3,4],[6,6,7,8]],grades=[0,1])
print(x(0))
print(x(0).Type())

print(x(0) - y(0))
print()
print(y(0) - x(0))
print()
print(x*y(0))
print()
print(x*y(0))
print()
print(y(0)-x)

print(x/y(0))
print()
print(y(0)/10)
print()
print(10/y(0))
print()
print(x/10)

x = ga.multivector([[1],[2]],basis=['e'],dtype='scalar')
# print(10/x) # This is not valid since x is not a scalar


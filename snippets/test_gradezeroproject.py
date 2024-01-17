import gasparsegen

ga = gasparsegen.GA(3,compute_mode="generated")
x = ga.multivector([[2,2,3,4],[5,6,7,8]],grades=[0,1])
y = ga.multivector([[2,2,3,4],[6,6,7,8]],grades=[0,1])
# print(x(0))
# print(x(0).Type())

# print(x(0) - y(0))
# print()
# print(y(0) - x(0))
# print()
# print(x*y(0))
# print()
print(x*y(0))
# print()
# print(y(0)-x)
import gasparse

ga = gasparse.GA(3)
basis = ga.basis()
locals().update(basis)
y = ga.multivector_array([[[[1.1,2.2,3.3],[4,5,6]],[[7,8,9],[10,11,12]]],
                          [[[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]]]],basis=[e1,e2,e3])
x = ga.multivector_array([1,2,3],grades=1)
print(x)
print()
print(y)
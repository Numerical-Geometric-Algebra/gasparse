import gasparse

ga = gasparse.GA(3)
basis = ga.basis()
locals().update(basis)

x = ga.multivector([[[1999,1,2,3,4,5,6,7],[0.1,1.1,2.2,3.3,4.4,5.5,6.6,7.7]],[[10,1.2,2.4,3.6,4.7,5.8,6.9,7.4],[0.11111,1.1,2.2,3.3,4.4,5.5,6.6,7.7]]])
# x = ga.multivector([[10,1,2,3,4,5,6,7],[0.1,1.1,2.2,3.3,4.4,5.5,6.6,7.7]])
# print(x[0])
# print(x[1])
# print(x[0])
print(x[:,1])
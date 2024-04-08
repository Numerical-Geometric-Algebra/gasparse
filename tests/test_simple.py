import gasparse
import numpy as np

ga = gasparse.GA(3,print_type_mv=1,print_type=1)
ga.set_precision(1e-30)

x = ga.multivector([1],basis=['e1'])
y = ga.multivector([1,1e-31],basis=['e2','e3'])

print(x*y)
print(x+y)
print((29*y) + x)

a = ga.multivector([[1e-32,0,0],[0,2,0],[0,0,3]],['e1','e2','e3'])
b = ga.multivector([[1,1e-31],[5,2]],['e2','e3'])

# print(a)
# print()
# print(b)
# print()
# print(a.sum())
# print(a.prod())
# print(b.sum())
# print(b.prod())

u = ga.multivector([[[[1.234567891234567891e10],[2],[5]],[[3],[4],[6]]],[[[0.1],[0.2],[0.5]],[[0.3],[0.4],[0.6]]]],['e1'])
# u = ga.multivector([[[1],[2],[5]],[[3],[4],[6]]],['e1'])
print()
print(u)
# print()
# print(u[0][0])
# print()
# print(u[0])


# test for sparse multivector type

# product by scalars
# addition of two multivectors
# product of two multivectors

# Which of the two is more efficient:
    # ternary products
    # ternary products nested


# binary_sparse_add
# binary_sparse_product
# ternary_sparse_product_nested
# ternary_sparse_product
# atomic_sparse_add
# atomic_sparse_product
# binary_mixed_product
# atomic_mixed_add
# atomic_mixed_product
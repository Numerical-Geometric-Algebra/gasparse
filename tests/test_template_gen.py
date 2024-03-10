import gasparse
from gasparse import multivector

vga = gasparse.GA(3,compute_mode="generated", print_type_mv=1)
cga = gasparse.GA(4,1,compute_mode="generic", print_type_mv=1)


basis = vga.basis()
locals().update(basis)

x = vga.multivector([1,23,4,4,8],[1.0,e1,e2,e12,e13])
y = vga.multivector([1,23,4,4],['e','e1','e2','e12'],dtype='blades')
z = vga.multivector([1,2,3,4,5,6,7,8])
l = vga.multivector([1,2,3],grades=1)
s = cga.multivector([1,2,3,4,5],grades=1)


'''
print(x.list())
print(x.list(bitmap=True))
print(x.list(1))
print(x.list([1,2]))
print(x.list(1,bitmap=True))
values,basis = x.list(1)
print(ga.size(1))
'''

t = x*s
r = multivector.geometric_product(x,s,l)



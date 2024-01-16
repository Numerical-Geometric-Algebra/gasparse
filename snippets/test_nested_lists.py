import gasparsegen

ga = gasparsegen.GA(3,compute_mode="generic")
cga = gasparsegen.GA(4,1,compute_mode="generic")
basis = ga.basis()
locals().update(basis)
xlst = [[[[1.1],[4.0]],
         [[7],[10]]],
        [[[1],[4]],
         [[7],[10]]],
        [[[1],[4]],
         [[7],[10]]],
        [[[1],[4]],
         [[7],[10.9999]]]]

ylst = [[[[0.0],[6.0]],
        [[8.0],[12.0]]],
       [[[2.0],[6.0]],
        [[9.0],[11.4]]],
       [[[1.0],[6.4532]],
        [[8.0],[10.123]]],
       [[[1.879],[0.0]],
        [[7.0],[10.89788]]]]

xbasis = ['e1']
ybasis = ['e23']
wbasis = ['e2']


x = ga.multivector(xlst,basis=xbasis)
y = ga.multivector(ylst,basis=ybasis)
w = ga.multivector([1.0],basis=wbasis)
z = w + y
z = w*y
z = x + y
z = x*y
z = 1*y
z = 1+y
print((x+y)(1))

y = cga.multivector(ylst,basis=ybasis)
w = cga.multivector([1.0],basis=wbasis)
z = x + y
z = x*y
z = 1*y
z = 1+y
z = w + y
z = w*y
print((x+y)(1))

vgagen = gasparsegen.GA(3,compute_mode="generated",print_type_mv=0)
cgagen = gasparsegen.GA(4,1,compute_mode="generated",print_type_mv=0)

y = vgagen.multivector(ylst,basis=ybasis)
x = vgagen.multivector(xlst,basis=xbasis)
z = x + y
z = x*y
z = 1*y
z = 1+y
z = e1 + y
z = e1*y
print((x+y)(1))
y = cgagen.multivector(ylst,basis=ybasis)
z = x + y
z = x*y
z = 1*y
z = 1+y
z = e1 + y
z = e1*y
print((x+y)(1))

vgagen = gasparsegen.GA(3,compute_mode="generated",print_type_mv=0)
x = vgagen.multivector(xlst,basis=xbasis)
y = cga.multivector(ylst,basis=ybasis)
z = x + y
z = x*y
z = 1*y
z = 1+y
z = e1 + y
z = e1*y
print((x+y)(1))

vgalarge = gasparsegen.GA(3,compute_mode="large",print_type_mv=0)
cgalarge = gasparsegen.GA(4,1,compute_mode="large",print_type_mv=0)

y = vgalarge.multivector(ylst,basis=ybasis)
x = vgalarge.multivector(xlst,basis=xbasis)
z = x + y
z = x*y
z = 1*y
z = 1+y
z = e1 + y
z = e1*y
print((x+y)(1))


#print(x.list())
#print(x)
#print()
#print(y)
#print(y.list()[0])
#print(y.list(1)[0])
#print()
#print(y.list(0)[0])



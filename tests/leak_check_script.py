#!/usr/bin/env python

import gasparse
from gasparse import multivector

cmodes = ["generated","large","generic"]
dtypes = ['blades','sparse','dense']
for cmode in cmodes:

    ga = gasparsegen.GA(3,compute_mode=cmode,print_type_mv=0)
    if cmode == 'generic' or cmode == 'large':
        dtypes = ['blades','sparse','dense']
    else:
        dtypes = ['blades','dense']

    alist = []
    blist = []

    for dtype in dtypes:
        print(dtype,cmode)
        x = ga.size()
        x = ga.size(1)
        ga.default(dtype)
        basis = ga.basis()
        locals().update(basis)
        a = e1+e2+e3
        b = 1.234+e3+e123
        c = e12
        alist.append(a)
        blist.append(b)
        x = ga.multivector([1,2,3],grades=1)
        x = a&b
        x = a^b
        x = a*b
        x = a|b
        x = a+b
        x = a-b
        x = a.dual()
        x = a.undual()
        x = b([0,1])
        x = ~b
        x = b(0)
        x = multivector.add(e1,e2,e3,e123)
        x = multivector.geometric_product(e1,a,b,c)
        x = multivector.geometric_product(a,b,c)
        x = 1.234+a
        x = 0.123*a
        x = b.list()
        x = b.list(1)
        x = b.list([1,2])
        x = a.grade()
        x = a.cast("blades")
    x = alist[0]*blist[1]
    if len(alist) >= 3:
        x = alist[0]*blist[2]
        x = alist[1]*blist[2]

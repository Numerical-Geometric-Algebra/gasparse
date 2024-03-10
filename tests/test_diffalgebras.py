#!/usr/bin/env python3
import gasparse
from gasparse import multivector
p = 3
# modes = ['generic','large','generated']
modes = ['generated']
# dtypes = ['sparse','dense','blades']
dtypes = ['blades']


blades_list = []
for mode in modes:
    for dtype in dtypes:
        if(mode == 'generated' and dtype == 'sparse'):
            continue
        ga = gasparsegen.GA(p,compute_mode=mode)
        ga.default(dtype)
        blades = ga.blades()
        blades_list.append(blades)
        locals().update(blades)
        print((e1+e12+e23+e123)([1,3]))


for blades0 in blades_list:
    for blades1 in blades_list:
        print()
        print(blades0["e1"],blades1["e2"])
        print(blades0["e1"]*blades1["e2"])
        print(blades0["e1"]&blades1["e2"])
        print(blades0["e1"]|blades1["e12"])
        print(blades0["e1"]+blades1["e2"])
        print(blades0["e1"]-blades1["e2"])
        print(blades0["e1"].dual())
        print(blades0["e1"].undual())
        print(blades0["e1"]|1)
        print(~blades0["e12"])
        print((blades0["e12"]+blades1["e123"])(3))
        print(1+blades0["e12"])
        print(2.345*blades0["e12"])
        print(multivector.add(blades0["e1"],blades1["e2"]))
        print(multivector.add(blades0["e1"],blades1["e2"],blades0["e3"]))
        print(multivector.geometric_product(blades0["e1"],blades1["e2"],blades0["e3"]))
        print(multivector.geometric_product(blades0["e1"],blades1["e2"]))
        print(multivector.geometric_product(blades0["e1"],blades1["e2"],blades0["e3"],blades1["e12"]))


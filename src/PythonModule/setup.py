#!/usr/bin/env python
from distutils.core import setup, Extension
import sys

genalgebras = False
if "--genalgebras" in sys.argv:
    genalgebras = True
    sys.argv.remove("--genalgebras")


if genalgebras:
    macros = [("INCLUDE_GENCODE",None)] # include generated code
    setup(name="gasparsegen", version="1.0",
          ext_modules=[Extension("gasparsegen",["gasparse.c","multivector.c","largemultivector.c","multivector_gen.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"],define_macros=macros)])
else:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse",["gasparse.c","multivector.c","largemultivector.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"])])

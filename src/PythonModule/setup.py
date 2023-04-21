#!/usr/bin/env python
from distutils.core import setup, Extension
import sys

genalgebras = False
if "--genalgebras" in sys.argv:
    genalgebras == True
    sys.argv.remove("--genalgebras")

if genalgebras:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse",["gasparse.c","multivector.c","multivector_gen.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt","-DINCLUDE_GENCODE"])])
else:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse",["gasparse.c","multivector.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"])])

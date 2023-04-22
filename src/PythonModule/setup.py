#!/usr/bin/env python
from distutils.core import setup, Extension
import sys

genalgebras = False
if "--genalgebras" in sys.argv:
    genalgebras = True
    sys.argv.remove("--genalgebras")

macros = [("INCLUDE_GENCODE",None)] # include generated code


if genalgebras:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse_gen",["gasparse.c","multivector.c","multivector_gen.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"],define_macros=macros)])
else:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse",["gasparse.c","multivector.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"])])

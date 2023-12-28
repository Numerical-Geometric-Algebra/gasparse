#!/usr/bin/env python
#from distutils.core import setup, Extension
from setuptools import Extension, setup
#from setuptools.command.build_ext import build_ext as _build_ext
import sys

genalgebras = False
if "--genalgebras" in sys.argv:
    genalgebras = True
    sys.argv.remove("--genalgebras")


if genalgebras:
    macros = [("INCLUDE_GENCODE",None)] # include generated code
    setup(name="gasparsegen", version="1.0",
          ext_modules=[Extension("gasparsegen",["src/gasparse.c","src/multivector.c","src/largemultivector.c","src/multivector_gen.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt"],define_macros=macros)])
else:
    setup(name="gasparse", version="1.0",
          ext_modules=[Extension("gasparse",["src/gasparse.c","src/multivector.c","src/largemultivector.c"],\
                                 extra_compile_args = ["-O0", "-mpopcnt","-g3"])])


# extra_compile_args = ["-O0", "-mpopcnt","-DCMAKE_EXPORT_COMPILE_COMMANDS=1"])])
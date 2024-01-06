#!/usr/bin/env python
#from distutils.core import setup, Extension
from setuptools import Extension, setup
#from setuptools.command.build_ext import build_ext as _build_ext
import sys
name = "gasparse"

macros = []
if "--genalgebras" in sys.argv:
    name = "gasparsegen"
    sys.argv.remove("--genalgebras")
    macros += [("INCLUDE_GENCODE",None)]
if "--numpy" in sys.argv:
    numpy = True
    macros += [("INCLUDE_NUMPY",None)]
    sys.argv.remove("--numpy")


setup(name=name, version="1.0",
    ext_modules=[Extension(name,["src/gasparse.c","src/multivector.c","src/largemultivector.c","src/multivector_gen.c","src/multivector_array.c"],\
    extra_compile_args = ["-O0", "-mpopcnt","-g3"],define_macros=macros)])

#!/usr/bin/env python

from distutils.core import setup, Extension
setup(name="gasparse", version="1.0",
      ext_modules=[Extension("gasparse",["gasparse.c","multivector.c"],extra_compile_args = ["-O0"])])

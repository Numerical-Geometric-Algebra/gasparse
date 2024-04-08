#!/usr/bin/env python
from setuptools import Extension, setup, find_packages
import sys
''' Used for installing the package with pip. For development purposes use the other script setup_dev.py
'''
name = "gasparse"
macros = [("RELEASE_BUILD",None)]
setup(  name=name, 
        version="0.0.0a1", # Should have used git tag for the version
        packages=find_packages(exclude=("tests",)),
        ext_modules=[Extension(name,["src/gasparse.c","src/multivector_object.c","src/multivector_types.c","src/multivector_large.c","src/multivector_gen.c","src/common.c"],\
        extra_compile_args = ["-O0", "-mpopcnt","-g3"],define_macros=macros)],
        description="A python library written entirely in C for Geometric Algebras to deal with sparse multivector arrays",
        long_description = open('README.md').read(),
        license='MIT',
        author="Francisco Mendia",
        author_email="francisco.mendia.99@gmail.com",
        python_requires=">=3.7"
)

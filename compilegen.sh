#!/usr/bin/env sh
python3 setup_dev.py build
cp build/lib.linux-x86_64-cpython-311/gasparse_dev.cpython-311-x86_64-linux-gnu.so .
python3 genalgebra.py
rm build/lib.linux-x86_64-cpython-311/gasparse_dev.cpython-311-x86_64-linux-gnu.so
python3 -m build
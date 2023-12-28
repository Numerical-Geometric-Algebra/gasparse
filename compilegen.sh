#!/usr/bin/env sh

python3 setup.py build
python3 genalgebra.py
python3 setup.py --genalgebras build

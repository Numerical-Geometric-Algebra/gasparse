python3 -m build
rm -r build
cibuildwheel --platform linux
mv dist/gasparse-0.0.5a0.tar.gz wheelhouse/gasparse-0.0.5a0.tar.gz
python3 -m twine upload wheelhouse/*
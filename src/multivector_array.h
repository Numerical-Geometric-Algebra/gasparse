#ifndef MULTIVECTOR_ARRAY_H_
#define MULTIVECTOR_ARRAY_H_
#include <Python.h>
#include "gasparse.h"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gasparse_ARRAY_API
//#include <numpy/arrayobject.h>

typedef struct PyMultivectorArrayObject{
    PyObject_HEAD
    void *data; // PyMultivectorObject or any other subtype such as dense and sparse
    PyMultivectorMixedMath_Funcs *mixed;
    PyAlgebraObject *GA;
    PyMultivectorSubType *type;
    Py_ssize_t ns; // Size of shapes and strides
    Py_ssize_t *strides; // Has ns + 1 elements
    Py_ssize_t *shapes; // Has ns elements
}PyMultivectorArrayObject;



#endif // MULTIVECTOR_ARRAY_H_
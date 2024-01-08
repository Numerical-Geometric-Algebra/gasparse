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
    Py_ssize_t *strides; // Has ns + 1 elements, the 1st element is the total nbr of mvs
    Py_ssize_t *shapes; // Has ns elements
}PyMultivectorArrayObject;

typedef struct PyMvArrayIter PyMvArrayIter;

typedef int (*mvarrayiternextfunc)(PyMvArrayIter*,Py_ssize_t);

typedef struct PyMvArrayIter{
    void *data;
    Py_ssize_t basic_size;
    Py_ssize_t ns; // Size of shapes and strides
    Py_ssize_t *strides; // Has ns + 1 elements
    Py_ssize_t *shapes; // Has ns elements
    Py_ssize_t *index; // Has ns elements, size of the shape
    mvarrayiternextfunc next;
}PyMvArrayIter;

typedef struct PyMultipleArrayIter{
    PyMvArrayIter *array_iter;
    Py_ssize_t **repeat; // The number of repeated symbols, also the size of dim[i][z]
    Py_ssize_t ***dims; // A dim array for each dimension
    Py_ssize_t ns; // The number of symbols
    Py_ssize_t nm; // The number of multivector_arrays
    Py_ssize_t *index; // An index for symbol
    Py_ssize_t *shapes; // A shape for symbol
    Py_ssize_t dim; // The last iterated dimension
}PyMultipleArrayIter;


typedef PyMvArrayIter PyMultivectorArrayIter;
typedef PyMultivectorArrayObject PyMvArrayObj;
typedef PyMultivectorArrayObject PyMvArrayObject;

extern PyTypeObject PyMultivectorArrayType;
PyObject *algebra_multivector_array(PyAlgebraObject *self, PyObject *args, PyObject *kwds);

#endif // MULTIVECTOR_ARRAY_H_
#ifndef MULTIVECTOR_OBJECT_H_
#define MULTIVECTOR_OBJECT_H_
// #include "pytypedefs.h"
#define PY_SSIZE_T_CLEAN
#include "types.h" 
#include <Python.h>

typedef struct PyMvBasicArray PyMvBasicArray;

//typedef int (*mvarrayiternextfunc)(PyMvArrayIter*,Py_ssize_t);

// Holds the data for a multivector array, to use in the array iterator
typedef struct PyMvBasicArray{   
    void *data; // The current iteration data
    void *data0; // The first element
    Py_ssize_t basic_size;
    Py_ssize_t ns; // Size of strides
    Py_ssize_t *strides; // Has ns + 1 elements, 1st element is the total size of the array
}PyMvBasicArray;

typedef struct PyMultipleArrayIter{
    PyMvBasicArray *arrays; 
    Py_ssize_t ns; // The number of symbols
    Py_ssize_t nm; // The number of arrays
    Py_ssize_t *index; // An index for symbol
    Py_ssize_t *shapes; // A shape for symbol
    Py_ssize_t dim; // The last iterated dimension
    char dflag; // flag for when dimensions are skipped 
}PyMultipleArrayIter;

struct ListGraph;

typedef struct ListGraph{
    struct ListGraph *parent;
    PyObject *element;
    PyObject *self;
    Py_ssize_t index;
}ListGraph;
typedef ScalarMultivector (*scalarop)(ScalarMultivector);
typedef ScalarMultivector (*scalarcomp)(ScalarMultivector,ScalarMultivector);

typedef int (*mixedoperation)(void*, PyMultivectorIter *, PyMultivectorIter * , PyAlgebraObject *,void *);

extern PyTypeObject PyMultivectorArrayType;
PyObject *algebra_multivector_array(PyAlgebraObject *self, PyObject *args, PyObject *kwds);

PyMultivectorObject *new_multivector_array(PyAlgebraObject *GA, char *type,  Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes);

PyMultivectorObject *new_multivector(PyAlgebraObject *ga, char *type);
PyMultivectorObject *new_multivector_inherit_type(PyAlgebraObject *GA, PyMultivectorSubType *type);

int get_multivector_type_table(PyAlgebraObject *ga, char *name, PyMultivectorSubType **subtype);
int check_multivector_mixed_type_table(PyMultivectorObject *mv, char *name);

//void free_multivector(PyMultivectorObject *self);
void free_multivector_data(PyMultivectorObject self);
Py_ssize_t parse_list_as_values(PyObject *values, ga_float **values_float);
int parse_list_as_bitmaps(PyObject *blades, int **bitmap);
int parse_list_as_basis_grades(PyAlgebraObject ga, int *grades, int **bitmaps, Py_ssize_t gsize);
//int parse_list_as_multivectors_1(PyObject *basis, ga_float **values, int **bitmaps);

//PyMultivectorObject *init_multivector(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga, PyTypeObject *obj_type, int type);
/*
// type methods
void multivector_dealloc(PyMultivectorObject *self);
PyObject *multivector_repr(PyMultivectorObject *self);
PyObject *multivector_grade_project(PyMultivectorObject *self, PyObject *args, PyObject *kwds);
PyObject *multivector_invert(PyMultivectorObject *self);
// number methods
PyObject *multivector_geometric_product(PyObject *left, PyObject *right);
PyObject *multivector_inner_product(PyObject *left, PyObject *right);
PyObject *multivector_outer_product(PyObject *left, PyObject *right);
PyObject *multivector_regressive_product(PyObject *left, PyObject *right);
PyObject *multivector_add(PyObject *left, PyObject *right);
PyObject *multivector_subtract(PyObject *left, PyObject *right);
PyObject *multivector_negative(PyMultivectorObject *self);
PyObject *multivector_positive(PyMultivectorObject *self);
// tp_methods class methods
PyObject* multivector_atomic_add(PyObject *self, PyObject *args);
PyObject* multivector_atomic_geometric_product(PyObject *self, PyObject *args);
PyObject* multivector_atomic_outer_product(PyObject *self, PyObject *args);
PyObject* multivector_exponential(PyObject *self, PyObject *args);
// other methods
PyObject* multivector_dual(PyMultivectorObject *self, PyObject *args);
PyObject* multivector_undual(PyMultivectorObject *self, PyObject *args);
PyObject* multivector_list(PyMultivectorObject *self, PyObject *args, PyObject *kwds);
PyObject* multivector_grade(PyMultivectorObject *self, PyObject *args);
PyObject* multivector_cast(PyMultivectorObject *self, PyObject *args);
*/


Py_ssize_t parse_list_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
//Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades);
char *bitmap_to_string(int bitmap);
int get_multivector_basis(PyAlgebraObject *self, PyObject *grades, PyMultivectorObject ***multivectors, ga_float **values, Py_ssize_t *mvsize);

PyMultivectorIter *init_multivector_iter(PyMultivectorObject *data, Py_ssize_t size);
void free_multivector_iter(PyMultivectorIter *iter, Py_ssize_t size);
//PyObject *type_iter_repr(PyMultivectorIter *iter, PrintTypeMV ptype, Py_ssize_t dsize);
//char *type_iter_repr_1(PyMultivectorIter *iter, PrintTypeMV ptype, Py_ssize_t dsize);

#endif
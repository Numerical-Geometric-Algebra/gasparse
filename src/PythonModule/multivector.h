#ifndef MULTIVECTOR_H_
#define MULTIVECTOR_H_
#define PY_SSIZE_T_CLEAN
#include "common.h"

typedef enum {
  MultivectorTypeMIN=-1,
  MultivectorType_sparse,
  MultivectorType_dense,
  MultivectorType_blades,
  MultivectorTypeMAX} MultivectorType;

// multivector iterator definitions
typedef struct _PyMultivectorIter PyMultivectorIter;
typedef struct _PyMultivectorObject PyMultivectorObject;

typedef int (*gaiternextfunc)(PyMultivectorIter *iter);
typedef PyMultivectorIter (*gaiterinitfunc)(PyMultivectorObject *data);
typedef struct _PyMultivectorIter{
    gaiternextfunc next;
    void *data;
    Py_ssize_t *index;
    Py_ssize_t size;
    Py_ssize_t niters; // number of iterations
    int bitmap;
    ga_float value;
    Py_ssize_t grade;
    int type;
}PyMultivectorIter;

typedef struct GradeProjectMap{
  int *grades0;
  int *grades1;
  int *grades;
  Py_ssize_t size0;
  Py_ssize_t size1;
  Py_ssize_t size;
}GradeProjectMap;

typedef PyMultivectorObject *(*gabinaryfunc)(PyMultivectorObject *left, PyMultivectorObject *right);
typedef PyMultivectorObject *(*gaaddfunc)(PyMultivectorObject *left, PyMultivectorObject *right, int sign);
typedef PyMultivectorObject *(*gaunarygradefunc)(PyMultivectorObject *self, int *grades, Py_ssize_t size);
typedef PyMultivectorObject *(*gaprodfunc)(PyMultivectorObject *left,PyMultivectorObject *right, ProductType ptype);
typedef PyMultivectorObject *(*gaternaryprodfunc)(PyMultivectorObject *data0,PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype);
typedef PyMultivectorObject *(*gaunaryfunc)(PyMultivectorObject *self);
typedef PyMultivectorObject *(*gaatomicfunc)(PyMultivectorObject *data, Py_ssize_t size);
typedef PyMultivectorObject *(*gaatomicprodfunc)(PyMultivectorObject *data, Py_ssize_t size,ProductType ptype);
typedef PyMultivectorObject *(*gacastfunc)(PyMultivectorObject *self);
typedef PyMultivectorObject *(*gabinarygradefunc)(PyMultivectorObject *left, PyMultivectorObject *right, ProductType ptype, GradeProjectMap gpmap);
typedef PyMultivectorObject *(*gascalarfunc)(PyMultivectorObject *self, ga_float value);
typedef PyMultivectorObject *(*gascalaraddfunc)(PyMultivectorObject *self, ga_float value, int sign);

typedef void *(*gainitfunc)(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga);
typedef void (*gafreefunc)(void *data);
typedef PyObject *(*gareprfunc)(void *data, PrintTypeMV ptype);


typedef struct PyMultivectorMath_Funcs{
    gaatomicfunc atomic_add;
    gaatomicprodfunc atomic_product;
    gaaddfunc add;
    gaprodfunc product;
    gabinarygradefunc graded_product;
    gaunarygradefunc grade_project;
    gascalarfunc scalar_product;
    gascalaraddfunc scalar_add;
    gaunaryfunc reverse;
    gaternaryprodfunc ternary_product;
}PyMultivectorMath_Funcs;

typedef struct PyMultivectorMixedMath_Funcs{
    gaaddfunc add;
    gaprodfunc product;
    gaatomicfunc atomic_add;
    gaatomicprodfunc atomic_product;
}PyMultivectorMixedMath_Funcs;


typedef struct PyMultivectorData_Funcs{
    gabinaryfunc copy;
    gafreefunc free;
    gainitfunc init;
    gacastfunc to_sparse;
    gacastfunc to_dense;
    gaiternextfunc iter_next;
    gaiterinitfunc iter_init;
}PyMultivectorData_Funcs;

typedef struct _PyMultivectorSubType{
    PyMultivectorMath_Funcs math_funcs;
    PyMultivectorData_Funcs data_funcs;
    int generated;
    int metric[32];
    char name[32];
    char type_name[32];
    Py_ssize_t metric_size;
    int ntype;// this is the type identifier use this to compare types
}PyMultivectorSubType;

extern PyMultivectorSubType multivector_subtypes_array[3];

typedef struct _PyMultivectorObject{
    PyObject_HEAD
    void *data;
    PyMultivectorMixedMath_Funcs mixed;
    PyAlgebraObject *GA;
    PyMultivectorSubType type;
}PyMultivectorObject;



typedef struct SparseMultivector{
    int *bitmap;
    ga_float *value;
    Py_ssize_t size;
}SparseMultivector;

typedef struct BladesMultivector{
    SparseMultivector *data;
    Py_ssize_t *grade;
    Py_ssize_t size;
}BladesMultivector;

typedef struct DenseMultivector{
    ga_float *value; // fixed size array for all the algebra
    Py_ssize_t size;
}DenseMultivector;

PyMultivectorObject *new_multivector(PyMultivectorObject *old, MultivectorType type);
void free_multivector(PyMultivectorObject *self);
PyMultivectorObject *init_multivector(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga, PyTypeObject *obj_type, int type);

// type methods
void multivector_dealloc(PyMultivectorObject *self);
PyObject *multivector_repr(PyMultivectorObject *self);
PyObject *multivector_grade_project(PyMultivectorObject *self, PyObject *args, PyObject *kwds);
PyObject *multivector_invert(PyMultivectorObject *self);
// number methods
PyObject *multivector_geometric_product(PyObject *left, PyObject *right);
PyObject *multivector_inner_product(PyObject *left, PyObject *right);
PyObject *multivector_outer_product(PyObject *left, PyObject *right);
PyObject *multivector_add(PyObject *left, PyObject *right);
PyObject *multivector_subtract(PyObject *left, PyObject *right);
PyObject *multivector_negative(PyMultivectorObject *self);
PyObject *multivector_positive(PyMultivectorObject *self);
// tp_methods
PyObject* multivector_atomic_add(PyObject *self, PyObject *args);
PyObject* multivector_atomic_geometric_product(PyObject *self, PyObject *args);
PyObject* multivector_atomic_outer_product(PyObject *self, PyObject *args);




#endif // MULTIVECTOR_H_

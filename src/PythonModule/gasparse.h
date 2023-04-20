#ifndef GASPARSE_H_
#define GASPARSE_H_
#define PY_SSIZE_T_CLEAN
#include "common.h"

#define Metric_SIZE(s) (s->p + s->q + s->r)
#define Is_NonZero(s) (s <= '9' && s >= '1')
#define MAX_SHAPE_SIZE 10

typedef enum {
  PrintTypeMIN=-1,
  PrintType_metric_array,
  PrintType_metric,
  PrintType_vectors,
  PrintTypeMAX} PrintType;

typedef enum {
    ProductTypeMIN=-1,
    ProductType_geometric,
    ProductType_inner,
    ProductType_outer,
    ProductType_regressive,
    ProductType_geometricinverted,
    ProductType_innerinverted,
    ProductType_outerinverted,
    ProductTypeMAX} ProductType;

typedef enum {
  PrintTypeMVMIN=-1,
  PrintTypeMV_reduced,
  PrintTypeMV_normal,
  PrintTypeMVMAX} PrintTypeMV;


typedef struct CliffordMap{
    char **sign;
    Py_ssize_t **bitmap;
    Py_ssize_t size;
}CliffordMap;


typedef struct GradeMap{ // used to map bitmap to position and grade
    Py_ssize_t *grade;
    Py_ssize_t *position;
    Py_ssize_t *grade_size; // number of basis vectors by grade
    Py_ssize_t max_grade;
    Py_ssize_t size;
}GradeMap;

typedef struct PyAlgebraObject{
    PyObject_HEAD
    GradeMap gm;
    CliffordMap product[ProductTypeMAX];
    Py_ssize_t p,q,r;
    char *metric;
    PrintType print_type;
    PrintTypeMV print_type_mv;
    ga_float precision;
}PyAlgebraObject;

/*
typedef struct PyArrayStrides{
    Py_ssize_t **strides; // n_tensors x n_symbols
    Py_ssize_t *n_strides;
    Py_ssize_t n_tensors;
    Py_ssize_t n_symbols;
}PyArrayStrides;

typedef struct PyGaArrayIterator{
    PyArrayStrides strides_s;
    void **data;
    int *depth;
    Py_ssize_t *index;
    Py_ssize_t sizeof_data;
}PyGaArrayIterator;
typedef struct _PyMultivectorArrayObject{
    PyObject_HEAD
    void *data;
    Py_ssize_t *shapes;
    Py_ssize_t *strides;
    Py_ssize_t shape_size;
    Py_ssize_t data_size;
    PyAlgebraObject *GA;
    PyMultivectorMath_Funcs math_funcs;
    PyMultivectorData_Funcs data_funcs;
    MultivectorType type;
}PyMultivectorArrayObject;
*/

#endif // GASPARSE_H_

#ifndef GASPARSE_H_
#define GASPARSE_H_
#define PY_SSIZE_T_CLEAN

#define Metric_SIZE(s) (s->p + s->q + s->r)
#define Is_NonZero(s) (s <= '9' && s >= '1')
#define MAX_SHAPE_SIZE 10

typedef enum {MultivectorTypeMIN=-1,MultivectorType_sparse,MultivectorType_dense,MultivectorType_blades,MultivectorType_scalar,MultivectorTypeMAX} MultivectorType;
typedef enum {PrintTypeMIN=-1,PrintType_metric_array,PrintType_metric,PrintType_vectors,PrintTypeMAX} PrintType;
typedef enum {ProductTypeMIN=-1,ProductType_geometric,ProductType_inner,ProductType_outer,ProductType_regressive,ProductType_inverted,ProductTypeMAX} ProductType;


typedef struct CliffordMap{
    char **sign;
    Py_ssize_t **bitmap;
    Py_ssize_t size;
}CliffordMap;


typedef struct DenseGradeMap{
    Py_ssize_t max_grade;
    Py_ssize_t grade_size;
    Py_ssize_t *grade; // grade of each basis blade
}DenseGradeMap;

typedef struct GradeMap{ // used to map bitmap to position and grade
    Py_ssize_t *grade;
    Py_ssize_t *position;
    Py_ssize_t *grade_size; // number of basis vectors by grade
    Py_ssize_t max_grade;
    Py_ssize_t size;
}GradeMap;

typedef struct PyGeometricAlgebraObject{
    PyObject_HEAD
    DenseGradeMap dgm;
    GradeMap gm;
    CliffordMap product[ProductTypeMAX];
    Py_ssize_t p,q,r;
    char *metric;
    PrintType print_type;
    float precision;
}PyGeometricAlgebraObject;

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


typedef struct _PyMultivectorObject PyMultivectorObject;


typedef void *(*gacastfunc)(void *data, PyGeometricAlgebraObject *ga);
typedef void *(*gaatomicfunc)(void *data, Py_ssize_t size, PyGeometricAlgebraObject *ga);
typedef PyMultivectorObject *(*gabinaryfunc)(PyMultivectorObject *left,PyMultivectorObject *right);
typedef PyMultivectorObject *(*gaproductfunc)(PyMultivectorObject *left,PyMultivectorObject *right, ProductType ptype);
typedef void *(*gaallocfunc)(PyGeometricAlgebraObject *ga);
typedef void *(*gabinaryextrafunc)(void *left, void *right, GradeMap grades, PyGeometricAlgebraObject *ga);
typedef void *(*gacopyfunc)(void *dest, void *src, PyGeometricAlgebraObject *ga);
typedef void *(*gainitfunc)(int *bitmap, float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size);
typedef void (*gafreefunc)(void *data);
typedef PyObject *(*gareprfunc)(void *data);

typedef struct PyMultivector_Funcs{
    gaatomicfunc atomic_add[MultivectorTypeMAX][MultivectorTypeMAX];
    gabinaryfunc add[MultivectorTypeMAX][MultivectorTypeMAX];
    gaproductfunc product[MultivectorTypeMAX][MultivectorTypeMAX];
    gacastfunc cast[MultivectorTypeMAX][MultivectorTypeMAX];
    gabinaryextrafunc graded_product[MultivectorTypeMAX][MultivectorTypeMAX];
    gacopyfunc copy[MultivectorTypeMAX];
    gaallocfunc alloc_dense[MultivectorTypeMAX];
    gaallocfunc alloc_sparse[MultivectorTypeMAX];
    gafreefunc free[MultivectorTypeMAX];
    gainitfunc init[MultivectorTypeMAX];
    gacastfunc to_sparse[MultivectorTypeMAX];
    gacastfunc to_dense[MultivectorTypeMAX];
    gareprfunc repr[MultivectorTypeMAX];
}PyMultivector_Funcs;


typedef struct PyMultivectorArrayObject{
    PyObject_HEAD
    void *data;
    Py_ssize_t *shapes;
    Py_ssize_t *strides;
    Py_ssize_t shape_size;
    Py_ssize_t data_size;
    PyGeometricAlgebraObject *GA;
    PyMultivector_Funcs funcs;
    MultivectorType type;
}PyMultivectorArrayObject;

typedef struct _PyMultivectorObject{
    PyObject_HEAD
    void *data;
    PyGeometricAlgebraObject *GA;
    PyMultivector_Funcs funcs; // These functions should be different from PyMultivectorArrayObject
    MultivectorType type;
}PyMultivectorObject;



/*
  Py_TYPE macro shouldn't be used in these type definitions
  PyObject_HEAD is used only to do reference counting with the Py_REFCNT macro
  I should implement my own referent counting for these types but I don't know how reference counting is implemented internally
*/

typedef struct SparseMultivector{
    PyObject_HEAD
    int *bitmap;
    float *value;
    Py_ssize_t size;
}SparseMultivector;

typedef struct BladesMultivector{
    PyObject_HEAD
    SparseMultivector *data;
    Py_ssize_t *grade;
    Py_ssize_t size;
}BladesMultivector;

typedef struct DenseMultivector{
    PyObject_HEAD
    float *value; // fixed size array for all the algebra
}DenseMultivector;


PyMultivectorObject *init_multivector(int *bitmap, float *value, Py_ssize_t size, PyGeometricAlgebraObject *ga, PyTypeObject *obj_type, int type);

// type methods
void multivector_dealloc(PyMultivectorObject *self);
PyObject *multivector_repr(PyMultivectorObject *self);

// number methods
PyObject *multivector_geometric_product(PyObject *left, PyObject *right);

Py_ssize_t grade(Py_ssize_t v);

#endif // GASPARSE_H_

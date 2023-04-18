#ifndef GASPARSE_H_
#define GASPARSE_H_
#define PY_SSIZE_T_CLEAN

#define Metric_SIZE(s) (s->p + s->q + s->r)
#define Is_NonZero(s) (s <= '9' && s >= '1')
#define GRADE(value) (__builtin_popcountll(value))
#define MAX_SHAPE_SIZE 10

typedef enum {MultivectorTypeMIN=-1,MultivectorType_sparse,MultivectorType_dense,MultivectorType_blades,MultivectorTypeMAX} MultivectorType;
typedef enum {PrintTypeMIN=-1,PrintType_metric_array,PrintType_metric,PrintType_vectors,PrintTypeMAX} PrintType;
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
typedef enum {PrintTypeMVMIN=-1,PrintTypeMV_reduced,PrintTypeMV_normal,PrintTypeMVMAX} PrintTypeMV;

typedef double ga_float;

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

typedef PyMultivectorObject *(*gabinaryfunc)(PyMultivectorObject *left, PyMultivectorObject *right);
typedef PyMultivectorObject *(*gaaddfunc)(PyMultivectorObject *left, PyMultivectorObject *right, int sign);
typedef PyMultivectorObject *(*gaunarygradefunc)(PyMultivectorObject *self, int *grades, Py_ssize_t size);
typedef PyMultivectorObject *(*gaprodfunc)(PyMultivectorObject *left,PyMultivectorObject *right, ProductType ptype);
typedef PyMultivectorObject *(*gaternaryprodfunc)(PyMultivectorObject *data0,PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype);
typedef PyMultivectorObject *(*gaunaryfunc)(PyMultivectorObject *self);
typedef PyMultivectorObject *(*gaatomicfunc)(PyMultivectorObject *data, Py_ssize_t size);
typedef PyMultivectorObject *(*gaatomicprodfunc)(PyMultivectorObject *data, Py_ssize_t size,ProductType ptype);
typedef PyMultivectorObject *(*gacastfunc)(PyMultivectorObject *self);
typedef PyMultivectorObject *(*gabinarygradefunc)(PyMultivectorObject *left, PyMultivectorObject *right, GradeMap grades);
typedef PyMultivectorObject *(*gascalarfunc)(PyMultivectorObject *self, ga_float value);
typedef PyMultivectorObject *(*gascalaraddfunc)(PyMultivectorObject *self, ga_float value, int sign);

typedef void *(*gainitfunc)(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga);
typedef void (*gafreefunc)(void *data);
typedef PyObject *(*gareprfunc)(void *data, PrintTypeMV ptype);

typedef struct PyMultivectorMath_Funcs{
    gaatomicfunc atomic_add[MultivectorTypeMAX];
    gaatomicprodfunc atomic_product[MultivectorTypeMAX];
    gaaddfunc add[MultivectorTypeMAX];
    gaprodfunc product[MultivectorTypeMAX];
    gabinarygradefunc graded_product[MultivectorTypeMAX];
    gaunarygradefunc grade_project[MultivectorTypeMAX];
    gascalarfunc scalar_product[MultivectorTypeMAX];
    gascalaraddfunc scalar_add[MultivectorTypeMAX];
    gaunaryfunc reverse[MultivectorTypeMAX];
    gaternaryprodfunc ternary_product[MultivectorTypeMAX];
    gaaddfunc mixed_add;
    gaprodfunc mixed_product;
    gaatomicfunc mixed_atomic_add;
    gaatomicprodfunc mixed_atomic_product;
}PyMultivectorMath_Funcs;


typedef struct PyMultivectorData_Funcs{
    gacastfunc cast[MultivectorTypeMAX][MultivectorTypeMAX];
    gabinaryfunc copy[MultivectorTypeMAX];
    gafreefunc free[MultivectorTypeMAX];
    gainitfunc init[MultivectorTypeMAX];
    gacastfunc to_sparse[MultivectorTypeMAX];
    gacastfunc to_dense[MultivectorTypeMAX];
}PyMultivectorData_Funcs;


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

typedef struct _PyMultivectorObject{
    PyObject_HEAD
    void *data;
    PyAlgebraObject *GA;
    PyMultivectorMath_Funcs math_funcs;
    PyMultivectorData_Funcs data_funcs;
    MultivectorType type;
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

// multivector iterator definitions
typedef struct _PyMultivectorIter PyMultivectorIter;
typedef int (*gaiternextfunc)(PyMultivectorIter *iter);
typedef struct _PyMultivectorIter{
    gaiternextfunc next;
    void *data;
    Py_ssize_t *index;
    Py_ssize_t size;
    int bitmap;
    ga_float value;
    Py_ssize_t grade;
    MultivectorType type;
}PyMultivectorIter;


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

// multivector methods
/*
// initialize empty multivectors
static SparseMultivector init_sparse_empty(Py_ssize_t size);
static DenseMultivector init_dense_empty(Py_ssize_t size);
static BladesMultivector init_blades_empty(Py_ssize_t size);

// convert dense representations to sparse
static SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense_y, Py_ssize_t size);
static BladesMultivector sparse_dense_to_blades_sparse(SparseMultivector dense, GradeMap gm);

// remove small values of a multivector
static void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size);

// initialize multivectors given an array of values and bitmaps
static SparseMultivector sparse_init_(int *bitmap, ga_float *value, Py_ssize_t size);
static BladesMultivector blades_init_(int *bitmap, ga_float *value, Py_ssize_t size, GradeMap gm);
static DenseMultivector dense_init_(int *bitmap, ga_float *value, Py_ssize_t size, Py_ssize_t algebra_size);

// free multivectors
static void sparse_free_(SparseMultivector sparse);
static void blades_free_(BladesMultivector blades);
static void dense_free_(DenseMultivector dense);
*/
// determine the grade of a bitmap

#endif // GASPARSE_H_

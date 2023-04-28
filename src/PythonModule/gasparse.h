#ifndef GASPARSE_H_
#define GASPARSE_H_
#define PY_SSIZE_T_CLEAN
#include "common.h"

#define METRIC_SIZE(s) (s->p + s->q + s->r)
#define MAX_GRADE(s) (s->p + s->q + s->r)
#define IS_NONZERO(s) (s <= '9' && s >= '1')
#define ABS(value) (((value) < 0) ? -(value): (value))
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
    ProductType_regressiveinverted,
    ProductTypeMAX} ProductType;

typedef enum {
  PrintTypeMVMIN=-1,
  PrintTypeMV_reduced,
  PrintTypeMV_normal,
  PrintTypeMVMAX} PrintTypeMV;

typedef enum {
  MultivectorTypeMIN=-1,
  MultivectorType_sparse,
  MultivectorType_dense,
  MultivectorType_blades,
  MultivectorTypeMAX} MultivectorType;

typedef enum {
  ComputationModeMIN = -1,
  ComputationMode_generated,
  ComputationMode_large,
  ComputationMode_generic,
  ComputationMode_devgeneration,
  ComputationModeMAX} ComputationMode;

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

typedef struct DualMap{
    char *sign;
    Py_ssize_t *bitmap;
    Py_ssize_t size;
}DualMap;

typedef struct PyMultivectorSubType PyMultivectorSubType;
typedef struct PyAlgebraObject PyAlgebraObject;
typedef struct PyMultivectorIter PyMultivectorIter;
typedef struct PyMultivectorObject PyMultivectorObject;
typedef struct PyMultivectorMixedMath_Funcs PyMultivectorMixedMath_Funcs;

typedef struct MultivectorDefaults{
    char *type_name;
}MultivectorDefaults;

typedef struct PyAlgebraObject{
    PyObject_HEAD
    GradeMap gm;
    DualMap dm;
    CliffordMap product[ProductTypeMAX];
    Py_ssize_t p,q,r;
    char *metric;
    PrintType print_type;
    PrintTypeMV print_type_mv;
    ga_float precision;
    PyMultivectorSubType *types;
    PyMultivectorMixedMath_Funcs *mixed;
    Py_ssize_t tsize; // number of types
    Py_ssize_t asize; // sizeof the algebra
    MultivectorDefaults mdefault;
}_PyAlgebraObject;



typedef int (*gaiternextfunc)(PyMultivectorIter *iter);
typedef PyMultivectorIter (*gaiterinitfunc)(PyMultivectorObject *data);
typedef struct PyMultivectorIter{
    gaiternextfunc next;
    void *data;
    Py_ssize_t *index;
    Py_ssize_t size;
    Py_ssize_t niters; // number of iterations
    int bitmap;
    ga_float value;
    Py_ssize_t grade;
    int type;
    char *type_name;
}_PyMultivectorIter;

typedef struct GradeProjectMap{
  int *grades0;
  int *grades1;
  int *grades;
  Py_ssize_t size0;
  Py_ssize_t size1;
  Py_ssize_t size;
}GradeProjectMap;

/* typedef int (*compbitmapfunc)(Py_ssize_t,Py_ssize_t,Py_ssize_t); */


typedef PyMultivectorObject *(*gabinaryfunc)(PyMultivectorObject *left, PyMultivectorObject *right);
typedef PyMultivectorObject *(*gaaddfunc)(PyMultivectorObject *left, PyMultivectorObject *right, int sign);
typedef PyMultivectorObject *(*gaunarygradefunc)(PyMultivectorObject *self, int *grades, Py_ssize_t size);
typedef PyMultivectorObject *(*gaprodfunc)(PyMultivectorObject *left,PyMultivectorObject *right, ProductType ptype);
typedef PyMultivectorObject *(*gaternaryprodfunc)(PyMultivectorObject *data0,PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype);
typedef PyMultivectorObject *(*gaunaryfunc)(PyMultivectorObject *self);
typedef PyMultivectorObject *(*gaatomicfunc)(PyMultivectorObject *data, Py_ssize_t size);
typedef PyMultivectorObject *(*gaatomicprodfunc)(PyMultivectorObject *data, Py_ssize_t size,ProductType ptype);
typedef PyMultivectorObject *(*gacastfunc)(PyMultivectorObject *from, PyMultivectorObject *to);
typedef PyMultivectorObject *(*gabinarygradefunc)(PyMultivectorObject *left, PyMultivectorObject *right, ProductType ptype, GradeProjectMap gpmap);
typedef PyMultivectorObject *(*gascalarfunc)(PyMultivectorObject *self, ga_float value);
typedef PyMultivectorObject *(*gascalaraddfunc)(PyMultivectorObject *self, ga_float value, int sign);

typedef void *(*gainitfunc)(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga);
typedef void (*gafreefunc)(void *data);

typedef PyMultivectorObject *(*gamixedaddfunc)(PyMultivectorObject *left, PyMultivectorObject *right, PyMultivectorObject *defdata, int sign);
typedef PyMultivectorObject *(*gamixedprodfunc)(PyMultivectorObject *left,PyMultivectorObject *right, PyMultivectorObject *defdata, ProductType ptype);
typedef PyMultivectorObject *(*gamixedatomicfunc)(PyMultivectorObject *data, Py_ssize_t size, PyMultivectorObject *defdata);
typedef PyMultivectorObject *(*gamixedatomicprodfunc)(PyMultivectorObject *data, Py_ssize_t size, PyMultivectorObject *defdata, ProductType ptype);

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
    gaunaryfunc dual;
    gaunaryfunc undual;
    gaternaryprodfunc ternary_product;
}PyMultivectorMath_Funcs;

typedef struct PyMultivectorMixedMath_Funcs{
    gamixedaddfunc add;
    gamixedprodfunc product;
    gamixedatomicfunc atomic_add;
    gamixedatomicprodfunc atomic_product;
    char *type_names[]; // null terminated
}_PyMultivectorMixedMath_Funcs;


typedef struct PyMultivectorData_Funcs{
    gabinaryfunc copy;
    gafreefunc free;
    gainitfunc init;
    gacastfunc cast;
    gaiternextfunc iter_next;
    gaiterinitfunc iter_init;
}PyMultivectorData_Funcs;

typedef struct PyMultivectorSubType{
    PyMultivectorMath_Funcs *math_funcs;
    PyMultivectorData_Funcs *data_funcs;
    int generated;
    char metric[32];
    char name[32];
    char type_name[32];
    Py_ssize_t msize; // metric array size number of basis vector
    Py_ssize_t asize; // algebra size number of basis blades
    int ntype;// this is the type identifier use this to compare types
}_PyMultivectorSubType;

// mixed type operation tables
extern PyMultivectorMixedMath_Funcs multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs cast_multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs largemultivector_mixed_fn;

// same type operation tables
extern PyMultivectorSubType largemultivector_subtypes_array[3];
extern PyMultivectorSubType multivector_subtypes_array[3];

void fill_missing_funcs(void);

typedef struct PyMultivectorObject{
    PyObject_HEAD
    void *data;
    PyMultivectorMixedMath_Funcs *mixed;
    PyAlgebraObject *GA;
    PyMultivectorSubType type;
}_PyMultivectorObject;



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

PyMultivectorObject *new_multivectorbyname(PyMultivectorObject *old, char *name);
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
PyObject *multivector_regressive_product(PyObject *left, PyObject *right);
PyObject *multivector_add(PyObject *left, PyObject *right);
PyObject *multivector_subtract(PyObject *left, PyObject *right);
PyObject *multivector_negative(PyMultivectorObject *self);
PyObject *multivector_positive(PyMultivectorObject *self);
// tp_methods class methods
PyObject* multivector_atomic_add(PyObject *self, PyObject *args);
PyObject* multivector_atomic_geometric_product(PyObject *self, PyObject *args);
PyObject* multivector_atomic_outer_product(PyObject *self, PyObject *args);
// other methods
PyObject* multivector_dual(PyMultivectorObject *self, PyObject *args);
PyObject* multivector_undual(PyMultivectorObject *self, PyObject *args);

// TYPE methods
int comp_abs(ga_float v, ga_float p);
int ga_check_value_precision(PyAlgebraObject *ga, ga_float v);
SparseMultivector init_sparse_empty(Py_ssize_t size);
DenseMultivector init_dense_empty(Py_ssize_t size);
BladesMultivector init_blades_empty(Py_ssize_t size);
void sparse_free_(SparseMultivector sparse);
void blades_free_(BladesMultivector blades);
void dense_free_(DenseMultivector dense);

Py_ssize_t parse_list_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size);
SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense, Py_ssize_t size);
Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades);
char *bitmap_to_string(int bitmap);

PyMultivectorIter *init_multivector_iter(PyMultivectorObject *data, Py_ssize_t size);
void free_multivector_iter(PyMultivectorIter *iter, Py_ssize_t size);

#endif // GASPARSE_H_

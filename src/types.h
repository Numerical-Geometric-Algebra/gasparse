#ifndef TYPES_H_
#define TYPES_H_
#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*
    Put here all the defines and all the typedefs of all the code
*/

#define METRIC_SIZE(s) (s->p + s->q + s->r)
#define MAX_GRADE(s) (s->p + s->q + s->r)
#define IS_NONZERO(s) (s <= '9' && s >= '1')
#define ABS(value) (((value) < 0) ? -(value): (value))
#define MAX_SHAPE_SIZE 10
#define GRADE(value) (__builtin_popcountll(value))
#define INDEX_DATA(s,i) ((s)->data + (i)*(s)->type->basic_size)

typedef double ga_float;

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
  MultivectorType_sparsevector,
  MultivectorType_densevector,
  MultivectorType_scalar,
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
    Py_ssize_t *grade_size; // number of basis multivectors by grade
    Py_ssize_t max_grade;
    Py_ssize_t size;
}GradeMap;

typedef struct GradeTable{ // For each grade stores the array of bitmaps of that grade
    Py_ssize_t **bitmaps;
    Py_ssize_t *grade_size; // number of basis elements by grade
    Py_ssize_t size;
}GradeTable;


typedef struct DualMap{
    char *sign;
    Py_ssize_t *bitmap;
    Py_ssize_t size;
}DualMap;

typedef struct SparseMapBase{
    Py_ssize_t position[4];
    char sign;
}SparseMapBase;

typedef struct SparseMap{// Map for the graded ternary products
    SparseMapBase *base;
    Py_ssize_t size;
}SparseMap;

// typedef struct _BasisElement BasisElement;

typedef struct BasisElement{
    int bitmap;
    ga_float value;
    struct BasisElement *next;
}BasisElement;

typedef SparseMap ****SparseTernaryMap;

// Structure for the base of the ternary map
typedef struct MapBase{
    Py_ssize_t position; // integer value representing the basis element
    Py_ssize_t grade; // The resulting grade
    char sign; // Sign of the resulting product (-1,0,1)
}MapBase;

typedef MapBase ******TernaryMap;
typedef MapBase ***PositionMap;

typedef struct PyMultivectorSubType PyMultivectorSubType;
typedef struct PyAlgebraObject PyAlgebraObject;
typedef struct PyMultivectorIter PyMultivectorIter;
typedef struct PyMultivectorObject PyMultivectorObject;
typedef struct PyMultivectorMixedMath_Funcs PyMultivectorMixedMath_Funcs;
typedef PyMultivectorObject PyMvObject;
typedef PyMultivectorSubType PyMvSubType;

typedef struct MultivectorDefaults{
    char *type_name;
}MultivectorDefaults;

typedef struct PyAlgebraObject{
    PyObject_HEAD
    GradeMap gm;
    GradeTable gt;
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

typedef struct PyMultivectorObject{
    PyObject_HEAD
    void *data;
    PyMultivectorMixedMath_Funcs *mixed;
    PyAlgebraObject *GA;
    PyMultivectorSubType *type;
    Py_ssize_t ns; // Size of shapes and strides
    Py_ssize_t *strides; // Has ns + 1 elements, the 1st element is the total nbr of mvs
    Py_ssize_t *shapes; // Has ns elements
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

typedef struct GradedDenseMultivector{// Use GradeMap to determine the size for each grade
    ga_float **values;
    Py_ssize_t *grades;
    Py_ssize_t size; // The number of different grades
}GradedDenseMultivector;

typedef ga_float ScalarMultivector;

typedef struct SparseVector{
    ga_float *value;
    int *basis; // The basis of the vectors e1,e2,e3,... as integers 1,2,3,...
    Py_ssize_t size;
}SparseVector;

typedef struct DenseVector{ // Fized size for a given geometric algebra 
    ga_float *value;
}DenseVector;

typedef GradedDenseMultivector GrdDenseMv;

typedef int (*gaiternextfunc)(PyMultivectorIter *iter);
//typedef PyMultivectorIter (*gaiterinitfunc)(PyMultivectorObject *data);
typedef PyMultivectorIter (*gaiterinitfunc)(void *data, PyMultivectorSubType *type);
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

// int (*operator)(outputs,inputs,GA,additional_args)

typedef int (*gainitfunc)(void*,PyAlgebraObject *ga, int *bitmap, ga_float *value, Py_ssize_t size);
typedef int (*gabinaryfunc)(void* out, void* left, void *right, PyAlgebraObject *GA);
typedef int (*gaaddfunc)(void* out, void* left, void *right, PyAlgebraObject *GA, int sign);
typedef int (*gaunarygradefunc)(void *out, void *self, PyAlgebraObject *GA, int *grades, Py_ssize_t size);
typedef int (*gaprodfunc)(void *out, void *left, void *right, PyAlgebraObject *GA, ProductType ptype);
typedef int (*gascalaraddfunc)(void* out, void *self, PyAlgebraObject *GA, ga_float value, int sign);
typedef int (*gaternaryprodfunc)(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *GA, ProductType ptype);
typedef int (*gaunaryfunc)(void *out, void *self, PyAlgebraObject *GA);
typedef int (*gaatomicfunc)(void *out,void *data, PyAlgebraObject *GA, Py_ssize_t size);
typedef int (*gaatomicprodfunc)(void *out, void *data, PyAlgebraObject *GA, Py_ssize_t size, ProductType ptype);
typedef int (*gabinarygradefunc)(void *out, void *left, void *right, PyAlgebraObject *GA, ProductType ptype, GradeProjectMap gpmap);
typedef int (*gascalarfunc)(void* out,void *self, PyAlgebraObject *GA, ga_float value);
typedef int (*gacastfunc)(PyMultivectorIter *from, void *to, PyAlgebraObject *GA);
typedef void (*gafreefunc)(void *data);

// Mixed type operations
typedef int (*gamixedaddfunc)(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, int sign);
typedef int (*gamixedprodfunc)(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, ProductType ptype);
typedef int (*gamixedatomicfunc)(void* out, PyMultivectorIter *iter, PyAlgebraObject *ga, Py_ssize_t dsize); 
typedef int (*gamixedatomicprodfunc)(void *out, PyMultivectorIter *iter, PyAlgebraObject *ga, Py_ssize_t size, ProductType ptype);

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
    gaunaryfunc dual;
    gaunaryfunc undual;
    gaunaryfunc exp;
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
    Py_ssize_t basic_size; // the size of the type
}_PyMultivectorSubType;

// mixed type operation tables
extern PyMultivectorMixedMath_Funcs multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs cast_multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs largemultivector_mixed_fn;

// same type operation tables
extern PyMultivectorSubType largemultivector_subtypes_array[3];
extern PyMultivectorSubType multivector_subtypes_array[4];
extern PyTypeObject PyMultivectorType;
void fill_missing_funcs(void);


PyObject *algebra_multivector(PyAlgebraObject *self, PyObject *args,PyObject *kwds);


#endif
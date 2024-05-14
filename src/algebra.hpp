#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include <Python.h>

typedef double ga_float;



namespace CliffordAlgebra {


class MultivectorType{
    // PyMultivectorMath_Funcs *math_funcs;
    // PyMultivectorData_Funcs *data_funcs;
    int generated;
    char metric[32];
    char name[32];
    char type_name[32];
    Py_ssize_t msize; // metric array size number of basis vector
    Py_ssize_t asize; // algebra size number of basis blades
    int ntype;// this is the type identifier use this to compare types
    Py_ssize_t basic_size; // the size of the type
};

// M is the multivector type
// T is the type of each component of the multivector

template<typename T, typename M>
class Algebra{

        typedef struct CliffordMap{
            char **sign;
            Py_ssize_t **bitmap;
            Py_ssize_t size;
        }CliffordMap;

private:
    typedef struct GradeMap{ // used to map bitmap to position and grade
        Py_ssize_t *grade;
        Py_ssize_t *position;
        Py_ssize_t *grade_size; // number of basis multivectors by grade
        Py_ssize_t max_grade;
        Py_ssize_t size;
    }GradeMap;

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
    PrintTypeMIN=-1,
    PrintType_metric_array,
    PrintType_metric,
    PrintType_vectors,
    PrintTypeMAX} PrintType;

    typedef enum {
    MultivectorTypeMIN=-1,
    MultivectorType_sparse,
    MultivectorType_dense,
    MultivectorType_blades,
    MultivectorType_sparsevector,
    MultivectorType_densevector,
    MultivectorType_scalar,
    MultivectorTypeMAX} MultivectorType;

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

    typedef struct MultivectorDefaults{
        char *type_name;
    }MultivectorDefaults;



    typedef struct MultivectorIter{
        int next();
        void *data;
        Py_ssize_t *index;
        Py_ssize_t size;
        Py_ssize_t niters; // number of iterations
        int bitmap;
        ga_float value;
        Py_ssize_t grade;
        int type;
        char *type_name;
    }_MultivectorIter;

    typedef struct GradeProjectMap{
    int *grades0;
    int *grades1;
    int *grades;
    Py_ssize_t size0;
    Py_ssize_t size1;
    Py_ssize_t size;
    }GradeProjectMap;

    GradeMap gm;
    GradeTable gt;
    DualMap dm;
    CliffordMap product_map[ProductTypeMAX];
    Py_ssize_t p,q,r;
    char *metric;
    PrintType print_type;
    PrintTypeMV print_type_mv;
    T precision;
    MultivectorType *types;
    // PyMultivectorMixedMath_Funcs *mixed;
    Py_ssize_t tsize; // number of types
    Py_ssize_t asize; // sizeof the algebra
    MultivectorDefaults mdefault;


    int mixed_add(M *out, MultivectorIter *iter0, MultivectorIter *iter1, int sign);
    int mixed_product(M *out, MultivectorIter *iter0, MultivectorIter *iter1, ProductType ptype);
    int mixed_atomic_add(M *out, MultivectorIter *iter, Algebra *ga, Py_ssize_t dsize); 
    int mixed_atomic_product(M *out, MultivectorIter *iter, Algebra *ga, Py_ssize_t size, ProductType ptype);
    char *mixed_type_names[]; // null terminated

    
    int init(M *out, int *bitmap, T *value, Py_ssize_t size);
    int add(M *out, M *left, M *right, int sign);
    int grade_project(M *out, M *self, int *grades, Py_ssize_t size);
    int product(M *out, M *left, M *right, ProductType ptype);
    int scalar_add(M *out, M *self, T value, int sign);
    int ternary_product(M *out, M *data0, M *data1, M *data2, ProductType ptype);
    int reverse(M *out, M *self);
    int dual(M *out, M *self);
    int undual(M *out, M *self);
    int exp(M *out, M *self);
    int atomic_add(M *out, M *data, Py_ssize_t size);
    int atomic_product(M *out, M *data, Py_ssize_t size, ProductType ptype);
    int graded_product(M *out, M *left, M *right, ProductType ptype, GradeProjectMap gpmap);
    int scalar_product(M* out,M *self, T value);


    void map_dealloc(CliffordMap*);
    void clifford_sub_algebra(Py_ssize_t k, char **s, int metric);


};

template <typename T,typename W>
class Multivector{

    Algebra<T,Multivector> *GA;
    W *MultivectorData;

    int alloc_mvarray_data();

    inline int init(Py_ssize_t i, int *bitmap, T *value, Py_ssize_t size);
    int add(Multivector *left, Multivector *right, int sign);
};

template <typename T>
struct SparseStruct{
    T *values;
    int *bitmap;
    Py_ssize_t size;
};

template <typename T>
struct DenseStruct{
    T* values;
};

template <typename T>
struct BladesStruct{
    SparseStruct<T> *data;
    Py_ssize_t *grade;
    Py_ssize_t size;
};

template <typename T>
class MultivectorSparse : public Multivector<T,SparseStruct<T>>{
    inline int init(Py_ssize_t i, int *bitmap, T *value, Py_ssize_t size);
};

template <typename T>
class MultivectorDense : Multivector<T,T*>{
    
};

template <typename T>
class MultivectorBlades : Multivector<T,BladesStruct<T>>{

};

template <typename T>
class MultivectorScalar : Multivector<T,T>{

};

// class MultivectorMixedMath{
//     // int add(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, Algebra *ga, int sign);
//     // gamixedaddfunc add;
//     // gamixedprodfunc product;
//     // gamixedatomicfunc atomic_add;
//     // gamixedatomicprodfunc atomic_product;
//     // char *type_names[]; // null terminated
// };

}
#endif

#ifndef MULTIVECTOR_GEN_H_
#define MULTIVECTOR_GEN_H_



typedef enum {
  gen_MultivectorTypeMIN = -1,
  gen_MultivectorType_dense,
  gen_MultivectorType_blades,
  gen_MultivectorTypeMAX} gen_MultivectorType;
/*
typedef struct gen_MathFuncs{
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
}gen_MathFuncs;

typedef struct gen_DataFuncs{
    gabinaryfunc copy;
    gafreefunc free;
    gainitfunc init;
    gacastfunc to_sparse;
    gacastfunc to_dense;
    gaiternextfunc iter_next;
    gaiterinitfunc iter_init;
}gen_DataFuncs;


typedef struct gen_AlgebraType{
    gen_MathFuncs math_funcs;
    gen_DataFuncs data_funcs;
    int generated;
    int metric[32];
    char name[32];
    char type_name[32];
    Py_ssize_t metric_size;
}gen_AlgebraType;
*/
extern PyMultivectorSubType gen_subtypes_array[4];

#endif // MULTIVECTOR_GEN_H_

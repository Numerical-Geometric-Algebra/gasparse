#ifndef MULTILINEAR_H_
#define MULTILINEAR_H_
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "types.h"

// typedef struct _PyOperatorObject PyOperatorObject;


// typedef int (*multilinearoperation)(void *out, void *data0, PyOperatorObject op);
typedef int (*gaspecialprodfunc)(void *out, void *data0, void *data1, void *data2, int *grades_out, Py_ssize_t gsize, PyAlgebraObject *ga, ProductType ptype1, ProductType ptype2);


typedef struct PyMultilinearMath_Funcs{
    // multilinearoperation multilinear_op; // The function that will call the multilinear operation
    gaspecialprodfunc ternary_prod; 
}PyMultilinearMath_Funcs;


typedef struct PyMultilinearObject{
    PyObject_HEAD
    PyAlgebraObject *GA;
    TernaryMap ternary_product[ProductTypeMAX][ProductTypeMAX];
    SparseTernaryMap sparse_ternary_product[ProductTypeMAX][ProductTypeMAX];
    PyMultilinearMath_Funcs f;
}PyMultilinearObject;


typedef enum OperatorFlags{
    OpFlag_gradepreserving = 1 << 0, // grade preserving flag
    OpFlag_precedence = 1 << 1, // right or left precedence, 0 for left, 1 for 
    OpFlag_reverse0 = 1 << 2, // reverse the argument of the multilinear function
    OpFlag_reverse1 = 1 << 3, // reverse the left_mv
    OpFlag_reverse2 = 1 << 4, // reverse the right_mv
    OpFlag_rightisleft = 1 << 5, // Ignore right and consider only the left mvs
    OpFlag_leftisright = 1 << 6, // Ignore left and consider only the right mvs
    OpFlag_binaryprod = 1 << 7, // Take a binary product 
}OperatorFlags;

typedef struct PyOperatorObject{
    PyObject_HEAD
    PyMultivectorObject left_mvs;
    PyMultivectorObject right_mvs;
    PyMultivectorObject scalar_mvs; // The scalar mvs that are multiplying by 
    uint8_t flags;
    int *grades; // If grade preserving, what grades?
    Py_ssize_t g_size; // size of the grade array
    int matrix; // Using a matrix to represent the linear operator
    ProductType left_ptype; // The type of the product for the left mvs
    ProductType right_ptype; // The type of the product for the right mvs
}PyOperatorObject;


#endif
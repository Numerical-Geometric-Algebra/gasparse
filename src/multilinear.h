#ifndef MULTILINEAR_H_
#define MULTILINEAR_H_
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "types.h"

typedef struct Operator{
    PyObject_HEAD
    PyMultivectorObject left_mvs;
    PyMultivectorObject right_mvs;
    int gp; // grade_preserving flag
    int *grades; // If grade preserving, what grades?
    Py_ssize_t g_size; // size of the grade array
    int matrix; // Using a matrix to represent the linear operator
    ProductType left_ptype; // The type of the product for the left mvs
    ProductType right_ptype; // The type of the product for the right mvs 
    int precedence; // right or left precedence, 0 for left, 1 for right
}Operator;


#endif
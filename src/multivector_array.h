#ifndef MULTIVECTOR_ARRAY_H_
#define MULTIVECTOR_ARRAY_H_
#include <python3.11/Python.h>
#include "gasparse.h"

typedef struct PyMultivectorArrayObject{
    PyObject_HEAD
    void *data;
    PyMultivectorMixedMath_Funcs *mixed;
    PyAlgebraObject *GA;
    PyMultivectorSubType *type;
}_PyMultivectorArrayObject;



#endif // MULTIVECTOR_ARRAY_H_
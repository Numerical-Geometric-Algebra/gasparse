#ifndef MULTIVECTOR_TYPES_H_
#define MULTIVECTOR_TYPES_H_
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "types.h"
//#include "common.h"

SparseMultivector init_sparse_empty(Py_ssize_t size);
SparseMultivector alloc_sparse(Py_ssize_t size);
PyMultivectorIter *init_multivector_iter(PyMultivectorObject *data, Py_ssize_t size);
void free_multivector_iter(PyMultivectorIter *iter, Py_ssize_t size);

void sparse_free_(SparseMultivector sparse);
void blades_free_(BladesMultivector blades);
void dense_free_(DenseMultivector dense);
void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size);
SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense, Py_ssize_t size);
int is_bigger_metric(PyAlgebraObject *ga0, PyAlgebraObject *ga1);

int comp_abs(ga_float v, ga_float p);
BladesMultivector init_blades_empty(Py_ssize_t size);
DenseMultivector init_dense_empty(Py_ssize_t size);

#endif
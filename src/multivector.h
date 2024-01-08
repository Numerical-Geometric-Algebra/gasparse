#ifndef MULTIVECTOR_H_
#define MULTIVECTOR_H_
#define PY_SSIZE_T_CLEAN
#include "types.h"

// mixed type operation tables
extern PyMultivectorMixedMath_Funcs multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs cast_multivector_mixed_fn;
extern PyMultivectorMixedMath_Funcs largemultivector_mixed_fn;

// same type operation tables
extern PyMultivectorSubType largemultivector_subtypes_array[3];
extern PyMultivectorSubType multivector_subtypes_array[3];

int comp_abs(ga_float v, ga_float p);
int ga_check_value_precision(PyAlgebraObject *ga, ga_float v);
SparseMultivector init_sparse_empty(Py_ssize_t size);
DenseMultivector init_dense_empty(Py_ssize_t size);
BladesMultivector init_blades_empty(Py_ssize_t size);
void sparse_free_(SparseMultivector sparse);
void blades_free_(BladesMultivector blades);
void dense_free_(DenseMultivector dense);

void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size);
SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense, Py_ssize_t size);


#endif // MULTIVECTOR_H_

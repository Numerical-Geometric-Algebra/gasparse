#ifndef EIN_EINSUM_H_
#define EIN_EINSUM_H_
#include "grade_einsum.h"
#include "parser.h"

typedef struct iterator_t{
    tensor_strides ts;
    void **data;
    size_t *index;
    int *depth;
    size_t sizeof_data;
}iterator_t;



tensor_multivectors append_out_tensor_t(
    subscript_shape,
    char*,
    size_t,
    tensor_multivectors);

void einsum_sum_prods_t(
    tensor_strides,
    tensor_multivectors,
    operator_functions,
    void*);

void sum_of_products_t(
    tensor_multivectors,
    operator_functions,
    void*,
    iterator_t);

iterator_t init_iterator_t(tensor_strides,void**,size_t);
size_t get_nbr_iters_t(iterator_t iter, int depth);

#endif // EIN_EINSUM_H_

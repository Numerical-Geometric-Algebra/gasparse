#ifndef EINSUM_H_
#define EINSUM_H_
#include <stdio.h>
#include <stdlib.h>

#include "parser.h"

typedef struct tensor{
    void *data;
    size_t *shapes;
    size_t shape_size;
    size_t data_size;
}tensor;

typedef struct tensor_multivectors{ // general
    void **data;
    size_t **shapes;
    size_t *shape_size;
    size_t *data_size; // total size of each tensor
    size_t size; // size of the shape size (number of multivector tensors)
    size_t type_size; // sizeof(type)
}tensor_multivectors;


typedef struct labels{
    size_t size;
    char *op_labels;
}labels;

typedef struct symbols{
    size_t *size;
    char **subscripts;
    size_t size_;
}symbols;

typedef struct tensor_strides{
    size_t **strides; // strides for each tensor for each symbol
    size_t *n_strides;
    size_t n_tensors;
    size_t n_symbols;
}tensor_strides;

typedef struct iterator{
    tensor_strides ts;
    void **data;
    size_t *index;
    int *depth;
    size_t sizeof_data;
    int *right;
}iterator;

typedef struct symbol_shape{
    size_t size;
    char *symbols; // list of symbols
    size_t *shape; // shape of each symbol
}symbol_shape;

typedef struct operator_functions{
    void *(*atomic_add)(void *data, size_t size, void *extra);
    void *(*add)(void *a,void *b, void *extra);
    void *(*product)(void *a,void *b, void *extra);
    void (*init)(void *data, size_t size);
    void (*assign)(void *data, void *temp);
    void (*free)(void *data, size_t size);
    /*void (*update_extra)(void *extra,
                         unsigned int *left_grades, size_t left_size,
                         unsigned int *right_grades, size_t right_size,
                         unsigned int *out_grades, size_t out_size);*/
}operator_functions;


labels parse_subscripts(char*,size_t,size_t);
symbols parse_args(char*,size_t);
void free_symbols(symbols);
symbols parse_all(char*,size_t,size_t*,size_t);
void free_symbol_shape(symbol_shape);

size_t get_nbr_inner_iters(iterator);
iterator init_iterator(tensor_strides,void**,size_t);

int general_iterator(iterator,int);
void free_tensor_strides(tensor_strides);

symbol_shape get_all_symbols(symbols,size_t**,size_t*,size_t);
tensor_strides compute_strides(size_t**,symbols,symbol_shape);


int main_einsum(
    tensor_multivectors,
    void*,
    operator_functions,
    symbols,
    tensor*);

tensor_multivectors append_out_tensor(
    symbol_shape,
    char*,
    size_t,
    tensor_multivectors);

void einsum_sum_prods(
    tensor_strides,
    tensor_multivectors,
    operator_functions,
    void*);

void sum_of_products(
    tensor_multivectors,
    operator_functions,
    void*,
    iterator);

#endif // EINSUM_H_

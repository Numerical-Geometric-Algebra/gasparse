#ifndef EINSUM_H_
#define EINSUM_H_

#include "dense.h"
#include "sparse.h"
#include "grade_sparse.h"

typedef struct dense_tensor_multivector{
    dense **data;
    size_t **shapes;
    size_t *shape_size;
    size_t size; // size of the shape_size
    map m;
    dense_grade_map dgm;
}dense_tensor_multivector;

typedef struct sparse_tensor_multivectors{
    sparse **data;
    size_t **shapes;
    size_t *shape_size;
    size_t size; // size of the shape_size
    map m;
    float precision;
    dense_grade_map dgm;
}sparse_tensor_multivectors;

typedef struct graded_tensor_multivectors{
    blades **data;
    size_t **shapes;
    size_t *shape_size;
    size_t *data_size; // total size of each tensor
    size_t size; // size of the shape size (number of multivector tensors)
    map m;
    grade_map gm;
    float precision;
}graded_tensor_multivectors;

typedef struct labels{
    size_t size;
    char *op_labels;
}labels;

typedef struct symbols{
    size_t *size;
    char **subscripts;
    size_t size_;
}symbols;

typedef struct graded_tensor{
    blades *data;
    size_t *shapes;
    size_t shape_size;
    size_t data_size;
}graded_tensor;

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
    size_t sizeof_data;
}iterator;

typedef struct symbol_shape{
    size_t size;
    char *symbols; // list of symbols
    size_t *shape; // shape of each symbol
}symbol_shape;

labels parse_subscripts(char*,size_t,size_t);
symbols parse_args(char*,size_t);
void free_symbols(symbols);
symbols parse_all(char*,size_t,size_t*,size_t);
void free_symbol_shape(symbol_shape);

void free_tensors_holder(graded_tensor_multivectors);
size_t get_nbr_inner_iters(iterator);
iterator init_iterator(tensor_strides,void**,size_t);
int outter_iterator(iterator);
int inner_iterator(iterator);
void einsum_sum_prods(tensor_strides,graded_tensor_multivectors);
void einsum_no_sum_prods(tensor_strides,graded_tensor_multivectors);
void sum_of_products(graded_tensor_multivectors,iterator);
graded_tensor vector_matrix_mult(graded_tensor_multivectors);

void free_tensor_strides(tensor_strides);

int main_einsum(graded_tensor_multivectors,symbols,graded_tensor*);
symbol_shape get_all_symbols(symbols,size_t**,size_t*,size_t);
graded_tensor_multivectors append_out_tensor(symbol_shape,char*,size_t,graded_tensor_multivectors);
tensor_strides compute_strides(size_t**,symbols,symbol_shape);

#endif // EINSUM_H_

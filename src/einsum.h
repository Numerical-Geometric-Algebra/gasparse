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


typedef struct general_tensor{
    void *data;
    size_t *shapes;
    size_t shape_size;
    size_t data_size;
}general_tensor;

typedef struct general_tensor_multivectors{ // general
    void **data;
    size_t **shapes;
    size_t *shape_size;
    size_t *data_size; // total size of each tensor
    size_t size; // size of the shape size (number of multivector tensors)
    size_t type_size; // sizeof(type)
}general_tensor_multivectors;

typedef struct general_extra{
    void *extra;
    size_t size;
    size_t type_size;
}general_extra;

typedef struct operator_functions{
    void *(*atomic_add)(void *data, size_t size, void *extra);
    void *(*add)(void *a,void *b, void *extra);
    void *(*product)(void *a,void *b,void *extra);
    void (*init)(void *data, size_t size);
    void (*assign)(void *data, void *temp);
    void (*free)(void *data, size_t size);
}operator_functions;


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


int general_main_einsum(general_tensor_multivectors,general_extra,operator_functions,symbols,general_tensor*);
general_tensor_multivectors general_append_out_tensor(symbol_shape,char*,size_t,general_tensor_multivectors);
void general_einsum_sum_prods(
    tensor_strides,
    general_tensor_multivectors,
    operator_functions,
    general_extra);
void general_sum_of_products(
    general_tensor_multivectors,
    operator_functions,
    general_extra,
    iterator);


#endif // EINSUM_H_

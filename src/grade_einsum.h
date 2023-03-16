#ifndef GRADE_EINSUM_H_
#define GRADE_EINSUM_H_

#include "parser.h"

enum order{LEFT, RIGHT, OUT, NONE=-1};
typedef enum order order;

typedef struct tensor_strides{
    size_t **strides; // n_tensors x n_symbols
    size_t *n_strides;
    size_t n_tensors;
    size_t n_symbols;
}tensor_strides;

typedef struct grades_strides{
    size_t **strides; // n_operations x n_symbols
    size_t *n_strides;
    size_t n_operations;
    size_t n_symbols;
}grades_strides;


typedef struct iterator{
    tensor_strides ts;
    grades_strides gs;
    void **data; // n_tensors
    int **grades; //  grade subscript by symbol : n_operations x n_symbols
    size_t *index; // index by symbol
    int *depth;
    size_t sizeof_data;
}iterator;


typedef struct subscript_shape{
    size_t size;
    char *subscripts; // all subscripts including grade subscripts
    size_t *shape; // for grade it will be max grade
}subscript_shape;

/*typedef struct symbol_shape{
    size_t size;
    char *symbols; // list of symbols
    char *subscripts; // list of subscripts
    char *grades; // list of grade subscripts
    int *operation_order; // left or right multiplication
    order *grade_order; // indicates left, right, output grade or none : n_symbols
    size_t *sub_shape; // shape of each subscript
    size_t max_grade; // shape of each grade
    size_t n_symbols; // number of symbols
}symbol_shape;*/

typedef struct subscripts_info{
    size_t sub_size_;
    size_t *sub_size;
    char **subscripts;
    size_t grade_size_;
    size_t *grade_size;
    char **grades;
}subscripts_info;



typedef struct strides_info{
    order *grade_order; // indicates left, right, output grade or none : n_symbols
    char *operator; // indicates the operation : n_operations
    int *operation_order; // left or right multiplication
    tensor_strides ts;
    grades_strides gs;
}strides_info;


typedef struct operator_functions{
    void *(*atomic_add)(void *data, size_t size, void *extra);
    void *(*add)(void *a,void *b, void *extra);
    void *(*product)(void *a,void *b, void *extra);
    void (*init)(void *data, size_t size);
    void (*assign)(void *data, void *temp);
    void (*free)(void *data, size_t size);
    void (*update_extra)(void *extra, int *grades, int *grade_order, size_t size, char operator);
}operator_functions;

strides_info compute_strides_info(expression_struct es);
#endif // GRADE_EINSUM_H_

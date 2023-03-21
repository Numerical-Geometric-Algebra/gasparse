#ifndef GRADE_EINSUM_H_
#define GRADE_EINSUM_H_

#include "parser.h"

#define MAX_SYMBOLS 8
#define MAX_SUBSCRIPTS 12

enum order{LEFT, RIGHT, OUT, NONE=-1};
typedef enum order order;

typedef struct tensor{
    void *data;
    size_t *shapes;
    size_t shape_size;
    size_t data_size;
}tensor;

typedef struct tensor_strides{
    size_t **strides; // n_tensors x n_subscripts
    size_t *n_strides;
    size_t n_tensors;
    size_t n_subscripts;
}tensor_strides;

typedef struct grades_strides{
    size_t **strides; // n_grades x n_subscripts
    size_t max_grade;
    size_t n_grades;
    size_t n_subscripts;
}grades_strides;


typedef struct iterator{
    tensor_strides ts;
    grades_strides gs;
    void ***data; // n_tensors
    size_t **grades; // grade subscript by symbol : n_grades
    size_t *index; // index by symbol
    int *depth;
    size_t out_index;
    size_t sizeof_data;
}iterator;

typedef struct tensor_multivectors{ // general
    void **data;
    size_t **shapes;
    size_t *shape_size;
    size_t *data_size; // total size of each tensor
    size_t size; // size of the shape size (number of multivector tensors)
    size_t type_size; // sizeof(type)
}tensor_multivectors;

typedef struct subscript_struct{
    size_t *size;
    char **subscripts;
    size_t size_;
}subscript_struct;


typedef struct subscript_shape{
    size_t size;
    char *subscripts;
    size_t *shape;
}subscript_shape;

typedef struct grades_struct{
    size_t **left,**right,**out;
    size_t left_size,right_size,out_size;
    char operator; // this doesn't change through iterations
}grades_struct;

typedef struct operation_tree{
    void **left; // is a pointer to a tensor
    void **right; // is a pointer to a tensor
    struct operation_tree *right_op;
    struct operation_tree *left_op;
    struct operation_tree *up_op;
    grades_struct grades;
    void *result; // the result is stored here
    int visited;
}operation_tree;

typedef struct symbol_iterator{
    size_t *repeated;
    size_t *pos; // position of the symbol
    char *symbols; // repeated symbols
    size_t size;
}symbol_iterator;

typedef struct operator_functions{
    void *(*atomic_add)(void *data, size_t size, void *extra);
    void *(*add)(void *a,void *b, void *extra);
    void *(*product)(void *a, void *b, void *extra, grades_struct grades);
    void (*init)(void *data, size_t size);
    void (*assign)(void *data, void *temp);
    void (*free)(void *data, size_t size);
}operator_functions;




void free_subscript_struct(subscript_struct sub);
void free_subscript_shape(subscript_shape sp);
void free_symbol_iterator(symbol_iterator iter);

void initialize_operation_tree(operation_tree *tree);

void einsum_sum_prods(
    operator_functions opfs,
    void* extra, iterator iter, operation_tree *op_tree);

void *recursive_products(operator_functions opfs,
                         void *extra,
                         operation_tree *op_tree,int *free_result);

void initialize_iterator(iterator *iter);

int set_grades(char *sub,
               size_t ***grades,
               size_t **grade_strides,
               char **grade_subscripts,
               size_t *size);


char *get_grade_subscripts(char **subscripts, size_t size);

tensor_multivectors append_out_tensor(tensor_multivectors tmvs,
                                      subscript_shape *sps,
                                      char *output_subscripts,
                                      char *grade_subscript,
                                      size_t max_grade, tensor *out);

int iterate_expression_operations(operation_tree **ot,
                                  expression_struct **es,
                                  int visited);

void reset_symbol_iterator(symbol_iterator iter);
int next_symbol(symbol_iterator iter, char symbol);


subscript_shape get_all_subscripts(subscript_struct *sub,
                                   size_t **shapes,
                                   size_t *shape_size,
                                   size_t size);

tensor_strides compute_tensor_strides(size_t **shapes,
                               subscript_struct sub,
                               subscript_shape sp);

int main_einsum(
    tensor_multivectors tmvs,
    void* extra,
    operator_functions opfs,
    expression_struct *es,
    char *output_subscripts,
    size_t max_grade, tensor* out);

#endif // GRADE_EINSUM_H_

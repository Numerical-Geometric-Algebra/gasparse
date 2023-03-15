#ifndef PARSER_H_
#define PARSER_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// for fixed grades I can set strides to 0 and shape of the corresponding index to 0
// these strides are for the output tensor only
// if stride is 0 but shape is not iterate in the inner iterator otherwise in the outer iterator

/*
typedef struct grade_index{
    unsigned int *l;
    unsigned int *r;
    unsigned int *k;
    size_t *l_strides;
    size_t *r_strides;
    size_t *k_strides;
    size_t *l_shape; // for each subscripts of the left multivector how many strides
    size_t *r_shape; // for each subscripts of the right multivector how many strides
    size_t *k_shape; // for each subscripts of the output multivector how many strides
    size_t l_size; // left multivector
    size_t r_size; // right multivector
    size_t k_size; // output multivector
}grade_index;

typedef struct operator_iterator{
    char *identifier;    // identifies the product to be computed
    char *prod_order;    // specifies if it is left or right multiplication
    grade_index *index;  // index for each product for each subscript
    char *input;         // indentifies what is the input in each product
    size_t prod_index;   // selects the identifier
    size_t n_products;   // size of the identifier
    size_t input_size;
}operator_iterator;
*/
/*
typedef struct operation_map{
    unsigned int *left_grades;  size_t left_size;
    unsigned int *right_grades; size_t right_size;
    unsigned int *out_grades;   size_t out_size;
    char operation; // defines the signs array to be choosen
}operation_map;*/

/* typedef struct operand{ */
/*     void *data; // this can point to a temporary variable */
/*     unsigned int *grades; */
/*     size_t grades_size; */
/* }operand; */

/* typedef struct operation{ */
/*     operand *left; */
/*     operand *right; */
/*     struct operation *up; */
/*     unsigned int *grades; size_t grades_size; */
/*     char operator; */
/* }operation; */

typedef struct sub_expression{
    char *grades; // grades per grade selection
    char symbol;
    char *subscripts;
}sub_expression;

typedef struct expression_struct{
    sub_expression left_sub;
    sub_expression right_sub;
    struct expression_struct *left;
    struct expression_struct *right;
    char *grades;
    char operator;
}expression_struct;

int parse_expression(char*,size_t);

void init_subexpression(sub_expression* sub);
void init_expression_struct(expression_struct *es);
int find_matching_angle_brackets(char *expression, size_t size, int *beg);
int sub_expression_parser(char *expression, size_t size, sub_expression *m);
int parse_expression_struct(char *expression, size_t size, sub_expression *right_sub,  sub_expression *left_sub, char *operator);
char *get_subscripts(char *expression, int beg, int end);
int recursive_parser(char *expression, size_t size, expression_struct *es);
#endif // PARSER_H_

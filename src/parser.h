#ifndef PARSER_H_
#define PARSER_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_SUBSCRIPTS 12

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
    struct expression_struct *up;
    int visited;
    char *grades;
    char operator;
}expression_struct;

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


int parse_simple_expression(char*,size_t,subscript_struct*);
expression_struct *parse_expression(char*,size_t,char**);
subscript_shape get_all_subscripts(subscript_struct*,size_t**,size_t*,size_t);
int separate_expression(char*,size_t,char**,char**);

int is_symbol(char c);
int is_operator(char c);
int is_grade_subscripts(char c);
int is_subscript(char c);
int is_constant_grade(char c);

void init_subexpression(sub_expression* sub);
void init_expression_struct(expression_struct *es);
int find_matching_angle_brackets(char *expression, size_t size, int *beg);
int sub_expression_parser(char *expression, size_t size, sub_expression *m);
int parse_expression_struct(char *expression, size_t size, sub_expression *right_sub,  sub_expression *left_sub, char *operator);
char *get_subscripts(char *expression, int beg, int end);
int recursive_parser(char *expression, size_t size, expression_struct *es);
void free_expression_struct_recursive(expression_struct *es);
#endif // PARSER_H_

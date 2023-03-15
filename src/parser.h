#ifndef PARSER_H_
#define PARSER_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

expression_struct *parse_expression(char*,size_t);

int is_symbol(char c);
void init_subexpression(sub_expression* sub);
void init_expression_struct(expression_struct *es);
int find_matching_angle_brackets(char *expression, size_t size, int *beg);
int sub_expression_parser(char *expression, size_t size, sub_expression *m);
int parse_expression_struct(char *expression, size_t size, sub_expression *right_sub,  sub_expression *left_sub, char *operator);
char *get_subscripts(char *expression, int beg, int end);
int recursive_parser(char *expression, size_t size, expression_struct *es);
#endif // PARSER_H_

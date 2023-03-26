#ifndef CAYLEY_H_
#define CAYLEY_H_

#include <stdlib.h>
#include <stdio.h>

typedef struct map{
    char **sign;
    size_t **bitmap;
    size_t size;
}map;

typedef struct grade_map{ // used to map bitmap to position and grade
    size_t *grade;
    size_t *position;
    size_t *grade_size; // number of basis vectors by grade
    size_t max_grade;
    size_t size;
}grade_map;

// used to encode the information to compute the general product
typedef struct project_map{
    size_t *left;
    size_t left_size; // left multivector
    size_t *right;
    size_t right_size; // right multivector
    size_t *out;
    size_t out_size; // output multivector
}project_map;


typedef struct dense_grade_map{
    size_t max_grade;
    size_t grade_size;
    size_t *grade; // grade of each basis blade
}dense_grade_map;


int comp_abs_eq(int,int);

map cayley_table(size_t,size_t,size_t);
map inner_cayley_table(map,grade_map);
map outer_cayley_table(map,grade_map);
void sub_algebra(size_t,char**,int);
void free_map(map);
map init_map(size_t);
size_t grade(size_t);
grade_map bitmap_grade_map(size_t);
void free_grade_map(grade_map);


size_t* get_grade_bool(size_t *,size_t,size_t);
map invert_map(map);

#endif // CAYLEY_H_

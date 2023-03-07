#ifndef CAYLEY_H_
#define CAYLEY_H_

#include <stdlib.h>
#include <stdio.h>

typedef struct map{
    int **sign;
    unsigned int **bitmap;
    size_t size;
}map;

typedef struct grade_map{ // used to map bitmap to position and grade
    unsigned int *grade;
    unsigned int *position;
    unsigned int *grade_size; // number of basis vectors by grade
    size_t max_grade;
    size_t size;
}grade_map;

// used to encode the information to compute the general product
typedef struct project_map{
    unsigned int *l;
    size_t l_size; // left multivector
    unsigned int *r;
    size_t r_size; // right multivector
    unsigned int *k;
    size_t k_size; // output multivector
}project_map;


typedef struct dense_grade_map{
    size_t max_grade;
    size_t grade_size;
    unsigned int *grade; // grade of each basis blade
}dense_grade_map;


int comp_abs_eq(size_t,size_t);

map cayley_table(size_t,size_t,size_t);
map inner_cayley_table(map,grade_map);
map outer_cayley_table(map,grade_map);
void sub_algebra(unsigned int,int**,int);
void free_map(map);
map init_map(size_t);
unsigned int grade(unsigned int);
grade_map bitmap_grade_map(size_t);
void free_grade_map(grade_map);


unsigned int* get_grade_bool(unsigned int *,size_t,size_t);
map invert_map(map);

#endif // CAYLEY_H_

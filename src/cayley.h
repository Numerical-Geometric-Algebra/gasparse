#ifndef CAYLEY_H_
#define CAYLEY_H_

#include <stdlib.h>
#include <stdio.h>

typedef struct map{
    int **sign;
    unsigned int **bitmap;
    size_t size;
}map;

typedef struct grade_map{ // used to map bitmap to position in the grade and grade
    unsigned int *grade;
    unsigned int *position;
    unsigned int *grade_size; // number of basis vectors by grade
    size_t max_grade;
    size_t size;
}grade_map;


map cayley_table(size_t,size_t,size_t);
void sub_algebra(unsigned int,int**,int);
void free_map(map);
void print_map(map);
unsigned int grade(unsigned int);
grade_map bitmap_grade_map(size_t);
void free_grade_map(grade_map);

#endif // CAYLEY_H_

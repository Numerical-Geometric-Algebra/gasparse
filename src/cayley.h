#ifndef CAYLEY_H_
#define CAYLEY_H_

#include <stdlib.h>
#include <stdio.h>

typedef struct map{
    int **sign;
    unsigned int **bitmap;
    size_t size;
}map;

typedef struct graded_map{
    int **sign;
    unsigned int **bitmap;
    int **grade;
    size_t size;
}graded_map;

map cayley_table(size_t,size_t,size_t);
void sub_algebra(unsigned int,int**,int);
void free_map(map);
void print_map(map);
unsigned int grade(unsigned int v);

#endif // CAYLEY_H_

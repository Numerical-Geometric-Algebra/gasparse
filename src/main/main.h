#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ALGEBRA_SIZE 10

typedef struct map{
    int **sign;
    int **bitmap;
    int **grade;
    size_t size;
}map;


// stores mixed grade bitmaps
typedef struct mixed_bitmap{
    unsigned int p; // identity-square
    unsigned int q; // negative-square
    unsigned int r; // zero-square
}mixed_bitmap;

typedef struct algebra_map{
    signed int sign:2; // This value can only be set to -1,0,1
    mixed_bitmap index;
    unsigned int grade;
}algebra_map;

typedef struct single_component{
    unsigned int index:MAX_ALGEBRA_SIZE;
    double value;
    unsigned int grade;
}single_component;

typedef struct components{
    single_component **sing_comp;
    unsigned int size;
}components;

typedef struct graded_components{
    components comps;
    unsigned int grade;
}graded_components;

typedef struct graded_multivector{
    graded_components **graded_comps;
    unsigned int size;
}graded_multivector;

typedef struct grade_index{
    unsigned int max_index:MAX_ALGEBRA_SIZE;
    unsigned int min_index:MAX_ALGEBRA_SIZE;
}grade_index;

typedef struct sparse_multivector{
    components comps;
    graded_multivector graded_mvector;
}sparse_multivector;

typedef struct multivectors{
    sparse_multivector *mvs;
    algebra_map **map; // this will be of size algebra_size by algebra_size
    unsigned int algebra_size:MAX_ALGEBRA_SIZE;
    unsigned int grades_size;
    unsigned int mvs_size;
    grade_index *grade_idx;
    mixed_bitmap metric;
}multivectors;

single_component **initialize_comps(int);
single_component **populate_comps(single_component**,int,int);
sparse_multivector grade_selection(sparse_multivector,unsigned int);
graded_multivector initialize_graded_mv(unsigned int);
graded_multivector populate_graded_mv(graded_multivector,single_component **,unsigned int,unsigned int);

// Algebra maping generation definitions
unsigned int compute_grade(unsigned int x);
int sign_reorder(mixed_bitmap,mixed_bitmap);
algebra_map geo_prod(mixed_bitmap,mixed_bitmap);
algebra_map **cayley_table(size_t,size_t,size_t,size_t*);
void free_map(algebra_map**,size_t);
void print_map(algebra_map**,size_t);

map cayley_table_map(size_t,size_t,size_t);
void sub_algebra_cayley(unsigned int,int**,int);
void free_map_map(map);
void print_map_map(map);

#endif // MAIN_H_

#ifndef SPARSE_H_
#define SPARSE_H_
#include "cayley.h"

typedef struct sparse{
    int *bitmap;
    float *value;
    unsigned int size;
}sparse;

typedef struct sparse_multivectors{
    sparse a;
    sparse b;
    map m;
    float precision;
    dense_grade_map dgm;
}sparse_multivectors;


sparse initialize_sparse(unsigned int);
int comp_abs(float,float);
sparse sparse_product(sparse_multivectors);
sparse sparse_product_(sparse,sparse,map,float);
sparse sparse_grade_project(sparse,unsigned int*,size_t,dense_grade_map);
sparse sparse_dense_to_sparse_sparse(sparse,unsigned int);
void sparse_remove_small(sparse,float,unsigned int*);
sparse sparse_general_product(sparse_multivectors,project_map);
sparse sparse_general_product_(sparse,sparse,map,project_map,dense_grade_map,float);
sparse sparse_scalar_multiply(float,sparse);
void free_sparse(sparse);
sparse sparse_copy(sparse);

sparse sparse_add_append(sparse,sparse);
sparse sparse_atomic_add_append(sparse**,size_t);
sparse sparse_add_add_(sparse,sparse,unsigned int,float);
sparse sparse_atomic_add_add_(sparse **,size_t,size_t,float);

sparse sparse_add_add(sparse_multivectors);
#endif // SPARSE_H_

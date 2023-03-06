#ifndef SPARSE_H_
#define SPARSE_H_
#include "cayley.h"

typedef struct sparse{
    int *bitmap;
    float *value;
    unsigned int size;
}sparse;

typedef struct sparse_multivectors{
    sparse *data;
    map m;
    size_t size;
    float precision;
}sparse_multivectors;




sparse initialize_sparse(unsigned int);
int comp_abs(float,float);
sparse sparse_product(sparse_multivectors);
sparse sparse_grade_project(sparse,unsigned int *,size_t,dense_grade_map);
sparse sparse_dense_to_sparse_sparse(sparse,unsigned int);
void sparse_remove_small(sparse, float,unsigned int *);
sparse sparse_general_product(sparse_multivectors,project_map,dense_grade_map);\

void free_sparse(sparse);
sparse sparse_copy(sparse);

#endif // SPARSE_H_

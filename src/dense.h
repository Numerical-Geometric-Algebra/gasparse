#ifndef DENSE_H_
#define DENSE_H_
#include "cayley.h"
#include "sparse.h"

typedef struct dense{
    float *value;
    unsigned int size;
}dense;

typedef struct dense_multivectors{
    dense *data;
    map m;
    size_t size;
}dense_multivectors;

dense initialize_dense(unsigned int);
dense dense_product(dense_multivectors);
dense inverse_dense_product(dense_multivectors);
dense sparse_to_dense(sparse,unsigned int);

dense dense_grade_project(dense,unsigned int *,size_t,dense_grade_map);
dense dense_general_product(dense_multivectors,project_map,dense_grade_map);
#endif // DENSE_H_

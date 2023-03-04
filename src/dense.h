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
map invert_map(map);
dense dense_grade_project(dense,unsigned int *,unsigned int *,size_t,size_t);

#endif // DENSE_H_

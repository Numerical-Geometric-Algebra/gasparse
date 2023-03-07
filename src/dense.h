#ifndef DENSE_H_
#define DENSE_H_
#include "cayley.h"
#include "sparse.h"

typedef struct dense{
    float *value;
    unsigned int size;
}dense;

typedef struct dense_multivectors{
    dense a;
    dense b;
    map m;
    dense_grade_map dgm;
}dense_multivectors;

typedef struct dense_tensor_multivector{
    dense **data;
    size_t **shapes;
    size_t *shape_size;
    map m;
    dense_grade_map dgm;
}dense_tensor_multivector;



dense initialize_dense(unsigned int);
dense dense_product(dense_multivectors);
dense dense_scalar_multiply(float,dense);
dense dense_product_(dense,dense,map);
dense inverse_dense_product(dense_multivectors);
dense sparse_to_dense(sparse,unsigned int);
dense dense_grade_project(dense,unsigned int *,size_t,dense_grade_map);
dense dense_general_product(dense_multivectors,project_map);
dense dense_general_product_(dense,dense,map,project_map,dense_grade_map);

dense dense_add(dense,dense);
dense dense_atomic_add(dense *,size_t);

#endif // DENSE_H_

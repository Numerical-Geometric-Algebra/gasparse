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
sparse geometric_product(sparse_multivectors);

#endif // SPARSE_H_

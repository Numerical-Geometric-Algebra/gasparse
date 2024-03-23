#include <Python.h>
#include "types.h"
#include "multilinear.h"
#include "multivector_types.h"

static int ternary_sparse_product(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse2 = (SparseMultivector*)data2;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign = m.sign[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0->value[i]*sparse1->value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    *sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return 1;
}


// Takes the sandwich product of data0 with data1: mult(data1,data0,data1)
static int ternary_sparse_sandwich_product(void *out, void *data0, void *data1, PyAlgebraObject *ga, ProductType ptype){
    
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap0, bitmap1;
    int sign0, sign1;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign0 = m.sign[sparse1->bitmap[j]][sparse0->bitmap[i]];
            if(!sign0) continue;
            bitmap0 = m.bitmap[sparse1->bitmap[j]][sparse0->bitmap[i]];
            for(Py_ssize_t k = 0; k < sparse1->size; k++){
                sign1 = m.sign[bitmap0][sparse1->bitmap[j]];
                if(!sign1) continue;
                bitmap1 = m.bitmap[bitmap0][sparse1->bitmap[j]];
                if(dense.bitmap[bitmap1] == -1) dense.bitmap[bitmap1] = bitmap1, size++;
                dense.value[bitmap1] += sparse0->value[i]*sparse1->value[j]*sign1;
            }
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    *sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return 1;
}
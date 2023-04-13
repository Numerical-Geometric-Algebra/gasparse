#include "multivector.c"

static SparseMultivector *sparse_array_init_(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize){
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(ssize*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < ssize; i++){
        sparse_array[i] = sparse_init_(bitmap[i],value[i],size[i]);
    }
    return sparse_array;
}

void *sparse_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize, GradeMap gm, Py_ssize_t algebra_size){
    return (void*)sparse_array_init_(bitmap,value,size,ssize);
}

static BladesMultivector *blades_array_init_(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize, GradeMap gm){
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(ssize*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < ssize; i++){
        blades_array[i] = blades_init_(bitmap[i],value[i],size[i],gm);
    }
    return blades_array;
}

void *blades_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize, GradeMap gm, Py_ssize_t algebra_size){
    return (void*)blades_array_init_(bitmap,value,size,ssize,gm);
}

static DenseMultivector *dense_array_init_(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize, Py_ssize_t algebra_size){
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(ssize*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < ssize; i++){
        dense_array[i] = dense_init_(bitmap[i],value[i],size[i],algebra_size);
    }
    return dense_array;
}

void *dense_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize, GradeMap gm, Py_ssize_t algebra_size){
    return (void*)dense_array_init_(bitmap,value,size,ssize,algebra_size);
}

static void sparse_array_free_(SparseMultivector *sparse_array, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        sparse_free_(sparse_array[i]);
    }
}

void sparse_array_free(void *sparse_array, Py_ssize_t size){
    sparse_array_free_((SparseMultivector*)sparse_array,size);
}

static void blades_array_free_(BladesMultivector *blades_array, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        blades_free_(blades_array[i]);
    }
}

void blades_array_free(void *blades_array, Py_ssize_t size){
    blades_array_free_((BladesMultivector*)blades_array,size);
}

static void dense_array_free_(DenseMultivector *dense_array, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        dense_free_(dense_array[i]);
    }
}

void dense_array_free(void *dense_array, Py_ssize_t size){
    dense_array_free_((DenseMultivector*)dense_array,size);
}

static SparseMultivector *sparse_sparse_array_add_(SparseMultivector *left_array, SparseMultivector *right_array, Py_ssize_t size, PyGeometricAlgebraObject ga, int sign){
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < size; i++){
        sparse_array[i] = sparse_sparse_add_(left_array[i],right_array[i],ga,sign);
    }
    return sparse_array;
}

static BladesMultivector *blades_blades_array_add_(BladesMultivector *left_array, BladesMultivector *right_array, Py_ssize_t size, PyGeometricAlgebraObject ga, int sign){
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++){
        blades_array[i] = blades_blades_add_(left_array[i],right_array[i],ga,sign);
    }
    return blades_array;
}

static DenseMultivector *dense_dense_array_add_(DenseMultivector *left_array, DenseMultivector *right_array, Py_ssize_t size, PyGeometricAlgebraObject ga, int sign){
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < size; i++){
        dense_array[i] = dense_dense_add_(left_array[i],right_array[i],ga,sign);
    }
    return dense_array;
}

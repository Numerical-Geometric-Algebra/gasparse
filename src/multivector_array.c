#include "multivector.c"

static PyMultivectorArrayObject *init_multivector_array_object(PyMultivectorArrayObject *src){
    PyMultivectorArrayObject *out = (PyMultivectorArrayObject*)PyMem_RawMalloc(sizeof(PyMultivectorArrayObject));
    out->shapes = (Py_ssize_t*)PyMem_RawMalloc(src->shape_size*sizeof(Py_ssize_t));
    out->strides = (Py_ssize_t*)PyMem_RawMalloc(src->shape_size*sizeof(Py_ssize_t));
    out->GA = src->GA;
    out->math_funcs = src->math_funcs;
    out->data_funcs = src->data_funcs;
    out->type = src->type;
    out->data_size = src->data_size;
    out->shape_size = src->shape_size;
    for(Py_ssize_t i = 0; i < src->shape_size; i++){
        out->shapes[i] = src->shapes[i];
        out->strides[i] = src->strides[i];
    }
    Py_SET_REFCNT((void*)out,1);
    Py_XINCREF((void*)src->GA);
    return out;
}

static void free_multivector_array_object(PyMultivectorArrayObject *self){
    Py_XDECREF((void*)self->GA);
    PyMem_RawFree(self->shapes);
    PyMem_RawFree(self->strides);
    PyMem_RawFree(self);
}





// templated for binary functions do not edit source code




static PyMultivectorArrayObject *unary_sparse_array_add(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = unary_sparse_add_(left_array[i], ga, sign);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_sparse_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    SparseMultivector *data1_array = (SparseMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = binary_sparse_add_(left_array[i], right_array[i], ga, sign);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_sparse_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    SparseMultivector *data1_array = (SparseMultivector*)data1->data;
    SparseMultivector *data2_array = (SparseMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = ternary_sparse_add_(left_array[i], ga, sign);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_sparse_array_product(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = unary_sparse_product_(left_array[i], ga, ptype);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_sparse_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    SparseMultivector *data1_array = (SparseMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = binary_sparse_product_(left_array[i], right_array[i], ga, ptype);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_sparse_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    SparseMultivector *data1_array = (SparseMultivector*)data1->data;
    SparseMultivector *data2_array = (SparseMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = ternary_sparse_product_(left_array[i], ga, ptype);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_sparse_array_reverse(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = unary_sparse_reverse_(left_array[i], ga);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_sparse_array_grade_project(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int *grades, Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *data0_array = (SparseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        sparse_array[i] = unary_sparse_grade_project_(left_array[i], ga, grades, size);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)sparse_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
    free_multivector_array_object(out_array);
    return NULL;
}
static void sparse_array_free(void *array, Py_ssize_t size){
    SparseMultivector *sparse_array = (SparseMultivector*)array;
    for(Py_ssize_t i = 0; i < size; i++)
        sparse_free_(sparse_array[i]);
}

static void *sparse_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize){
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(ssize*sizeof(SparseMultivector));
    Py_ssize_t i = 0;
    for(i = 0; i < ssize; i++){
        sparse_array[i] = sparse_init_(bitmap[i],value[i],size[i]);
        if(sparse_array[i].size == -1)
            goto failure;
    }
    return (void*)sparse_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        sparse_free_(sparse_array[i]);
    PyMem_RawFree(sparse_array);
}


static PyMultivectorArrayObject *unary_blades_array_add(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = unary_blades_add_(left_array[i], ga, sign);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_blades_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    BladesMultivector *data1_array = (BladesMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = binary_blades_add_(left_array[i], right_array[i], ga, sign);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_blades_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    BladesMultivector *data1_array = (BladesMultivector*)data1->data;
    BladesMultivector *data2_array = (BladesMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = ternary_blades_add_(left_array[i], ga, sign);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_blades_array_product(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = unary_blades_product_(left_array[i], ga, ptype);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_blades_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    BladesMultivector *data1_array = (BladesMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = binary_blades_product_(left_array[i], right_array[i], ga, ptype);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_blades_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    BladesMultivector *data1_array = (BladesMultivector*)data1->data;
    BladesMultivector *data2_array = (BladesMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = ternary_blades_product_(left_array[i], ga, ptype);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_blades_array_reverse(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = unary_blades_reverse_(left_array[i], ga);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_blades_array_grade_project(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int *grades, Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *data0_array = (BladesMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        blades_array[i] = unary_blades_grade_project_(left_array[i], ga, grades, size);
        if(blades_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)blades_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
    free_multivector_array_object(out_array);
    return NULL;
}
static void blades_array_free(void *array, Py_ssize_t size){
    BladesMultivector *blades_array = (BladesMultivector*)array;
    for(Py_ssize_t i = 0; i < size; i++)
        blades_free_(blades_array[i]);
}

static void *blades_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize){
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(ssize*sizeof(BladesMultivector));
    Py_ssize_t i = 0;
    for(i = 0; i < ssize; i++){
        blades_array[i] = blades_init_(bitmap[i],value[i],size[i]);
        if(blades_array[i].size == -1)
            goto failure;
    }
    return (void*)blades_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        blades_free_(blades_array[i]);
    PyMem_RawFree(blades_array);
}


static PyMultivectorArrayObject *unary_dense_array_add(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = unary_dense_add_(left_array[i], ga, sign);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_dense_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    DenseMultivector *data1_array = (DenseMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = binary_dense_add_(left_array[i], right_array[i], ga, sign);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_dense_array_add(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, int sign){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    DenseMultivector *data1_array = (DenseMultivector*)data1->data;
    DenseMultivector *data2_array = (DenseMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = ternary_dense_add_(left_array[i], ga, sign);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_dense_array_product(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = unary_dense_product_(left_array[i], ga, ptype);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *binary_dense_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    DenseMultivector *data1_array = (DenseMultivector*)data1->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = binary_dense_product_(left_array[i], right_array[i], ga, ptype);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *ternary_dense_array_product(
    PyMultivectorArrayObject *data0,
    PyMultivectorArrayObject *data1,
    PyMultivectorArrayObject *data2,
    Py_ssize_t size, ProductType ptype){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    DenseMultivector *data1_array = (DenseMultivector*)data1->data;
    DenseMultivector *data2_array = (DenseMultivector*)data2->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = ternary_dense_product_(left_array[i], ga, ptype);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_dense_array_reverse(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = unary_dense_reverse_(left_array[i], ga);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}


static PyMultivectorArrayObject *unary_dense_array_grade_project(
    PyMultivectorArrayObject *data0,
    Py_ssize_t size, int *grades, Py_ssize_t size){
    PyGeometricAlgebraObject *ga = left->GA;
    PyMultivectorArrayObject *out_array = init_multivector_array_object(left);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *data0_array = (DenseMultivector*)data0->data;
    Py_ssize_t i = 0;
    for(i = 0; i < size; i++){
        dense_array[i] = unary_dense_grade_project_(left_array[i], ga, grades, size);
        if(dense_array[i].size == -1)
            goto failure;
    }
    out_array->data = (void*)dense_array;
    return out_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
    free_multivector_array_object(out_array);
    return NULL;
}
static void dense_array_free(void *array, Py_ssize_t size){
    DenseMultivector *dense_array = (DenseMultivector*)array;
    for(Py_ssize_t i = 0; i < size; i++)
        dense_free_(dense_array[i]);
}

static void *dense_array_init(int **bitmap, ga_float **value, Py_ssize_t *size, Py_ssize_t ssize){
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(ssize*sizeof(DenseMultivector));
    Py_ssize_t i = 0;
    for(i = 0; i < ssize; i++){
        dense_array[i] = dense_init_(bitmap[i],value[i],size[i]);
        if(dense_array[i].size == -1)
            goto failure;
    }
    return (void*)dense_array;
failure:
    for(Py_ssize_t j = 0; j < i; j++) // free the other multivectors
        dense_free_(dense_array[i]);
    PyMem_RawFree(dense_array);
}

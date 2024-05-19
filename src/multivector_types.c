#include <Python.h>
#include <math.h>
#include "pyport.h"
#include "types.h"
#include "common.h"
#include "multivector_types.h"


// returns true if abs(v) < p
int comp_abs(ga_float v, ga_float p){
    ga_float r = (v < 0) ? -v : v;
    return r < p;
}

// returns true if abs(v) < ga->precision
int ga_check_value_precision(PyAlgebraObject *ga, ga_float v){
    ga_float r = (v < 0) ? -v : v;
    return r < ga->precision;
}

PyMultivectorIter *init_multivector_iter(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorIter *iter = (PyMultivectorIter*)PyMem_RawMalloc(size*sizeof(PyMultivectorIter));
    for(Py_ssize_t i = 0; i < size; i++){
        gaiterinitfunc iter_init = data[i].type->data_funcs->iter_init;
        iter[i] = iter_init(data[i].data,data[i].type);
    }
    return iter;
}

void free_multivector_iter(PyMultivectorIter *iter, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++)
        if(iter[i].index)
            free(iter[i].index);
    free(iter);
}

SparseMultivector alloc_sparse(Py_ssize_t size){
    SparseMultivector sparse = {.size = -1};
    sparse.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    sparse.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    if(!sparse.bitmap || !sparse.value){
       PyMem_RawFree(sparse.bitmap);
       PyMem_RawFree(sparse.value);
       PyErr_SetString(PyExc_MemoryError,"Error allocating memory for a sparse multivector");
       return sparse;
    }
    sparse.size = size;
    return sparse;
}

SparseMultivector init_sparse_empty(Py_ssize_t size){
    SparseMultivector sparse = alloc_sparse(size);
    // Set all values to zero
    for(Py_ssize_t i = 0; i < size; i++){
        sparse.bitmap[i] = -1;
        sparse.value[i] = 0;
    }

    return sparse;
}

DenseMultivector init_dense_empty(Py_ssize_t size){
    DenseMultivector dense = {.size = -1};
    dense.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    if(!dense.value){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for a dense multivector");
        return dense;
    }

    dense.size = size;
    for(Py_ssize_t i = 0; i < size; i++)
        dense.value[i] = 0;

    return dense;
}


BladesMultivector init_blades_empty(Py_ssize_t size){
    BladesMultivector blades = {.size = -1};
    blades.data = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    blades.grade = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    if(!blades.data || !blades.grade){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for a blades multivector");
        return blades;
    }
    blades.size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        blades.data[i].bitmap = NULL;
        blades.data[i].value = NULL;
        blades.grade[i] = i;
    }
    return blades;
}

static Py_ssize_t *init_grade_size(GradeMap gm){
    Py_ssize_t *gsize = (Py_ssize_t*)PyMem_RawMalloc((gm.max_grade+1)*sizeof(Py_ssize_t));
    if(!gsize){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for grade size array");
        return NULL;
    }
    for(Py_ssize_t i = 0; i <= gm.max_grade; i++)
        gsize[i] = 0;
    return gsize;
}

SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense, Py_ssize_t size){
    Py_ssize_t k = 0;
    SparseMultivector sparse = init_sparse_empty(size);
    if(sparse.size == -1)
        return sparse;
    for(Py_ssize_t i = 0; i < dense.size; i++){
        if(k < size){
            if(dense.bitmap[i] != -1){
                sparse.bitmap[k] = dense.bitmap[i];
                sparse.value[k] = dense.value[i];
                k++;
            }
        }
    }
    sparse.size = size;
    return sparse;
}

static Py_ssize_t sparse_remove_rel_small(SparseMultivector x, ga_float percentage){
    ga_float max = 0;
    Py_ssize_t size = 0;
    // Find the maximum and the size of x
    for(Py_ssize_t i = 0; i < x.size; i++){
        if(x.bitmap[i] != -1) size++;
        if(max < fabs(x.value[i]))
            max = fabs(x.value[i]);
    }

    // Compare with the maximum
    for(Py_ssize_t i = 0; i < x.size; i++){
        if(fabs(x.value[i]) < max*percentage && x.bitmap[i] != -1){
            x.bitmap[i] = -1; // remove basis element from the multivector
            size--;
        }
    }

    return size;
}

SparseMultivector sparse_remove_relative_small(SparseMultivector x, ga_float percentage){
    Py_ssize_t size = sparse_remove_rel_small(x,percentage);
    return sparse_dense_to_sparse_sparse(x,size);
}



void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size){
     // Remove if value is too small
    for(Py_ssize_t i = 0; i < y.size; i++){
        // Check if value was set
        if(y.bitmap[i] > -1){
            // Check if value is too small
            if(comp_abs(y.value[i],precision)){
                y.bitmap[i] = -1;
                (*size)--;
            }
        }
    }
}




static BladesMultivector sparse_dense_to_blades_sparse(SparseMultivector dense, GradeMap gm){
    BladesMultivector sparse = {.size = -1};
    Py_ssize_t ssize = 0, grade = -1;
    Py_ssize_t *gsize = init_grade_size(gm);
    Py_ssize_t *gindex = init_grade_size(gm);
    if(!gsize || !gindex){
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        return sparse;
    }
    int bitmap;

    for(Py_ssize_t i = 0; i < dense.size; i++){
        if(dense.bitmap[i] == -1) continue;
        grade = gm.grade[dense.bitmap[i]];
        if(!gsize[grade]) gindex[grade] = ssize++; // first time incrementing
        gsize[grade]++;
    }

    if(!ssize){
        sparse.data = NULL;
        sparse.grade = NULL;
        sparse.size = 0;
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        return sparse;
    }

    sparse.data = (SparseMultivector*)PyMem_RawMalloc(ssize*sizeof(SparseMultivector));
    sparse.grade =  (Py_ssize_t*)PyMem_RawMalloc(ssize*sizeof(Py_ssize_t));
    if(!sparse.data || !sparse.grade){
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        sparse.size = -1;
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory");
        return sparse;
    }
    sparse.size = ssize;

    // initialize each grade
    for(Py_ssize_t i = 0; i <= gm.max_grade; i++){ // iterate over grades
        if(!gsize[i]) continue;
        sparse.data[gindex[i]] = init_sparse_empty(gsize[i]);
        sparse.grade[gindex[i]] = i;
    }

    for(Py_ssize_t i = 0; i < dense.size; i++){
        bitmap = dense.bitmap[i];
        if(bitmap == -1) continue;
        // if(fabs(dense.value[i]) < sparse_max*percentage) continue; // ignore relatively small values

        grade = gm.grade[bitmap]; gsize[grade]--;
        sparse.data[gindex[grade]].bitmap[gsize[grade]] = bitmap;
        sparse.data[gindex[grade]].value[gsize[grade]] = dense.value[i];
    }

    PyMem_RawFree(gsize);
    PyMem_RawFree(gindex);
    return sparse;
}


void sparse_free_(SparseMultivector sparse){
    PyMem_RawFree(sparse.value);
    PyMem_RawFree(sparse.bitmap);
}

void blades_free_(BladesMultivector blades){
    if(blades.data){
        for(Py_ssize_t i = 0; i < blades.size; i++){
            PyMem_RawFree(blades.data[i].bitmap);
            PyMem_RawFree(blades.data[i].value);
        }
        PyMem_RawFree(blades.data);
    }
    PyMem_RawFree(blades.grade);
}

void dense_free_(DenseMultivector dense){
    PyMem_RawFree(dense.value);
}

static int scalar_init(void *out, PyAlgebraObject *GA, int *bitmap, ga_float *value, Py_ssize_t size){
    
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    if(!value) *scalar = 0;
    else *scalar = *value;
    return 1;
}

static int sparse_init(void *out, PyAlgebraObject *GA, int *bitmap, ga_float *value, Py_ssize_t size){
    SparseMultivector *sparse = (SparseMultivector *)out;
    if(!size){
        sparse->size = 0;
        sparse->bitmap = NULL;
        sparse->value = NULL;
        return 1;
    }
    
    sparse->value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    sparse->bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    sparse->size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        sparse->value[i] = value[i];
        sparse->bitmap[i] = bitmap[i];
    }
    return 1;
}

static int blades_init(void *out,PyAlgebraObject *ga,int *bitmap, ga_float *value, Py_ssize_t size){
    GradeMap gm = ga->gm;
    BladesMultivector *blades = out;
    if(!size){
        blades->size = 0;
        blades->grade = NULL;
        blades->data = NULL;
        return 1;
    }
    SparseMultivector ssparse = {.bitmap = bitmap, .value = value, .size = size};
    *blades = sparse_dense_to_blades_sparse(ssparse,gm);
    return 1;
}

static int dense_init(void *out, PyAlgebraObject *ga,int *bitmap, ga_float *value, Py_ssize_t size){
    Py_ssize_t algebra_size = ga->asize;
    DenseMultivector *dense = (DenseMultivector *)out;
    if(!size){
        dense->size = 0;
        dense->value = NULL;
        return 1;
    }
    dense->value = (ga_float*)PyMem_RawMalloc(algebra_size*sizeof(ga_float));
    dense->size = algebra_size;
    // set all values to 0
    for(Py_ssize_t i = 0; i < algebra_size; i++)
        dense->value[i] = 0;
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= algebra_size){
            PyMem_RawFree(dense->value);
            dense->value = NULL;
            dense->size = -1;
            return 0; // raise error
        }
        dense->value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return 1;
}

static int cast_to_sparse(PyMultivectorIter *from, void *to, PyAlgebraObject *GA){
    SparseMultivector *psparse = (SparseMultivector*)to;
    if(!from || !psparse){
        return 0;
    }
    
    SparseMultivector sparse;
    sparse.value = (ga_float*)PyMem_RawMalloc(from->niters*sizeof(ga_float));
    sparse.bitmap = (int*)PyMem_RawMalloc(from->niters*sizeof(int));
    Py_ssize_t i = 0;
    while(from->next(from)){
        sparse.value[i] = from->value;
        sparse.bitmap[i] = from->bitmap;
        i++;
    }
    sparse.size = from->niters;
    *psparse = sparse;
    return 1;
}

static int cast_to_dense(PyMultivectorIter *from, void *to, PyAlgebraObject *GA){
    DenseMultivector *pdense = (DenseMultivector*)to;
    if(!from || !pdense){
        return 0;
    }

    DenseMultivector dense = {.size = GA->asize, .value = NULL};
    dense.value = (ga_float*)PyMem_RawMalloc(dense.size*sizeof(ga_float));
    for(Py_ssize_t i = 0; i < dense.size; i++)
        dense.value[i] = 0;

    while(from->next(from))
        dense.value[from->bitmap] += from->value;

    *pdense = dense;
    return 1;
}

static int cast_to_blades(PyMultivectorIter *from, void *to,PyAlgebraObject *GA){
    BladesMultivector *pblades = (BladesMultivector*)to;
    
    if(!from || !pblades ){
        return 0;
    }

    SparseMultivector sparse = {.size = from->niters, .value = NULL, .bitmap = NULL};
    sparse.value = (ga_float*)PyMem_RawMalloc(from->niters*sizeof(ga_float));
    sparse.bitmap = (int*)PyMem_RawMalloc(from->niters*sizeof(int));
    Py_ssize_t i = 0;
    while(from->next(from)){
        sparse.value[i] = from->value;
        sparse.bitmap[i] = from->bitmap;
        i++;
    }
    *pblades = sparse_dense_to_blades_sparse(sparse,GA->gm);
    sparse_free_(sparse);
    return 1;
}

static int scalar_iter_next(PyMultivectorIter *iter){
    if(*iter->index >= 1){
        *iter->index = 0;
        return 0;
    }
    iter->bitmap = 0;
    iter->value = *((ScalarMultivector*)iter->data);
    (*iter->index)++;
    return 1;
}

static int sparse_iter_next(PyMultivectorIter *iter){
    SparseMultivector *sparse = (SparseMultivector*)iter->data;
    if(*iter->index >= sparse->size){
        *iter->index = 0;
        return 0;
    }
    iter->bitmap = sparse->bitmap[*iter->index];
    iter->value = sparse->value[(*iter->index)++];
    iter->grade = GRADE(iter->bitmap);
    return 1;
}

static int dense_iter_next(PyMultivectorIter *iter){
    DenseMultivector *dense = (DenseMultivector*)iter->data;
    if(*iter->index >= dense->size){
        *iter->index = 0;
        return 0;
    }
    iter->bitmap = *iter->index;
    iter->value = dense->value[(*iter->index)++];
    iter->grade = GRADE(iter->bitmap);
    return 1;
}

static int blades_iter_next(PyMultivectorIter *iter){
    BladesMultivector *blades = (BladesMultivector*)iter->data;
    if(*iter->index >= blades->size){
        iter->index[1] = 0;
        iter->index[0] = 0;
        return 0;
    }
    iter->bitmap = blades->data[*iter->index].bitmap[iter->index[1]];
    iter->value = blades->data[*iter->index].value[iter->index[1]++];
    iter->grade = blades->grade[*iter->index];
    if(iter->index[1] >= blades->data[*iter->index].size){
        iter->index[1] = 0;
        (*iter->index)++;
    }
    return 1;
}

static PyMultivectorIter scalar_iter_init(void *data, PyMultivectorSubType *type){
    PyMultivectorIter iter;
    iter.data = data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = type->ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = 1;
    iter.next = type->data_funcs->iter_next;
    iter.type_name = type->type_name;
    return iter;
}

static PyMultivectorIter sparse_iter_init(void *data, PyMultivectorSubType *type){
    PyMultivectorIter iter;
    SparseMultivector *sparse = (SparseMultivector*)data;
    iter.data= data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = type->ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = sparse->size;
    iter.next = type->data_funcs->iter_next;
    iter.type_name = type->type_name;
    return iter;
}

static PyMultivectorIter dense_iter_init(void *data, PyMultivectorSubType *type){
    PyMultivectorIter iter;
    DenseMultivector *dense = (DenseMultivector*)data;
    iter.data= data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = type->ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = dense->size;
    iter.next = type->data_funcs->iter_next;
    iter.type_name = type->type_name;
    return iter;
}

static PyMultivectorIter blades_iter_init(void *data, PyMultivectorSubType *type){
    PyMultivectorIter iter;
    BladesMultivector *blades = (BladesMultivector*)data;
    iter.data= data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = type->ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(2*sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.index[1] = 0;
    iter.size = 2;
    iter.niters = 0;
    iter.type_name = type->type_name;
    for(Py_ssize_t i = 0; i < blades->size; i++)
        iter.niters += blades->data[i].size;

    iter.next = type->data_funcs->iter_next;
    return iter;
}


char *bitmap_to_string(int bitmap){
    Py_ssize_t size = GRADE((Py_ssize_t)bitmap) + 2;
    char *str = (char*)PyMem_RawMalloc(size*sizeof(char));
    unsigned int x = (unsigned int)bitmap;
    str[0] = 'e';
    Py_ssize_t c = 1;
    while(x){
        str[c++] = (char)__builtin_ctz(x) + '1'; // count the number of trailing zeros
        x &= x - 1; // clear the least significant bit set
    }
    str[c] = '\0';
    return str;
}

/*
static int sparse_copy(void* out, void* input){
    SparseMultivector *sparse = (SparseMultivector *)input;
    SparseMultivector *copy = (SparseMultivector *)out;// = {.size = -1, .bitmap = NULL, .value = NULL};
    if(sparse->size < 0) return 0;
    if(sparse->size == 0){
        copy->size = 0;
        return 1;
    }

    copy->bitmap = (int*)PyMem_RawMalloc(sparse->size*sizeof(int));
    copy->value = (ga_float*)PyMem_RawMalloc(sparse->size*sizeof(ga_float));
    if(!copy->bitmap || !copy->value){
        PyMem_RawFree(copy->bitmap);
        PyMem_RawFree(copy->value);
        return 0;
    }
    copy->size = sparse->size;
    memcpy(copy->bitmap,sparse->bitmap,sparse->size*sizeof(int));
    memcpy(copy->value,sparse->value,sparse->size*sizeof(ga_float));
    return 1;
}
*/

SparseMultivector sparse_new_(Py_ssize_t size){
    SparseMultivector sparse = {.size = -1};
    sparse.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    sparse.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    if(!sparse.bitmap || !sparse.value){
       PyMem_RawFree(sparse.bitmap);
       PyMem_RawFree(sparse.value);
       PyErr_SetString(PyExc_MemoryError,"Error allocating memory for a sparse multivector");
       return sparse;
    }
    sparse.size = size;

    return sparse;
}

static int unary_scalar_scalaradd(void *out, void* data0,PyAlgebraObject *GA,ga_float value, int sign){
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    ScalarMultivector *scalar0 = (ScalarMultivector *)data0;
    
    *scalar = sign*(*scalar0) + value;
    return 1;
}

static int unary_sparse_scalaradd(void *out, void *data0, PyAlgebraObject *GA, ga_float value, int sign){
    SparseMultivector *sparse0 = (SparseMultivector *)data0;
    SparseMultivector *sparse = (SparseMultivector *)out;
    sparse->size = -1;
    Py_ssize_t scalarindex = -1;
    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        if(sparse0->bitmap[i] == 0){
            scalarindex = i;
            break;
        }
    }

    if(scalarindex != -1){ // Already in the multivector 
        *sparse = init_sparse_empty(sparse0->size);
        if(sparse->size == -1)
            return -1;

        for(Py_ssize_t i = 0; i < sparse0->size; i++){
            sparse->value[i] = sign*sparse0->value[i];
            sparse->bitmap[i] = sparse0->bitmap[i];
        }
        sparse->value[scalarindex] += value;
    }
    else{ // A new element
        *sparse = init_sparse_empty(sparse0->size + 1);
        if(sparse->size == -1)
            return -1;
        sparse->value[0] = value;
        sparse->bitmap[0] = 0;

        for(Py_ssize_t i = 0; i < sparse0->size; i++){
            sparse->value[i+1] = sign*sparse0->value[i];
            sparse->bitmap[i+1] = sparse0->bitmap[i];
        }
    }
    return 1;
}

static int unary_scalar_scalarproduct(void *out, void* data0,PyAlgebraObject *GA,ga_float value){
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    ScalarMultivector *scalar0 = (ScalarMultivector *)data0;
    
    *scalar = value*(*scalar0);
    return 1;
}

static int unary_sparse_scalarproduct(void* out, void* data0, PyAlgebraObject *ga, ga_float value){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = out;

     *sparse = alloc_sparse(sparse0->size);
    if(sparse->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        sparse->bitmap[i] = sparse0->bitmap[i];// copy the bitmap
        sparse->value[i] = value*sparse0->value[i]; // multiply by a scalar
    }

    return 1;
}

static int binary_scalar_add(void *out, void* data0, void* data1,PyAlgebraObject *GA, int sign){
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    ScalarMultivector *scalar0 = (ScalarMultivector *)data0;
    ScalarMultivector *scalar1 = (ScalarMultivector *)data1;
    *scalar = *scalar0 + sign*(*scalar1);
    return 1;
}

static int binary_sparse_add(void *out, void *data0, void *data1, PyAlgebraObject *ga, int sign){
    SparseMultivector *sparse0 = (SparseMultivector*)data0; 
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse = (SparseMultivector*)out;
    
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        dense.bitmap[sparse0->bitmap[i]] = sparse0->bitmap[i];
        dense.value[sparse0->bitmap[i]] += sparse0->value[i];
    }

    for(Py_ssize_t i = 0; i < sparse1->size; i++){
        dense.bitmap[sparse1->bitmap[i]] = sparse1->bitmap[i];
        dense.value[sparse1->bitmap[i]] += sign*sparse1->value[i];
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);
    sparse_free_(dense);
    return 1;
}

// Use this only for geometric or outer product between scalars
static int binary_scalar_product(void *out, void* data0, void* data1, PyAlgebraObject *GA, ProductType ptype){
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    ScalarMultivector *scalar0 = (ScalarMultivector *)data0;
    ScalarMultivector *scalar1 = (ScalarMultivector *)data1;
    *scalar = (*scalar0)*(*scalar1);
    return 1;
}

static int binary_sparse_product(void *out, void *data0, void *data1, PyAlgebraObject *ga, ProductType ptype){
    
    SparseMultivector *sparse0 = (SparseMultivector*)data0; 
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign = m.sign[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0->bitmap[i]][sparse1->bitmap[j]];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sparse0->value[i]*sparse1->value[j]*sign;
        }
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);

    sparse_free_(dense);
    return 1;
}

// Use this only for geometric or outer product between scalars
static int ternary_scalar_product(void *out, void* data0, void* data1,void* data2, PyAlgebraObject *GA, ProductType ptype){
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    ScalarMultivector *scalar0 = (ScalarMultivector *)data0;
    ScalarMultivector *scalar1 = (ScalarMultivector *)data1;
    ScalarMultivector *scalar2 = (ScalarMultivector *)data2;
    *scalar = (*scalar0)*(*scalar1)*(*scalar2);
    return 1;
}

/*

static int ternary_sparse_product_nested(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    SparseMultivector *sparse0 = (SparseMultivector*)data0; 
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse2 = (SparseMultivector*)data2;
    SparseMultivector *sparse = (SparseMultivector*)out;

    CliffordMap m = ga->product[ptype];
    SparseMultivector dense = init_sparse_empty(m.size);
    int sign0, sign1;
    Py_ssize_t bitmap0, bitmap1;

    ga_float value;
    // Compute the ternary product in three nested for loops
    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign0 = m.sign[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(!sign0) continue;
            bitmap0 = m.bitmap[sparse0->bitmap[i]][sparse1->bitmap[j]];
            value = sparse0->value[i]*sparse1->value[j]*sign0;
            for(Py_ssize_t k = 0; k < sparse1->size; k++){
                sign1 = m.sign[bitmap0][sparse2->bitmap[k]];
                if(!sign1) continue;
                bitmap1 = m.bitmap[bitmap0][sparse2->bitmap[j]];
                dense.bitmap[bitmap1] = bitmap1;
                dense.value[bitmap1] = sign1*value*sparse2->value[j];
            }
        }
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);
    if(sparse->size == -1){
        sparse_free_(dense);
        return 0;
    }

    sparse_free_(dense);
    return 1;
}
*/

static int ternary_sparse_product(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    SparseMultivector *sparse0 = (SparseMultivector*)data0; 
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse2 = (SparseMultivector*)data2;
    SparseMultivector *sparse = (SparseMultivector*)out;
    
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return 0;

    Py_ssize_t bitmap;
    Py_ssize_t j;
    int sign;

    // dense0 <- geometric_product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign = m.sign[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0->bitmap[i]][sparse1->bitmap[j]];
            dense0.bitmap[bitmap] = bitmap;
            dense0.value[bitmap] += sparse0->value[i]*sparse1->value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(m.size);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return 0;
    }
    
    // dense1 <- copy(dense0)
    // dense0 <- reset(dense0)
    j = 0;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && j < m.size){
            dense1.value[j] = dense0.value[i];
            dense1.bitmap[j] = dense0.bitmap[i];
            j++;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }
    
    dense1.size = j;
    // dense0 <- geometric_product(dense1,sparse2)
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2->size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2->bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[dense1.bitmap[i]][sparse2->bitmap[j]];
            dense0.bitmap[bitmap] = bitmap;
            dense0.value[bitmap] += dense1.value[i]*sparse2->value[j]*sign;
        }
    }

    *sparse = sparse_remove_relative_small(dense0, ga->precision);
    if(sparse->size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return 0;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return 1;
}


static int unary_scalar_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    ScalarMultivector *scalar0 = (ScalarMultivector*)data0;
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    GradeMap gm = ga->gm;
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g) return 0;
    if(*g) *scalar = *scalar0; // Projecting to scalar
    else *scalar = 0;
    return 1;
}

static int unary_sparse_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;

    GradeMap gm = ga->gm;
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g)
        return 0;

    int size = 0;
    for(Py_ssize_t i = 0; i < sparse0->size; i++)
        if(g[gm.grade[sparse0->bitmap[i]]])
            size++;

    *sparse = init_sparse_empty(size--);
    if(sparse->size == -1){
        PyMem_RawFree(g);
        return 0;
    }

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        if(g[gm.grade[sparse0->bitmap[i]]]){
            sparse->value[size] = sparse0->value[i];
            sparse->bitmap[size] = sparse0->bitmap[i];
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    return 1;
}

static int unary_scalar_reverse(void *out, void *data0, PyAlgebraObject *ga){
    ScalarMultivector *scalar0 = (ScalarMultivector*)data0;
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    *scalar = *scalar0;
    return 1;
}

static int unary_sparse_reverse(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    GradeMap gm = ga->gm;
    *sparse = init_sparse_empty(sparse0->size);
    if(sparse->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        int grade = gm.grade[sparse0->bitmap[i]];
        int sign = (grade & 2) ? -1 : 1;
        sparse->value[i] = sign*sparse0->value[i];
        sparse->bitmap[i] = sparse0->bitmap[i];
    }
    return 1;
}

static int unary_sparse_dual(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    DualMap dm = ga->dm;
    *sparse = init_sparse_empty(sparse0->size);
    if(sparse->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        Py_ssize_t bitmap = sparse0->bitmap[i];
        sparse->value[i] = dm.sign[bitmap]*sparse0->value[i];
        sparse->bitmap[i] = dm.bitmap[bitmap];
    }

    return 1;
}

static int unary_sparse_undual(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    DualMap dm = ga->dm;
    *sparse = init_sparse_empty(sparse0->size);
    if(sparse->size == -1)
        return 0;
    int sign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        Py_ssize_t bitmap = sparse0->bitmap[i];
        sparse->value[i] = sign*dm.sign[bitmap]*sparse0->value[i];
        sparse->bitmap[i] = dm.bitmap[bitmap];
    }

    return 1;
}


static int unary_blades_scalaradd(void *out, void *data0, PyAlgebraObject *ga, ga_float value, int sign){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    SparseMultivector scalar = {.size = -1};
    Py_ssize_t scalarindex = -1;
    for(Py_ssize_t i = 0; i < blades0->size; i++){
        if(blades0->grade[i] == 0){
            scalarindex = i;
            break;
        }
    }

    if((scalarindex != -1 && ga_check_value_precision(ga,sign*blades0->data[scalarindex].value[0] + value)) ||
        (scalarindex == -1 && ga_check_value_precision(ga,value))){
        Py_ssize_t size = blades0->size;
        if(scalarindex != -1)
            size -= 1;

        *blades = init_blades_empty(size);
        if(blades->size == -1)
            return 0;

        Py_ssize_t k = 0;
        for(Py_ssize_t i = 0; i < blades0->size; i++){
            if(i == scalarindex) i++; // skip scalar index
            Py_ssize_t sizei = blades0->data[i].size;
            blades->data[k] = init_sparse_empty(sizei);
            blades->grade[k] = blades0->grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades->data[k].value[j] = sign*blades0->data[i].value[j];
                blades->data[k].bitmap[j] = blades0->data[i].bitmap[j];
            }
            k++;
        }
        return 1;
    }

    if(scalarindex != -1){
        *blades = init_blades_empty(blades0->size);
        if(blades->size == -1)
            return 0;

        for(Py_ssize_t i = 0; i < blades0->size; i++){
            Py_ssize_t sizei = blades0->data[i].size;
            blades->data[i] = init_sparse_empty(sizei);
            blades->grade[i] = blades0->grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades->data[i].value[j] = sign*blades0->data[i].value[j];
                blades->data[i].bitmap[j] = blades0->data[i].bitmap[j];
            }
        }
        *blades->data[scalarindex].value += value;
    }else{

        scalar = init_sparse_empty(1);
        if(scalar.size == -1)
            return 0;
        *blades = init_blades_empty(blades0->size + 1);
        if(blades->size == -1)
            return 0;

        *scalar.value = value;
        *scalar.bitmap = 0;
        *blades->data = scalar;

        for(Py_ssize_t i = 0; i < blades0->size; i++){
            Py_ssize_t sizei = blades0->data[i].size;
            blades->data[i+1] = init_sparse_empty(sizei);
            blades->grade[i+1] = blades0->grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades->data[i+1].value[j] = sign*blades0->data[i].value[j];
                blades->data[i+1].bitmap[j] = blades0->data[i].bitmap[j];
            }
        }
    }

    return 1;
}


static int unary_blades_scalarproduct(void *out, void *data0, PyAlgebraObject *ga, ga_float value){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *sparse = (BladesMultivector*)out;
    *sparse = init_blades_empty(blades0->size);
    if(sparse->size == -1) return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        sparse->data[i] = init_sparse_empty(blades0->data[i].size);
        for(Py_ssize_t j = 0; j < blades0->data[i].size; j++){
            sparse->data[i].bitmap[j] = blades0->data[i].bitmap[j];
            sparse->data[i].value[j] = value*blades0->data[i].value[j];
        }
    }

    return 1;
}

static int binary_blades_add(void *out, void *data0, void *data1,  PyAlgebraObject *ga, int sign){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades1 = (BladesMultivector*)data1;
    BladesMultivector *sparse = (BladesMultivector*)out;
    
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1)
        return 0;

    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga->precision;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        sub = blades0->data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sub.value[j];
        }
    }

    for(Py_ssize_t i = 0; i < blades1->size; i++){
        sub = blades1->data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sign*sub.value[j];
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    *sparse = sparse_dense_to_blades_sparse(dense,ga->gm);
    if(sparse->size == -1){
        sparse_free_(dense);
        return 0;
    }

    sparse_free_(dense);
    return 1;
}

static int binary_blades_product(void *out, void *data0, void *data1, PyAlgebraObject *ga, ProductType ptype){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades1 = (BladesMultivector*)data1;
    BladesMultivector *sparse = (BladesMultivector*)out;
    
    CliffordMap m = ga->product[ptype];
    GradeMap gm = ga->gm;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;

    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        ssparse0 = blades0->data[i];
        for(Py_ssize_t j = 0; j < blades1->size; j++){
            ssparse1 = blades1->data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = m.bitmap[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    *sparse = sparse_dense_to_blades_sparse(dense,gm);
    if(sparse->size == -1){
        sparse_free_(dense);
        return 0;
    }

    sparse_free_(dense);
    return 1;
}

static int ternary_blades_product(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades1 = (BladesMultivector*)data1;
    BladesMultivector *blades2 = (BladesMultivector*)data2;
    BladesMultivector *sparse = (BladesMultivector*)out;

    CliffordMap m = ga->product[ptype];
    GradeMap gm = ga->gm;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;

    if(dense0.size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        ssparse0 = blades0->data[i];
        for(Py_ssize_t j = 0; j < blades1->size; j++){
            ssparse1 = blades1->data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = m.bitmap[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return 0;
    }
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.bitmap[size] = dense0.bitmap[i];
            dense1.value[size] = dense0.value[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < blades2->size; j++){
            ssparse2 = blades2->data[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = m.bitmap[dense1.bitmap[i]][ssparse2.bitmap[k]];
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    *sparse = sparse_dense_to_blades_sparse(dense0,gm);
    if(sparse->size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return 0;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return 1;
}


static int unary_blades_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    Py_ssize_t *g = get_grade_bool(grades,grade_size,MAX_GRADE(ga) + 1);
    if(!g)
        return 0;

    int size = 0;
    for(Py_ssize_t i = 0; i < blades0->size; i++)
        if(g[blades0->grade[i]])
            size++;

    *blades = init_blades_empty(size--);
    if(blades->size == -1){
        PyMem_RawFree(g);
        return 0;
    }

    // copy the values of the selected grades
    for(Py_ssize_t i = 0; i < blades0->size; i++){
        int grade = blades0->grade[i];
        if(g[grade]){
            Py_ssize_t bsize = blades0->data[i].size;
            blades->data[size].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
            blades->data[size].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
            if(!blades->data[size].bitmap || !blades->data[size].value){
                blades_free_(*blades);
                PyMem_RawFree(g);
                return 0;
            }
            blades->data[size].size = bsize;
            blades->grade[size] = grade;
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades->data[size].bitmap[j] = blades0->data[i].bitmap[j];
                blades->data[size].value[j] = blades0->data[i].value[j];
            }
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    return 1;
}



static int unary_blades_reverse(void *out, void *data0, PyAlgebraObject *ga){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    *blades = init_blades_empty(blades0->size);
    if(blades->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        int grade = blades0->grade[i];
        Py_ssize_t bsize = blades0->data[i].size;
        blades->data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades->data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades->data[i].bitmap || !blades->data[i].value){
            blades_free_(*blades);
            return 0;
        }
        blades->data[i].size = bsize;
        blades->grade[i] = grade;

        if(grade & 2){
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades->data[i].bitmap[j] = blades0->data[i].bitmap[j];
                blades->data[i].value[j] = -blades0->data[i].value[j];
            }
        }else{
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades->data[i].bitmap[j] = blades0->data[i].bitmap[j];
                blades->data[i].value[j] = blades0->data[i].value[j];
            }
        }
    }

    return 1;
}

static int unary_blades_dual(void *out, void *data0, PyAlgebraObject *ga){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    DualMap dm = ga->dm;
    *blades = init_blades_empty(blades0->size);
    if(blades->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        Py_ssize_t bsize = blades0->data[i].size;
        blades->data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades->data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades->data[i].bitmap || !blades->data[i].value){
            blades_free_(*blades);
            return 0;
        }
        blades->data[i].size = bsize;
        blades->grade[i] = METRIC_SIZE(ga) - blades0->grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0->data[i].bitmap[j];
            blades->data[i].bitmap[j] = dm.bitmap[bitmap];
            blades->data[i].value[j] = dm.sign[bitmap]*blades0->data[i].value[j];
        }
    }

    return 1;
}

static int unary_blades_undual(void *out, void *data0, PyAlgebraObject *ga){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    DualMap dm = ga->dm;
    *blades = init_blades_empty(blades0->size);
    if(blades->size == -1)
        return 0;

    int sign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < blades0->size; i++){
        Py_ssize_t bsize = blades0->data[i].size;
        blades->data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades->data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades->data[i].bitmap || !blades->data[i].value){
            blades_free_(*blades);
            return 0;
        }
        blades->data[i].size = bsize;
        blades->grade[i] = METRIC_SIZE(ga) - blades0->grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0->data[i].bitmap[j];
            blades->data[i].bitmap[j] = dm.bitmap[bitmap];
            blades->data[i].value[j] = sign*dm.sign[bitmap]*blades0->data[i].value[j];
        }
    }

    return 1;
}


static int unary_dense_scalaradd(void *out, void *data0, PyAlgebraObject *ga, ga_float value, int sign){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;
    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[i] = sign*dense0->value[i];

    dense->value[0] += value; // increment the scalar part

    return 1;
}


static int unary_dense_scalarproduct(void *out, void *data0, PyAlgebraObject *ga, ga_float value){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1) return 0;
    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[i] = value*dense0->value[i];

    return 1;
}

static int binary_dense_add(void *out, void *data0, void *data1, PyAlgebraObject *ga, int sign){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense1 = (DenseMultivector*)data1;
    DenseMultivector *dense = (DenseMultivector*)out;
    
    Py_ssize_t ssize = -1;
    *dense = init_dense_empty(ga->asize);
    if(dense0->size == -1) return 0;

    // determine which resides in the smaller algebra
    if(dense0->size > dense1->size){
        ssize = dense1->size;
        for(Py_ssize_t i = ssize; i < dense0->size; i++)
            dense->value[i] = dense0->value[i];
    }else{
        ssize = dense0->size;
        for(Py_ssize_t i = ssize; i <  dense1->size; i++)
            dense->value[i] = sign*dense1->value[i];
    }

    for(Py_ssize_t i = 0; i < ssize; i++)
        dense->value[i] += dense0->value[i] + sign*dense1->value[i];

    return 1;
}

static int binary_dense_product(void *out, void *data0, void *data1, PyAlgebraObject *ga, ProductType ptype){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense1 = (DenseMultivector*)data1;
    DenseMultivector *dense = (DenseMultivector*)out;

    CliffordMap m = ga->product[ptype];
    *dense = init_dense_empty(m.size);
    if(dense->size == -1) return 0;
    int sign;
    for(Py_ssize_t i = 0; i < dense0->size; i++){
        for(Py_ssize_t j = 0; j < dense1->size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            dense->value[m.bitmap[i][j]] += dense0->value[i]*dense1->value[j]*sign;
        }
    }

    return 0;
}

static int ternary_dense_product(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense1 = (DenseMultivector*)data1;
    DenseMultivector *dense2 = (DenseMultivector*)data2;
    DenseMultivector *dense = (DenseMultivector*)out;

    CliffordMap m = ga->product[ptype];
    
    DenseMultivector temp = {.size = -1};
    *dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return 0;
    if(dense->size == -1) return 0;
    int sign;
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            temp.value[m.bitmap[i][j]] += dense0->value[i]*dense1->value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            dense->value[m.bitmap[i][j]] += temp.value[i]*dense2->value[j]*sign;
        }
    }
    dense_free_(temp);
    return 1;
}

static int unary_dense_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    GradeMap gm = ga->gm;
    
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g) return 0;
    *dense = init_dense_empty(dense0->size);
    for(Py_ssize_t i = 0; i < dense->size; i++)
        if(g[gm.grade[i]])
            dense->value[i] = dense0->value[i];

    PyMem_RawFree(g);
    return 1;
}

static int unary_dense_reverse(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    GradeMap gm = ga->gm;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < dense0->size; i++){
        int sign = (gm.grade[i] & 2) ? -1 : 1;
        dense->value[i] = sign*dense0->value[i];
    }

    return 1;
}

static int unary_dense_dual(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    DualMap dm = ga->dm;

    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[dm.bitmap[i]] = dm.sign[i]*dense0->value[i];

    return 1;
}

static int unary_dense_undual(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;

    DualMap dm = ga->dm;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    int sign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[dm.bitmap[i]] = sign*dm.sign[i]*dense0->value[i];

    return 1;
}

static int atomic_scalar_add(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize){
    ScalarMultivector *data = (ScalarMultivector*)data0;
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    *scalar = 0;
    for(Py_ssize_t j = 0; j < dsize; j++) *scalar += data[j];
    return 1;
}

static int atomic_sparse_add(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize){
    SparseMultivector *data = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;

    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return 0;

    for(Py_ssize_t j = 0; j < dsize; j++){
        for(Py_ssize_t i = 0; i < data[j].size; i++){
            dense.bitmap[data[j].bitmap[i]] = data[j].bitmap[i];
            dense.value[data[j].bitmap[i]] += data[j].value[i];
        }
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);
    if(sparse->size == -1){
        sparse_free_(dense);
        return 0;
    }

    sparse_free_(dense);
    return 1;
}

static int atomic_scalar_product(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize, ProductType ptype){
    ScalarMultivector *data = (ScalarMultivector*)data0;
    ScalarMultivector *scalar = (ScalarMultivector*)out;
    *scalar = 0;
    for(Py_ssize_t j = 0; j < dsize; j++) *scalar *= data[j];
    return 1;
}

static int atomic_sparse_product(void *out, void *data0, PyAlgebraObject *ga,Py_ssize_t dsize, ProductType ptype){
    SparseMultivector *data = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 1;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return 0;
    }

    Py_ssize_t tsize = 1;
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = m.bitmap[temp.bitmap[k]][data[i].bitmap[j]];
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    *sparse = sparse_remove_relative_small(temp,ga->precision);

    if(sparse->size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return 0;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return 1;
}
/*
static SparseMultivector unary_sparse_exponential_(SparseMultivector sparse0, PyAlgebraObject *ga){
    SparseMultivector exp = sparse_new_(sparse0.size + 1);
    SparseMultivector pow = sparse_copy_(sparse0);
    SparseMultivector prev;
    Py_ssize_t niters = 50;
    if(pow.size == -1 || exp.size == -1){
        sparse_free_(pow);
        sparse_free_(exp);
        pow.size = -1;
        return pow;
    }

    *exp.bitmap = 0;
    *exp.value = 1;
    memcpy(exp.bitmap+1,sparse0.bitmap,sparse0.size*sizeof(int));
    memcpy(exp.value+1,sparse0.value,sparse0.size*sizeof(ga_float));

    for(Py_ssize_t k = 2; k < niters+2; k++){
        prev = pow;
        pow = binary_sparse_product_(pow,sparse0,ga,ProductType_geometric);
        sparse_free_(prev); prev = pow;
        pow = unary_sparse_scalarproduct_(pow,ga,1.0/(ga_float)k);
        sparse_free_(prev); prev = exp;
        exp = binary_sparse_add_(exp,pow,ga,1);
        sparse_free_(prev);
    }

    sparse_free_(pow);
    return exp;
}

*/
static int atomic_blades_add(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t size){
    BladesMultivector *data = (BladesMultivector*)data0;
    BladesMultivector *sparse = (BladesMultivector*)out;
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return 0;
    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga->precision;

    for(Py_ssize_t k = 0; k < size; k++){
        for(Py_ssize_t i = 0; i < data[k].size; i++){
            sub = data[k].data[i];
            for(Py_ssize_t j = 0; j < sub.size; j++){
                bitmap = sub.bitmap[j];
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += sub.value[j];
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;
    
    sparse_remove_rel_small(dense, ga->precision);
    *sparse = sparse_dense_to_blades_sparse(dense,ga->gm);
    sparse_free_(dense);
    return 1;
}

static int atomic_blades_product(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize, ProductType ptype){
    BladesMultivector *data = (BladesMultivector*)data0;
    BladesMultivector *sparse = (BladesMultivector*)out;
    CliffordMap m = ga->product[ptype];
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return 0;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = m.bitmap[temp.bitmap[j]][sdata.bitmap[l]];
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += temp.value[j]*sdata.value[l]*sign;
                }
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_rel_small(temp,ga->precision);
    *sparse = sparse_dense_to_blades_sparse(temp,ga->gm);
    sparse_free_(dense);
    sparse_free_(temp);
    return 1;
}

static int atomic_dense_add(void *out, void *data0,PyAlgebraObject *ga, Py_ssize_t size){
    DenseMultivector *data = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    *dense = init_dense_empty(ga->asize);
    if(dense->size == -1) return 0;
    ga_float value;
    for(Py_ssize_t i = 0; i < dense->size; i++){
        value = 0;
        for(Py_ssize_t k = 0; k < size; k++)
            value += data[k].value[i];
        dense->value[i] += value;
    }

    return 1;
}

static int atomic_dense_product(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize, ProductType ptype){
    DenseMultivector *data = (DenseMultivector*)data0;
    DenseMultivector *dense_out = (DenseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return 0;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return 0;
    }

    *temp.value = 1; // initialize temp to unit scalar
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                dense.value[m.bitmap[k][j]] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    *dense_out = temp;
    return 1;
}


static int binary_mixed_add(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, int sign){
    SparseMultivector *sparse = (SparseMultivector*)out;
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return 0;

    Py_ssize_t size = 0;

    while(iter0->next(iter0)){
        if(dense.bitmap[iter0->bitmap] == -1){
            dense.bitmap[iter0->bitmap] = iter0->bitmap;
            size++;
        }
        dense.value[iter0->bitmap] += iter0->value;
    }

    while(iter1->next(iter1)){
        if(dense.bitmap[iter1->bitmap] == -1){
            dense.bitmap[iter1->bitmap] = iter1->bitmap;
            size++;
        }
        dense.value[iter1->bitmap] += sign*iter1->value;
    }

    size = sparse_remove_rel_small(dense,ga->precision);
    *sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return 1;
}

static int binary_mixed_product(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    SparseMultivector *sparse = (SparseMultivector*)out;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = m.bitmap[iter0->bitmap][iter1->bitmap];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);
    sparse_free_(dense);
    return 1;
}

int is_bigger_metric(PyAlgebraObject *ga0, PyAlgebraObject *ga1){
    Py_ssize_t size = METRIC_SIZE(ga0) < METRIC_SIZE(ga1) ?  METRIC_SIZE(ga0) :  METRIC_SIZE(ga1);
    for(Py_ssize_t i = 0; i < size; i++)
        if(ga0->metric[i] != ga1->metric[i])
            return -1;
    return METRIC_SIZE(ga0) > METRIC_SIZE(ga1);
}


static int atomic_mixed_add(void* out, PyMultivectorIter *iter, PyAlgebraObject *ga, Py_ssize_t dsize){
    SparseMultivector *sparse = (SparseMultivector*)out;
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return 0;
    
    for(Py_ssize_t j = 0; j < dsize; j++){
        while(iter->next(iter)){
            dense.bitmap[iter->bitmap] = iter->bitmap;
            dense.value[iter->bitmap] += iter->value;
        }iter++;
    }

    *sparse = sparse_remove_relative_small(dense, ga->precision);
    sparse_free_(dense);
    return 1;
}



static int atomic_mixed_product(void *out, PyMultivectorIter *iter, PyAlgebraObject *ga, Py_ssize_t size, ProductType ptype){
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return 0;
    }

    
    Py_ssize_t tsize = 1;
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = m.bitmap[temp.bitmap[k]][iter->bitmap];
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*iter->value*sign;
            }
        }iter++;
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    *sparse = sparse_remove_relative_small(temp, ga->precision);
    sparse_free_(dense);
    sparse_free_(temp);
    return 1;
}

void sparse_free(void *sparse){
    if(!sparse)
        return;
    sparse_free_(*((SparseMultivector*)sparse));
}

 void blades_free(void *blades){
    if(!blades)
        return;
    blades_free_(*((BladesMultivector*)blades));
}

void dense_free(void *dense){
    if(!dense)
        return;
    dense_free_(*((DenseMultivector*)dense));
}

void scalar_free(void *scalar){
    return;
}

static PyMultivectorData_Funcs multivector_sparse_data_fn = {
    .free = sparse_free,
    .init = sparse_init,
    .iter_init = sparse_iter_init,
    .iter_next = sparse_iter_next,
    .cast = cast_to_sparse,
};

static PyMultivectorData_Funcs multivector_blades_data_fn = {
    .free = blades_free,
    .init = blades_init,
    .iter_init = blades_iter_init,
    .iter_next = blades_iter_next,
    .cast = cast_to_blades,
};
 
static PyMultivectorData_Funcs multivector_dense_data_fn = {
    .free = dense_free,
    .init = dense_init,
    .iter_init = dense_iter_init,
    .iter_next = dense_iter_next,
    .cast = cast_to_dense,
};

static PyMultivectorData_Funcs multivector_scalar_data_fn = {
    .free = scalar_free,
    .init = scalar_init,
    .iter_init = scalar_iter_init,
    .iter_next = scalar_iter_next,
    .cast = NULL,
};

static PyMultivectorMath_Funcs multivector_sparse_math_fn = {
    .product = binary_sparse_product,
    .atomic_product = atomic_sparse_product,
    .ternary_product = ternary_sparse_product,
    .add = binary_sparse_add,
    .atomic_add = atomic_sparse_add,
    .grade_project = unary_sparse_gradeproject,
    .scalar_product = unary_sparse_scalarproduct,
    .scalar_add = unary_sparse_scalaradd,
    .reverse = unary_sparse_reverse,
    .dual = unary_sparse_dual,
    .undual = unary_sparse_undual,
};

 static PyMultivectorMath_Funcs multivector_blades_math_fn = {
    .product = binary_blades_product,
    .atomic_product = atomic_blades_product,
    .ternary_product = ternary_blades_product,
    .add = binary_blades_add,
    .atomic_add = atomic_blades_add,
    .grade_project = unary_blades_gradeproject,
    .scalar_product = unary_blades_scalarproduct,
    .scalar_add = unary_blades_scalaradd,
    .reverse = unary_blades_reverse,
    .dual = unary_blades_dual,
    .undual = unary_blades_undual,
};
 static PyMultivectorMath_Funcs multivector_dense_math_fn = {
    .product = binary_dense_product,
    .atomic_product = atomic_dense_product,
    .ternary_product = ternary_dense_product,
    .add = binary_dense_add,
    .atomic_add = atomic_dense_add,
    .grade_project = unary_dense_gradeproject,
    .scalar_product = unary_dense_scalarproduct,
    .scalar_add = unary_dense_scalaradd,
    .reverse = unary_dense_reverse,
    .dual = unary_dense_dual,
    .undual = unary_dense_undual,
};

static PyMultivectorMath_Funcs multivector_scalar_math_fn = {
    .product = binary_scalar_product,
    .atomic_product = atomic_scalar_product,
    .ternary_product = ternary_scalar_product,
    .add = binary_scalar_add,
    .atomic_add = atomic_scalar_add,
    .grade_project = unary_scalar_gradeproject,
    .scalar_product = unary_scalar_scalarproduct,
    .scalar_add = unary_scalar_scalaradd,
    .reverse = unary_scalar_reverse,
    .dual = NULL,
    .undual = NULL,
};

static const PyMultivectorSubType sparse_subtype = {
    .math_funcs = &multivector_sparse_math_fn,
    .data_funcs = &multivector_sparse_data_fn,
    .name = "",
    .type_name = "sparse",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_sparse,
    .basic_size = sizeof(SparseMultivector),
};
static const PyMultivectorSubType blades_subtype = {
    .math_funcs = &multivector_blades_math_fn,
    .data_funcs = &multivector_blades_data_fn,
    .name = "",
    .type_name = "blades",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_blades,
    .basic_size = sizeof(BladesMultivector),
};
static const PyMultivectorSubType dense_subtype = {
    .math_funcs = &multivector_dense_math_fn,
    .data_funcs = &multivector_dense_data_fn,
    .name = "",
    .type_name = "dense",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_dense,
    .basic_size = sizeof(DenseMultivector),
};

static const PyMultivectorSubType scalar_subtype = {
    .math_funcs = &multivector_scalar_math_fn,
    .data_funcs = &multivector_scalar_data_fn,
    .name = "",
    .type_name = "scalar",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_scalar,
    .basic_size = sizeof(ScalarMultivector),
};

PyMultivectorMixedMath_Funcs multivector_mixed_fn = {
  .add = binary_mixed_add,
  .product = binary_mixed_product,
  .atomic_add = atomic_mixed_add,
  .atomic_product = atomic_mixed_product,
  .type_names = {"blades","sparse","dense","scalar",NULL},
};


// PyMultivectorSubType multivector_subtypes_array[4] = {sparse_subtype,dense_subtype,blades_subtype,scalar_subtype};

PyMultivectorSubType multivector_subtypes_array[4] = {{
    .math_funcs = &multivector_sparse_math_fn,
    .data_funcs = &multivector_sparse_data_fn,
    .name = "",
    .type_name = "sparse",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_sparse,
    .basic_size = sizeof(SparseMultivector),
},{
    .math_funcs = &multivector_dense_math_fn,
    .data_funcs = &multivector_dense_data_fn,
    .name = "",
    .type_name = "dense",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_dense,
    .basic_size = sizeof(DenseMultivector),
},{
    .math_funcs = &multivector_blades_math_fn,
    .data_funcs = &multivector_blades_data_fn,
    .name = "",
    .type_name = "blades",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_blades,
    .basic_size = sizeof(BladesMultivector),
},{
    .math_funcs = &multivector_scalar_math_fn,
    .data_funcs = &multivector_scalar_data_fn,
    .name = "",
    .type_name = "scalar",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_scalar,
    .basic_size = sizeof(ScalarMultivector),
}};

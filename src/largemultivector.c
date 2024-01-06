#include <python3.11/Python.h>
#include "gasparse.h"


static Py_ssize_t *init_grade_size(PyAlgebraObject *ga){
    Py_ssize_t *gsize = (Py_ssize_t*)PyMem_RawMalloc((MAX_GRADE(ga)+1)*sizeof(Py_ssize_t));
    if(!gsize){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for grade size array");
        return NULL;
    }
    for(Py_ssize_t i = 0; i <= MAX_GRADE(ga); i++)
        gsize[i] = 0;
    return gsize;
}


static BladesMultivector sparse_dense_to_blades_sparse(SparseMultivector dense, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    Py_ssize_t ssize = 0, grade = -1;
    Py_ssize_t *gsize = init_grade_size(ga);
    Py_ssize_t *gindex = init_grade_size(ga);
    if(!gsize || !gindex){
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        return sparse;
    }
    int bitmap;
    for(Py_ssize_t i = 0; i < dense.size; i++){
        if(dense.bitmap[i] == -1) continue;
        grade = GRADE(dense.bitmap[i]);
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
    for(Py_ssize_t i = 0; i <= MAX_GRADE(ga); i++){ // iterate over grades
        if(!gsize[i]) continue;
        sparse.data[gindex[i]] = init_sparse_empty(gsize[i]);
        sparse.grade[gindex[i]] = i;
    }

    for(Py_ssize_t i = 0; i < dense.size; i++){
        bitmap = dense.bitmap[i];
        if(bitmap == -1) continue;
        grade = GRADE(bitmap); gsize[grade]--;
        sparse.data[gindex[grade]].bitmap[gsize[grade]] = bitmap;
        sparse.data[gindex[grade]].value[gsize[grade]] = dense.value[i];
    }

    PyMem_RawFree(gsize);
    PyMem_RawFree(gindex);
    return sparse;
}



static BladesMultivector blades_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    if(!size){
        BladesMultivector blades = {.size = 0,.grade = NULL, .data = NULL};
        return blades;
    }
    SparseMultivector ssparse = {.bitmap = bitmap, .value = value, .size = size};
    return sparse_dense_to_blades_sparse(ssparse,ga);
}

static int cast_to_blades(PyMultivectorObject *data, PyMultivectorObject *to){
    PyMultivectorIter *iter = init_multivector_iter(data,1);
    BladesMultivector *pblades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    
    if(!iter || !pblades || !to){
        free_multivector_iter(iter,1);
        PyMem_RawFree(pblades);
        return 0;
    }

    SparseMultivector sparse = {.size = iter->niters, .value = NULL, .bitmap = NULL};
    sparse.value = (ga_float*)PyMem_RawMalloc(iter->niters*sizeof(ga_float));
    sparse.bitmap = (int*)PyMem_RawMalloc(iter->niters*sizeof(int));
    Py_ssize_t i = 0;
    while(iter->next(iter)){
        sparse.value[i] = iter->value;
        sparse.bitmap[i] = iter->bitmap;
        i++;
    }
    *pblades = sparse_dense_to_blades_sparse(sparse,data->GA);
    sparse_free_(sparse);
    to->data = (void*)pblades;
    free_multivector_iter(iter,1);
    return 1;
}




static SparseMultivector binary_sparse_geometricproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector ternary_sparse_geometricproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_outerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector ternary_sparse_outerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            if(GRADE(dense1.bitmap[i])+GRADE(sparse2.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_innerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(abs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector ternary_sparse_innerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(abs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            if(abs((_grade0=GRADE(dense1.bitmap[i]))-(_grade1=GRADE(sparse2.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_regressiveproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    Py_ssize_t _grade0;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap,inner_bitmap;
    int undualsign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    Py_ssize_t pss = ga->asize - 1;
    Py_ssize_t l,r; int lsign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        l = pss ^ sparse0.bitmap[i];
        _grade0 = GRADE(l);
        lsign = undualsign*dm.sign[sparse0.bitmap[i]];
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            r = pss ^ sparse1.bitmap[j];
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*m.sign[l][r]*dm.sign[sparse1.bitmap[j]]*lsign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector unary_sparse_gradeproject_(SparseMultivector sparse0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    SparseMultivector sparse = {.size = -1};
    Py_ssize_t *g = get_grade_bool(grades,grade_size,MAX_GRADE(ga) + 1);
    if(!g)
        return sparse;

    int size = 0;
    for(Py_ssize_t i = 0; i < sparse0.size; i++)
        if(g[GRADE(sparse0.bitmap[i])])
            size++;

    sparse = init_sparse_empty(size--);
    if(sparse.size == -1){
        PyMem_RawFree(g);
        return sparse;
    }

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        if(g[GRADE(sparse0.bitmap[i])]){
            sparse.value[size] = sparse0.value[i];
            sparse.bitmap[size] = sparse0.bitmap[i];
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    return sparse;
}



static SparseMultivector unary_sparse_reverse_(SparseMultivector sparse0, PyAlgebraObject *ga){
    SparseMultivector sparse = init_sparse_empty(sparse0.size);
    if(sparse.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        int sign = (GRADE(sparse0.bitmap[i]) & 2) ? -1 : 1;
        sparse.value[i] = sign*sparse0.value[i];
        sparse.bitmap[i] = sparse0.bitmap[i];
    }

    return sparse;
}

static SparseMultivector unary_sparse_dual_(SparseMultivector sparse0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    SparseMultivector sparse = init_sparse_empty(sparse0.size);
    if(sparse.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        Py_ssize_t bitmap = sparse0.bitmap[i];
        sparse.value[i] = dm.sign[bitmap]*sparse0.value[i];
        sparse.bitmap[i] = pss ^ bitmap;
    }

    return sparse;
}

static SparseMultivector unary_sparse_undual_(SparseMultivector sparse0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    SparseMultivector sparse = init_sparse_empty(sparse0.size);
    if(sparse.size == -1)
        return sparse;
    int sign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        Py_ssize_t bitmap = sparse0.bitmap[i];
        sparse.value[i] = sign*dm.sign[bitmap]*sparse0.value[i];
        sparse.bitmap[i] = pss ^ bitmap;
    }

    return sparse;
}


static BladesMultivector binary_blades_add_(BladesMultivector blades0, BladesMultivector blades1,  PyAlgebraObject *ga, int sign){
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1)
        return sparse;

    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga->precision;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        sub = blades0.data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sub.value[j];
        }
    }

    for(Py_ssize_t i = 0; i < blades1.size; i++){
        sub = blades1.data[i];
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

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}



static BladesMultivector binary_blades_geometricproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    
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

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_geometricproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
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
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_outerproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(grade0+grade1!=GRADE(bitmap)) continue;
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

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_outerproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    Py_ssize_t grade0,grade1,grade2;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(grade0+grade1!=GRADE(bitmap)) continue;
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
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
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            grade2 = blades2.grade[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                if(GRADE(dense1.bitmap[i])+grade2!=GRADE(bitmap)) continue;
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_innerproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(abs(grade0-grade1)!=GRADE(bitmap)||!grade0||!grade1) continue;
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

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_innerproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    Py_ssize_t grade0,grade1,grade2;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    Py_ssize_t _grade0;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(abs(grade0-grade1)!=GRADE(bitmap)||!grade0||!grade1) continue;
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
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
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            grade2 = blades2.grade[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                if(abs((_grade0=GRADE(dense1.bitmap[i]))-grade2)!=GRADE(bitmap)||!_grade0||!grade2) continue;
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_regressiveproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize - 1;
    ga_float precision = ga->precision;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    Py_ssize_t l,r; int lsign;
    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    Py_ssize_t bitmap,inner_bitmap;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = max_grade - blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = max_grade - blades1.grade[j];
            for(Py_ssize_t n = 0; n < ssparse0.size; n++){
                l = pss ^ ssparse0.bitmap[n]; // product with the pseudoscalar
                lsign = undualsign*dm.sign[ssparse0.bitmap[n]];
                for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                    r = pss ^ ssparse1.bitmap[k];
                    bitmap = pss^(inner_bitmap = l ^ r);
                    if(grade0 + grade1 != GRADE(inner_bitmap)) continue;
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[n]*ssparse1.value[k]*m.sign[l][r]*dm.sign[ssparse1.bitmap[k]]*lsign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}


static BladesMultivector unary_blades_dual_(BladesMultivector blades0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    BladesMultivector blades = init_blades_empty(blades0.size);
    if(blades.size == -1)
        return blades;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        Py_ssize_t bsize = blades0.data[i].size;
        blades.data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades.data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades.data[i].bitmap || !blades.data[i].value){
            blades_free_(blades);
            blades.size = -1;
            return blades;
        }
        blades.data[i].size = bsize;
        blades.grade[i] = max_grade - blades0.grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0.data[i].bitmap[j];
            blades.data[i].bitmap[j] = pss^bitmap;
            blades.data[i].value[j] = dm.sign[bitmap]*blades0.data[i].value[j];
        }
    }

    return blades;
}

static BladesMultivector unary_blades_undual_(BladesMultivector blades0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    BladesMultivector blades = init_blades_empty(blades0.size);
    if(blades.size == -1)
        return blades;
    Py_ssize_t pss = ga->asize-1;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    int sign = max_grade & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < blades0.size; i++){
        Py_ssize_t bsize = blades0.data[i].size;
        blades.data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades.data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades.data[i].bitmap || !blades.data[i].value){
            blades_free_(blades);
            blades.size = -1;
            return blades;
        }
        blades.data[i].size = bsize;
        blades.grade[i] = max_grade - blades0.grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0.data[i].bitmap[j];
            blades.data[i].bitmap[j] = pss^bitmap;
            blades.data[i].value[j] = sign*dm.sign[bitmap]*blades0.data[i].value[j];
        }
    }

    return blades;
}



static DenseMultivector binary_dense_geometricproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_geometricproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_outerproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_outerproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_innerproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    Py_ssize_t _grade0, _grade1;

    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(abs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_innerproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    Py_ssize_t _grade0, _grade1;

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(abs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(abs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_regressiveproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap,inner_bitmap;
    Py_ssize_t _grade0;
    Py_ssize_t pss = ga->asize - 1;
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;

    Py_ssize_t l,r; int lsign;
    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    for(Py_ssize_t i = 0; i < dense0.size; i++){
        l = pss ^ i; // dual of the left multivector
        lsign = undualsign*dm.sign[i];
        _grade0 = GRADE(l);
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            r = pss ^ j; // dual of the right multivector
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*m.sign[l][r]*dm.sign[j]*lsign;
        }
    }

    return dense;
}

static DenseMultivector unary_dense_gradeproject_(DenseMultivector dense0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    DenseMultivector dense = {.size = -1};
    Py_ssize_t *g = get_grade_bool(grades,grade_size,MAX_GRADE(ga)+1);
    if(!g) return dense;
    dense = init_dense_empty(dense0.size);
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(g[GRADE(i)])
            dense.value[i] = dense0.value[i];

    PyMem_RawFree(g);
    return dense;
}

static DenseMultivector unary_dense_reverse_(DenseMultivector dense0, PyAlgebraObject *ga){
    DenseMultivector dense = init_dense_empty(dense0.size);
    if(dense.size == -1)
        return dense;

    for(Py_ssize_t i = 0; i < dense0.size; i++){
        int sign = (GRADE(i) & 2) ? -1 : 1;
        dense.value[i] = sign*dense0.value[i];
    }

    return dense;
}

static DenseMultivector unary_dense_dual_(DenseMultivector dense0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    DenseMultivector dense = init_dense_empty(dense0.size);
    if(dense.size == -1)
        return dense;

    for(Py_ssize_t i = 0; i < dense0.size; i++)
        dense.value[pss^i] = dm.sign[i]*dense0.value[i];

    return dense;
}

static DenseMultivector unary_dense_undual_(DenseMultivector dense0, PyAlgebraObject *ga){
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    DenseMultivector dense = init_dense_empty(dense0.size);
    if(dense.size == -1)
        return dense;

    int sign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        dense.value[pss^i] = sign*dm.sign[i]*dense0.value[i];

    return dense;
}




static SparseMultivector atomic_sparse_geometricproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static SparseMultivector atomic_sparse_outerproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                if(GRADE(temp.bitmap[k])+GRADE(data[i].bitmap[j])!=GRADE(bitmap)) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static SparseMultivector atomic_sparse_innerproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    Py_ssize_t _grade0, _grade1;

    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                if(abs((_grade0=GRADE(temp.bitmap[k]))-(_grade1=GRADE(data[i].bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static BladesMultivector atomic_blades_add_(BladesMultivector *data, Py_ssize_t size, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return sparse;
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

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    sparse_free_(dense);
    return sparse;
}



static BladesMultivector atomic_blades_geometricproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
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
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static BladesMultivector atomic_blades_outerproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    Py_ssize_t sgrade;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                sgrade = data[i].grade[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        if(GRADE(temp.bitmap[j])+sgrade!=GRADE(bitmap)) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static BladesMultivector atomic_blades_innerproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    Py_ssize_t sgrade;
    Py_ssize_t _grade0;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                sgrade = data[i].grade[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        if(abs((_grade0=GRADE(temp.bitmap[j]))-sgrade)!=GRADE(bitmap)||!_grade0||!sgrade) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}

static DenseMultivector atomic_dense_geometricproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    
    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}
static DenseMultivector atomic_dense_outerproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    
    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                if(GRADE(k)+GRADE(j)!=GRADE(bitmap)) continue;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}
static DenseMultivector atomic_dense_innerproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                if(abs((_grade0=GRADE(k))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}

static SparseMultivector binary_mixed_regressiveproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize - 1;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    Py_ssize_t _grade0;
    SparseMultivector sparse;
    Py_ssize_t size = 0;
    Py_ssize_t l,r; int lsign;

    Py_ssize_t bitmap, inner_bitmap;
    while(iter0->next(iter0)){
        l = pss^iter0->bitmap;
        lsign = undualsign*dm.sign[iter0->bitmap];
        _grade0 = GRADE(l);
        while(iter1->next(iter1)){
            r = pss^iter1->bitmap;
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*m.sign[l][r]*lsign*dm.sign[iter1->bitmap];
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}



static SparseMultivector binary_mixed_geometricproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    
    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_geometricproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static SparseMultivector binary_mixed_outerproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    
    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            if(GRADE(iter0->bitmap)+GRADE(iter1->bitmap)!=GRADE(bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_outerproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                if(GRADE(temp.bitmap[k])+GRADE(iter->bitmap)!=GRADE(bitmap)) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static SparseMultivector binary_mixed_innerproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    Py_ssize_t _grade0, _grade1;

    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            if(abs((_grade0=GRADE(iter0->bitmap))-(_grade1=GRADE(iter1->bitmap)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_innerproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    Py_ssize_t _grade0, _grade1;

    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                if(abs((_grade0=GRADE(temp.bitmap[k]))-(_grade1=GRADE(iter->bitmap)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


   
static PyMultivectorObject *unary_sparse_gradeproject(PyMultivectorObject *data0, int *grades, Py_ssize_t size){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_gradeproject_(*sparse0,data0->GA, grades, size);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_reverse(PyMultivectorObject *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_reverse_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_dual(PyMultivectorObject *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_dual_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_undual(PyMultivectorObject *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_undual_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

   
static PyMultivectorObject *unary_dense_gradeproject(PyMultivectorObject *data0, int *grades, Py_ssize_t size){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_gradeproject_(*dense0,data0->GA, grades, size);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_reverse(PyMultivectorObject *data0){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_reverse_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_dual(PyMultivectorObject *data0){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_dual_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_undual(PyMultivectorObject *data0){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_undual_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

     
static PyMultivectorObject *unary_blades_dual(PyMultivectorObject *data0){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_dual_(*blades0,data0->GA);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_blades_undual(PyMultivectorObject *data0){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_undual_(*blades0,data0->GA);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}



  static PyMultivectorObject *binary_sparse_product(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)data0->data;
    SparseMultivector *psparse1 = (SparseMultivector*)data1->data;
    SparseMultivector *psparse  = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!psparse0 ||!psparse1 ||!psparse || !out){
        PyMem_RawFree(psparse);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = binary_sparse_geometricproduct_(*psparse0,*psparse1,data0->GA);
            break;
        case ProductType_inner:
            *psparse = binary_sparse_innerproduct_(*psparse0,*psparse1,data0->GA);
            break;
        case ProductType_outer:
            *psparse = binary_sparse_outerproduct_(*psparse0,*psparse1,data0->GA);
            break;
        case ProductType_regressive:
            *psparse = binary_sparse_regressiveproduct_(*psparse0,*psparse1,data0->GA);
            break;
        default:
            PyMem_RawFree(psparse);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)psparse;
    Py_SET_REFCNT(out,1);
    return out;
}
 static PyMultivectorObject *ternary_sparse_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)data0->data;
    SparseMultivector *psparse1 = (SparseMultivector*)data1->data;
    SparseMultivector *psparse2 = (SparseMultivector*)data2->data;
    SparseMultivector *psparse  = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!psparse0 ||!psparse1 ||!psparse2 ||!psparse || !out){
        PyMem_RawFree(psparse);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = ternary_sparse_geometricproduct_(*psparse0,*psparse1,*psparse2,data0->GA);
            break;
        case ProductType_inner:
            *psparse = ternary_sparse_innerproduct_(*psparse0,*psparse1,*psparse2,data0->GA);
            break;
        case ProductType_outer:
            *psparse = ternary_sparse_outerproduct_(*psparse0,*psparse1,*psparse2,data0->GA);
            break;
        default:
            PyMem_RawFree(psparse);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)psparse;
    Py_SET_REFCNT(out,1);
    return out;
}
  static PyMultivectorObject *binary_dense_product(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)data0->data;
    DenseMultivector *pdense1 = (DenseMultivector*)data1->data;
    DenseMultivector *pdense  = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pdense0 ||!pdense1 ||!pdense || !out){
        PyMem_RawFree(pdense);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = binary_dense_geometricproduct_(*pdense0,*pdense1,data0->GA);
            break;
        case ProductType_inner:
            *pdense = binary_dense_innerproduct_(*pdense0,*pdense1,data0->GA);
            break;
        case ProductType_outer:
            *pdense = binary_dense_outerproduct_(*pdense0,*pdense1,data0->GA);
            break;
        case ProductType_regressive:
            *pdense = binary_dense_regressiveproduct_(*pdense0,*pdense1,data0->GA);
            break;
        default:
            PyMem_RawFree(pdense);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
 static PyMultivectorObject *ternary_dense_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)data0->data;
    DenseMultivector *pdense1 = (DenseMultivector*)data1->data;
    DenseMultivector *pdense2 = (DenseMultivector*)data2->data;
    DenseMultivector *pdense  = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pdense0 ||!pdense1 ||!pdense2 ||!pdense || !out){
        PyMem_RawFree(pdense);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = ternary_dense_geometricproduct_(*pdense0,*pdense1,*pdense2,data0->GA);
            break;
        case ProductType_inner:
            *pdense = ternary_dense_innerproduct_(*pdense0,*pdense1,*pdense2,data0->GA);
            break;
        case ProductType_outer:
            *pdense = ternary_dense_outerproduct_(*pdense0,*pdense1,*pdense2,data0->GA);
            break;
        default:
            PyMem_RawFree(pdense);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
  static PyMultivectorObject *binary_blades_product(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)data0->data;
    BladesMultivector *pblades1 = (BladesMultivector*)data1->data;
    BladesMultivector *pblades  = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pblades0 ||!pblades1 ||!pblades || !out){
        PyMem_RawFree(pblades);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = binary_blades_geometricproduct_(*pblades0,*pblades1,data0->GA);
            break;
        case ProductType_inner:
            *pblades = binary_blades_innerproduct_(*pblades0,*pblades1,data0->GA);
            break;
        case ProductType_outer:
            *pblades = binary_blades_outerproduct_(*pblades0,*pblades1,data0->GA);
            break;
        case ProductType_regressive:
            *pblades = binary_blades_regressiveproduct_(*pblades0,*pblades1,data0->GA);
            break;
        default:
            PyMem_RawFree(pblades);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}
 static PyMultivectorObject *ternary_blades_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)data0->data;
    BladesMultivector *pblades1 = (BladesMultivector*)data1->data;
    BladesMultivector *pblades2 = (BladesMultivector*)data2->data;
    BladesMultivector *pblades  = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pblades0 ||!pblades1 ||!pblades2 ||!pblades || !out){
        PyMem_RawFree(pblades);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = ternary_blades_geometricproduct_(*pblades0,*pblades1,*pblades2,data0->GA);
            break;
        case ProductType_inner:
            *pblades = ternary_blades_innerproduct_(*pblades0,*pblades1,*pblades2,data0->GA);
            break;
        case ProductType_outer:
            *pblades = ternary_blades_outerproduct_(*pblades0,*pblades1,*pblades2,data0->GA);
            break;
        default:
            PyMem_RawFree(pblades);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_blades_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;
    BladesMultivector *blades1 = (BladesMultivector*)data1->data;
    PyMultivectorObject *out = new_multivector_inherit_type(data0);// pass -1 to inherit type
    if(!blades || !blades0 || !blades1 || !out){
        multivector_dealloc(out);
        PyMem_RawFree(blades);
    }
    *blades = binary_blades_add_(*blades0,*blades1,data0->GA, sign);

    if(blades->size == -1){
        PyMem_RawFree(blades);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = blades;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static PyMultivectorObject *atomic_sparse_product(PyMultivectorObject *data0, Py_ssize_t size, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    SparseMultivector *psparse  = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!psparse0 || !psparse || !out){
        PyMem_RawFree(psparse);
        PyMem_RawFree(psparse0);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    for(Py_ssize_t i = 0; i < size; i++)
        psparse0[i] = *((SparseMultivector*)data0[i].data);

    switch(ptype){
        case ProductType_geometric:
            *psparse = atomic_sparse_geometricproduct_(psparse0,size,data0->GA);
            break;
        case ProductType_inner:
            *psparse = atomic_sparse_innerproduct_(psparse0,size,data0->GA);
            break;
        case ProductType_outer:
            *psparse = atomic_sparse_outerproduct_(psparse0,size,data0->GA);
            break;
        default:
            PyMem_RawFree(psparse);
            PyMem_RawFree(psparse0);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)psparse;
    Py_SET_REFCNT(out,1);
    PyMem_RawFree(psparse0);
    return out;
}
static PyMultivectorObject *atomic_dense_product(PyMultivectorObject *data0, Py_ssize_t size, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    DenseMultivector *pdense  = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pdense0 || !pdense || !out){
        PyMem_RawFree(pdense);
        PyMem_RawFree(pdense0);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    for(Py_ssize_t i = 0; i < size; i++)
        pdense0[i] = *((DenseMultivector*)data0[i].data);

    switch(ptype){
        case ProductType_geometric:
            *pdense = atomic_dense_geometricproduct_(pdense0,size,data0->GA);
            break;
        case ProductType_inner:
            *pdense = atomic_dense_innerproduct_(pdense0,size,data0->GA);
            break;
        case ProductType_outer:
            *pdense = atomic_dense_outerproduct_(pdense0,size,data0->GA);
            break;
        default:
            PyMem_RawFree(pdense);
            PyMem_RawFree(pdense0);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    PyMem_RawFree(pdense0);
    return out;
}
static PyMultivectorObject *atomic_blades_product(PyMultivectorObject *data0, Py_ssize_t size, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    BladesMultivector *pblades  = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data0);
    if(!pblades0 || !pblades || !out){
        PyMem_RawFree(pblades);
        PyMem_RawFree(pblades0);
        multivector_dealloc(out);
        return NULL; // raise error
    }

    for(Py_ssize_t i = 0; i < size; i++)
        pblades0[i] = *((BladesMultivector*)data0[i].data);

    switch(ptype){
        case ProductType_geometric:
            *pblades = atomic_blades_geometricproduct_(pblades0,size,data0->GA);
            break;
        case ProductType_inner:
            *pblades = atomic_blades_innerproduct_(pblades0,size,data0->GA);
            break;
        case ProductType_outer:
            *pblades = atomic_blades_outerproduct_(pblades0,size,data0->GA);
            break;
        default:
            PyMem_RawFree(pblades);
            PyMem_RawFree(pblades0);
            multivector_dealloc(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    PyMem_RawFree(pblades0);
    return out;
}
  static PyMultivectorObject *atomic_blades_add(PyMultivectorObject *data, Py_ssize_t size){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector_inherit_type(data); // pass -1 to inherit type
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        blades_array[i] = *((BladesMultivector*)data[i].data);

    *blades = atomic_blades_add_(blades_array,size,data->GA);

    if(blades->size == -1){
        PyMem_RawFree(blades);
        PyMem_RawFree(blades_array);
        multivector_dealloc(out);
        return NULL;
    }

    out->data = blades;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(blades_array);
    return out;
}



static PyMultivectorObject *atomic_mixed_product(PyMultivectorObject *data, Py_ssize_t size,PyMultivectorObject *def, ProductType ptype){
    PyMultivectorIter *iter = init_multivector_iter(data,size);
    SparseMultivector *psparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(def->GA,"sparselarge");
    if(!iter || !psparse || !out || !def){
        PyMem_RawFree(psparse);
        multivector_dealloc(out);
        free_multivector_iter(iter,size);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = atomic_mixed_geometricproduct_(iter,size,def->GA);
            break;
        case ProductType_inner:
            *psparse = atomic_mixed_innerproduct_(iter,size,def->GA);
            break;
        case ProductType_outer:
            *psparse = atomic_mixed_outerproduct_(iter,size,def->GA);
            break;
        default:
            PyMem_RawFree(psparse);
            multivector_dealloc(out);
            free_multivector_iter(iter,size);
            return NULL;
    }

    free_multivector_iter(iter,size);
    out->data = (void*)psparse;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_mixed_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *def, ProductType ptype){
    PyMultivectorIter *iter0 = init_multivector_iter(data0,1);
    PyMultivectorIter *iter1 = init_multivector_iter(data1,1);
    SparseMultivector *psparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(def->GA,"sparselarge");
    if(!iter1 || !iter0 || !psparse || !out || !def){
        free_multivector_iter(iter1,1);
        free_multivector_iter(iter0,1);
        PyMem_RawFree(psparse);
        multivector_dealloc(out);
        return NULL;
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = binary_mixed_geometricproduct_(iter0,iter1,def->GA);
            break;
        case ProductType_inner:
            *psparse = binary_mixed_innerproduct_(iter0,iter1,def->GA);
            break;
        case ProductType_outer:
            *psparse = binary_mixed_outerproduct_(iter0,iter1,def->GA);
            break;
        case ProductType_regressive:
            *psparse = binary_mixed_regressiveproduct_(iter0,iter1,def->GA);
            break;
        default:
            PyMem_RawFree(psparse);
            multivector_dealloc(out);
            free_multivector_iter(iter1,1);
            free_multivector_iter(iter0,1);
            return NULL;
    }

    free_multivector_iter(iter0,1);
    free_multivector_iter(iter1,1);

    if(psparse->size == -1){
        PyMem_RawFree(psparse);
        multivector_dealloc(out);
        return NULL;
    }
    out->data = (void*)psparse;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static void* blades_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    *blades = blades_init_(bitmap,value,size,ga);
    return (void*)blades;
}


PyMultivectorMixedMath_Funcs largemultivector_mixed_fn = {
  .add = NULL,
  .product = binary_mixed_product,
  .atomic_add = NULL,
  .atomic_product = atomic_mixed_product,
  .type_names = {"sparselarge","denselarge","bladeslarge",NULL},
};


 PyMultivectorMath_Funcs largemultivector_sparse_math_fn = {
    .product = binary_sparse_product,
    .atomic_product = atomic_sparse_product,
    .ternary_product = ternary_sparse_product,
    .grade_project = unary_sparse_gradeproject,
    .reverse = unary_sparse_reverse,
    .add = NULL,
    .atomic_add = NULL,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_sparse_dual,
    .undual = unary_sparse_undual,
};
 PyMultivectorMath_Funcs largemultivector_dense_math_fn = {
    .product = binary_dense_product,
    .atomic_product = atomic_dense_product,
    .ternary_product = ternary_dense_product,
    .grade_project = unary_dense_gradeproject,
    .reverse = unary_dense_reverse,
    .add = NULL,
    .atomic_add = NULL,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_dense_dual,
    .undual = unary_dense_undual,
};

PyMultivectorMath_Funcs largemultivector_blades_math_fn = {
    .product = binary_blades_product,
    .atomic_product = atomic_blades_product,
    .ternary_product = ternary_blades_product,
    .grade_project = NULL,
    .reverse = NULL,
    .add = binary_blades_add,
    .atomic_add = atomic_blades_add,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_blades_dual,
    .undual = unary_blades_undual,
};


 
PyMultivectorData_Funcs largemultivector_sparse_data_fn = {
    .free = NULL,
    .init = NULL,
    .iter_next = NULL,
    .iter_init = NULL,
};
 
PyMultivectorData_Funcs largemultivector_dense_data_fn = {
    .free = NULL,
    .init = NULL,
    .iter_next = NULL,
    .iter_init = NULL,
};

PyMultivectorData_Funcs largemultivector_blades_data_fn = {
    .free = NULL,
    .init = blades_init,
    .cast = cast_to_blades,
    .iter_next = NULL,
    .iter_init = NULL,
};


 const PyMultivectorSubType largesparse_subtype = {
    .math_funcs = &largemultivector_sparse_math_fn,
    .data_funcs = &largemultivector_sparse_data_fn,
    .name = "",
    .type_name = "sparselarge",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_sparse,
    .basic_size = sizeof(SparseMultivector),
};
 const PyMultivectorSubType largedense_subtype = {
    .math_funcs = &largemultivector_dense_math_fn,
    .data_funcs = &largemultivector_dense_data_fn,
    .name = "",
    .type_name = "denselarge",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_dense,
    .basic_size = sizeof(DenseMultivector),
};
 const PyMultivectorSubType largeblades_subtype = {
    .math_funcs = &largemultivector_blades_math_fn,
    .data_funcs = &largemultivector_blades_data_fn,
    .name = "",
    .type_name = "bladeslarge",
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_blades,
    .basic_size = sizeof(BladesMultivector),
};

PyMultivectorSubType largemultivector_subtypes_array[3] = {largesparse_subtype,largedense_subtype,largeblades_subtype};


void fill_missing_funcs(void){
     for(Py_ssize_t i = 0; i < 3; i++){
         if(largemultivector_subtypes_array[i].data_funcs->free == NULL)
            largemultivector_subtypes_array[i].data_funcs->free = multivector_subtypes_array[i].data_funcs->free;
         if(largemultivector_subtypes_array[i].data_funcs->init == NULL)
            largemultivector_subtypes_array[i].data_funcs->init = multivector_subtypes_array[i].data_funcs->init;
         if(largemultivector_subtypes_array[i].data_funcs->iter_next == NULL)
            largemultivector_subtypes_array[i].data_funcs->iter_next = multivector_subtypes_array[i].data_funcs->iter_next;
         if(largemultivector_subtypes_array[i].data_funcs->iter_init == NULL)
            largemultivector_subtypes_array[i].data_funcs->iter_init = multivector_subtypes_array[i].data_funcs->iter_init;
         if(largemultivector_subtypes_array[i].data_funcs->cast == NULL)
            largemultivector_subtypes_array[i].data_funcs->cast = multivector_subtypes_array[i].data_funcs->cast;
         if(largemultivector_subtypes_array[i].math_funcs->product == NULL)
            largemultivector_subtypes_array[i].math_funcs->product = multivector_subtypes_array[i].math_funcs->product;
         if(largemultivector_subtypes_array[i].math_funcs->atomic_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->atomic_product = multivector_subtypes_array[i].math_funcs->atomic_product;
         if(largemultivector_subtypes_array[i].math_funcs->ternary_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->ternary_product = multivector_subtypes_array[i].math_funcs->ternary_product;
         if(largemultivector_subtypes_array[i].math_funcs->grade_project == NULL)
            largemultivector_subtypes_array[i].math_funcs->grade_project = multivector_subtypes_array[i].math_funcs->grade_project;
         if(largemultivector_subtypes_array[i].math_funcs->reverse == NULL)
            largemultivector_subtypes_array[i].math_funcs->reverse = multivector_subtypes_array[i].math_funcs->reverse;
         if(largemultivector_subtypes_array[i].math_funcs->add == NULL)
            largemultivector_subtypes_array[i].math_funcs->add = multivector_subtypes_array[i].math_funcs->add;
         if(largemultivector_subtypes_array[i].math_funcs->atomic_add == NULL)
            largemultivector_subtypes_array[i].math_funcs->atomic_add = multivector_subtypes_array[i].math_funcs->atomic_add;
         if(largemultivector_subtypes_array[i].math_funcs->scalar_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->scalar_product = multivector_subtypes_array[i].math_funcs->scalar_product;
         if(largemultivector_subtypes_array[i].math_funcs->scalar_add == NULL)
            largemultivector_subtypes_array[i].math_funcs->scalar_add = multivector_subtypes_array[i].math_funcs->scalar_add;
         if(largemultivector_subtypes_array[i].math_funcs->dual == NULL)
            largemultivector_subtypes_array[i].math_funcs->dual = multivector_subtypes_array[i].math_funcs->dual;
         if(largemultivector_subtypes_array[i].math_funcs->undual == NULL)
            largemultivector_subtypes_array[i].math_funcs->undual = multivector_subtypes_array[i].math_funcs->undual;
     }

     largemultivector_mixed_fn.add = multivector_mixed_fn.add;
     largemultivector_mixed_fn.atomic_add = multivector_mixed_fn.atomic_add;
}
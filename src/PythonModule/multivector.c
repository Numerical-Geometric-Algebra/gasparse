#include <Python.h>
#include "gasparse.h"


static PyMultivectorObject *new_multivector(PyMultivectorObject *old, MultivectorType type){ // return a multivector of type type
    PyMultivectorObject *self = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    Py_SET_REFCNT((PyObject*)self,1);
    self->math_funcs = old->math_funcs;
    self->data_funcs = old->data_funcs;
    self->type = type;
    Py_SET_TYPE(self,Py_TYPE(old));
    Py_XINCREF(Py_TYPE(self));
    self->GA = old->GA;
    Py_XINCREF((PyObject*)self->GA);
    self->data = NULL;
    return self;
}

static void free_multivector(PyMultivectorObject *self){
    Py_XDECREF((PyObject*)self->GA);
    Py_XDECREF(Py_TYPE(self));
    PyMem_RawFree(self);
}

// returns true if abs(v) < p
static int comp_abs(ga_float v, ga_float p){
    ga_float r = (v < 0) ? -v : v;
    return r < p;
}

// returns true if abs(v) < ga->precision
static int ga_check_value_precision(PyAlgebraObject *ga, ga_float v){
    ga_float r = (v < 0) ? -v : v;
    return r < ga->precision;
}

static SparseMultivector init_sparse_empty(Py_ssize_t size){
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
    for(Py_ssize_t i = 0; i < size; i++){
        sparse.bitmap[i] = -1;
        sparse.value[i] = 0;
    }

    return sparse;
}

static DenseMultivector init_dense_empty(Py_ssize_t size){
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


static BladesMultivector init_blades_empty(Py_ssize_t size){
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

static void sparse_remove_small(SparseMultivector y, ga_float precision, Py_ssize_t *size){
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


static SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense, Py_ssize_t size){
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
        grade = gm.grade[bitmap]; gsize[grade]--;
        sparse.data[gindex[grade]].bitmap[gsize[grade]] = bitmap;
        sparse.data[gindex[grade]].value[gsize[grade]] = dense.value[i];
    }

    PyMem_RawFree(gsize);
    PyMem_RawFree(gindex);
    return sparse;
}

static Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades){
    Py_ssize_t *g = (Py_ssize_t*)PyMem_RawMalloc(n_grades*sizeof(Py_ssize_t));
    if(!g){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory");
        return NULL;
    }
    if(size == 0){ // if size is 0 project to all grades
        for(Py_ssize_t i = 0; i < n_grades; i++)
            g[i] = 1;
    }else{
        for(Py_ssize_t i = 0; i < n_grades; i++)
            g[i] = 0;
        for(Py_ssize_t i = 0; i < size; i++)
            g[grades[i]] = 1;
    }
    return g;
}


static void sparse_free_(SparseMultivector sparse){
    PyMem_RawFree(sparse.value);
    PyMem_RawFree(sparse.bitmap);
}

static void blades_free_(BladesMultivector blades){
    if(blades.data){
        for(Py_ssize_t i = 0; i < blades.size; i++){
            PyMem_RawFree(blades.data[i].bitmap);
            PyMem_RawFree(blades.data[i].value);
        }
        PyMem_RawFree(blades.data);
    }
    PyMem_RawFree(blades.grade);
}

static void dense_free_(DenseMultivector dense){
    PyMem_RawFree(dense.value);
}

static SparseMultivector sparse_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    SparseMultivector sparse;
    sparse.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    sparse.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    sparse.size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        sparse.value[i] = value[i];
        sparse.bitmap[i] = bitmap[i];
    }
    return sparse;
}

static BladesMultivector blades_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    GradeMap gm = ga->gm;
    SparseMultivector ssparse = {.bitmap = bitmap, .value = value, .size = size};
    return sparse_dense_to_blades_sparse(ssparse,gm);
}

static DenseMultivector dense_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    Py_ssize_t algebra_size = ga->product[ProductType_geometric].size;
    DenseMultivector dense = {.value = NULL, .size = -1};
    dense.value = (ga_float*)PyMem_RawMalloc(algebra_size*sizeof(ga_float));
    dense.size = algebra_size;
    // set all values to 0
    for(Py_ssize_t i = 0; i < algebra_size; i++)
        dense.value[i] = 0;
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= algebra_size){
            PyMem_RawFree(dense.value);
            dense.value = NULL;
            dense.size = -1;
            return dense; // raise error
        }
        dense.value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return dense;
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


static Py_ssize_t get_multivector_size(PyMultivectorObject *data){
    switch(data->type){
        case MultivectorType_sparse:;
            SparseMultivector *sparse = (SparseMultivector*)data->data;
            return sparse->size;
            break;
        case MultivectorType_dense:;
            DenseMultivector *dense = (DenseMultivector*)data->data;
            return dense->size;
        case MultivectorType_blades:;
            BladesMultivector *blades = (BladesMultivector*)data->data;
            Py_ssize_t size = 0;
            for(Py_ssize_t i = 0; i < blades->size; i++)
                size += blades->data[i].size;
            return size;
        default:
            return -1;
    }
}

static void free_multivector_iter(PyMultivectorIter *iter, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++)
        if(iter[i].index)
            free(iter[i].index);
    free(iter);
}


static PyMultivectorIter *init_multivector_iter(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorIter *iter = (PyMultivectorIter*)PyMem_RawMalloc(size*sizeof(PyMultivectorIter));
    for(Py_ssize_t i = 0; i < size; i++){
        /* iter[i].next = {sparse_iter_next,dense_iter_next,blades_iter_next,NULL}; */
        iter[i].data = data[i].data;
        iter[i].bitmap = -1;
        iter[i].value = 0;
        iter[i].grade = -1;
        iter[i].type = data[i].type;
        if(data[i].type == MultivectorType_blades){
            iter[i].index = (Py_ssize_t*)PyMem_RawMalloc(2*sizeof(Py_ssize_t));
            iter[i].index[0] = 0;
            iter[i].index[1] = 0;
            iter[i].size = 2;
        }else{
            iter[i].index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
            iter[i].index[0] = 0;
            iter[i].size = 1;
        }
        switch(data[i].type){
            case MultivectorType_sparse:
                iter[i].next = sparse_iter_next;
                break;
            case MultivectorType_dense:
                iter[i].next = dense_iter_next;
                break;
            case MultivectorType_blades:
                iter[i].next = blades_iter_next;
                break;
            default:
                iter[i].next = NULL;
                break;
        }
    }

    return iter;
}


static char *bitmap_to_string(int bitmap){
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

static PyObject *type_iter_repr(PyMultivectorIter *iter, PrintTypeMV ptype, Py_ssize_t dsize){
    if(ptype == PrintTypeMV_reduced){
        if(dsize){
            char **str_blade = (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            Py_ssize_t len = 0;
            PyObject *out = NULL;
            char *out_str;
            char sep[] = " + ";

            Py_ssize_t i = 0;
            while(iter->next(iter)){
                char *value = PyOS_double_to_string((double)iter->value,'f',6,0,NULL);
                if(iter->bitmap){
                    char *bitmap = bitmap_to_string(iter->bitmap);
                    Py_ssize_t size = strlen(bitmap) + strlen(value)+3;
                    str_blade[i] = (char*)PyMem_RawMalloc(size*sizeof(char));
                    PyOS_snprintf(str_blade[i],size,"%s^%s",value,bitmap);
                    PyMem_RawFree(bitmap);
                    PyMem_Free(value);
                }else{
                    Py_ssize_t size = strlen(value) + 1;
                    str_blade[i] = (char*)PyMem_RawMalloc(size*sizeof(char));
                    strcpy(str_blade[i],value);
                    PyMem_Free(value);
                }
                len += strlen(str_blade[i]);
                i++;
            }

            len += dsize*3 + 3;
            out_str = (char*)PyMem_RawMalloc(len*sizeof(char));

            Py_ssize_t j = 0;
            for(i = 0; i < dsize-1; i++){
                strcpy(out_str + j, str_blade[i]);
                j += strlen(str_blade[i]);
                strcpy(out_str + j, sep);
                j += strlen(sep);
            }
            strcpy(out_str + j, str_blade[i]);
            j += strlen(str_blade[i]);

            out = Py_BuildValue("s",out_str);

            for(Py_ssize_t i = 0; i < dsize; i++)
                PyMem_RawFree(str_blade[i]);

            PyMem_RawFree(out_str);

            return out;
        }else return Py_BuildValue("s","0.0");
    }else if(ptype == PrintTypeMV_normal){
        if(dsize){
            char **str_bitmap = (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            char **str_value =  (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            Py_ssize_t len_bitmap = 0, len_value = 0;
            PyObject *out = NULL;
            Py_ssize_t i = 0;
            while(iter->next(iter)){
                str_bitmap[i] = bitmap_to_string(iter->bitmap);
                len_bitmap += strlen(str_bitmap[i]);
                str_value[i] = PyOS_double_to_string((double)iter->value,'f',6,0,NULL);
                len_value += strlen(str_value[i]);
                i++;
            }
            len_value += dsize + 3;
            len_bitmap += dsize + 3;
            char *format_value = (char*)PyMem_RawMalloc((len_value)*sizeof(char));
            char *format_bitmap = (char*)PyMem_RawMalloc((len_bitmap)*sizeof(char));

            Py_ssize_t j = 0, k = 0; i = 0;
            for(i = 0; i < dsize - 1; i++){
                strcpy(format_bitmap + j, str_bitmap[i]);
                j += strlen(str_bitmap[i]);
                format_bitmap[j++] = ',';

                strcpy(format_value + k, str_value[i]);
                k += strlen(str_value[i]);
                format_value[k++] = ',';
            }
            strcpy(format_bitmap + j, str_bitmap[i]);
            j += strlen(str_bitmap[i]);
            format_bitmap[j] = '\0';


            strcpy(format_value + k, str_value[i]);
            k += strlen(str_value[i]);
            format_value[k] = '\0';

            char type_str[10];
            switch(iter->type){
                case MultivectorType_sparse:
                    strcpy(type_str,"sparse");
                    break;
                case MultivectorType_blades:
                    strcpy(type_str,"blades");
                    break;
                case MultivectorType_dense:
                    strcpy(type_str,"dense");
                    break;
                default:
                    strcpy(type_str,"???");
                    break;
            }

            char *format = ".multivector([%s],blades=[%s],dtype=%s)";
            size_t size = len_value + len_bitmap + strlen(format);
            char *format_out = (char*)PyMem_RawMalloc(size*sizeof(char));

            PyOS_snprintf(format_out,size,format,format_value,format_bitmap,type_str);
            out = Py_BuildValue("s",format_out);

            for(Py_ssize_t i = 0; i < dsize; i++){
                PyMem_RawFree(str_bitmap[i]);
                PyMem_Free(str_value[i]);
            }
            PyMem_RawFree(str_bitmap);
            PyMem_RawFree(str_value);
            PyMem_RawFree(format_value);
            PyMem_RawFree(format_bitmap);
            PyMem_RawFree(format_out);
            return out;
        }else return Py_BuildValue("s",".multivector(0.0)");
    }
    return NULL; // raise error
}


static SparseMultivector unary_sparse_scalaradd_(SparseMultivector sparse0, PyAlgebraObject *ga, ga_float value, int sign){
    SparseMultivector sparse = {.size = -1};
    Py_ssize_t scalarindex = -1;
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        if(sparse0.bitmap[i] == 0){
            scalarindex = i;
            break;
        }
    }


    if((scalarindex != -1 && ga_check_value_precision(ga,sign*sparse0.value[scalarindex] + value)) ||
       (scalarindex == -1 && ga_check_value_precision(ga,value))){
        Py_ssize_t size = sparse0.size;
        if(scalarindex != -1)
            size -= 1;

        sparse = init_sparse_empty(size);
        if(sparse.size == -1)
            return sparse;

        Py_ssize_t j = 0;
        for(Py_ssize_t i = 0; i < sparse0.size; i++){
            if(i == scalarindex) i++; // skip scalar index
            sparse.value[j] = sign*sparse0.value[i];
            sparse.bitmap[j] = sparse0.bitmap[i];
            j++;
        }
        return  sparse;
    }


    if(scalarindex != -1){
        sparse = init_sparse_empty(sparse0.size);
        if(sparse.size == -1)
            return sparse;

        for(Py_ssize_t i = 0; i < sparse0.size; i++){
            sparse.value[i] = sign*sparse0.value[i];
            sparse.bitmap[i] = sparse0.bitmap[i];
        }
        sparse.value[scalarindex] += value;
    }
    else{
        sparse = init_sparse_empty(sparse0.size + 1);
        if(sparse.size == -1)
            return sparse;
        sparse.value[0] = value;
        sparse.bitmap[0] = 0;

        for(Py_ssize_t i = 0; i < sparse0.size; i++){
            sparse.value[i+1] = sign*sparse0.value[i];
            sparse.bitmap[i+1] = sparse0.bitmap[i];
        }
    }

    return sparse;
}

static SparseMultivector unary_sparse_scalarproduct_(SparseMultivector sparse0, PyAlgebraObject *ga, ga_float sign){
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1)
        return sparse;
    Py_ssize_t size = 0;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        if(dense.bitmap[sparse0.bitmap[i]] == -1){
            dense.bitmap[sparse0.bitmap[i]] = sparse0.bitmap[i];
            size++;
        }
        dense.value[sparse0.bitmap[i]] += sign*sparse0.value[i];
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_add_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga, int sign){
    Py_ssize_t size = 0;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        if(dense.bitmap[sparse0.bitmap[i]] == -1){
            dense.bitmap[sparse0.bitmap[i]] = sparse0.bitmap[i];
            size++;
        }
        dense.value[sparse0.bitmap[i]] += sparse0.value[i];
    }
    for(Py_ssize_t i = 0; i < sparse1.size; i++){
        if(dense.bitmap[sparse1.bitmap[i]] == -1){
            dense.bitmap[sparse1.bitmap[i]] = sparse1.bitmap[i];
            size++;
        }
        dense.value[sparse1.bitmap[i]] += sign*sparse1.value[i];
    }
    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_product_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector ternary_sparse_product_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = geometric_product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0.bitmap[i]][sparse1.bitmap[j]];
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

    // dense0 = geometric_product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[dense1.bitmap[i]][sparse2.bitmap[j]];
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
    Py_SET_REFCNT((PyObject*)&sparse,1);
    return sparse;
}

static SparseMultivector unary_sparse_gradeproject_(SparseMultivector sparse0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    SparseMultivector sparse = {.size = -1};
    GradeMap gm = ga->gm;
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g)
        return sparse;

    int size = 0;
    for(Py_ssize_t i = 0; i < sparse0.size; i++)
        if(g[gm.grade[sparse0.bitmap[i]]])
            size++;

    sparse = init_sparse_empty(size--);
    if(sparse.size == -1){
        PyMem_RawFree(g);
        return sparse;
    }

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        if(g[gm.grade[sparse0.bitmap[i]]]){
            sparse.value[size] = sparse0.value[i];
            sparse.bitmap[size] = sparse0.bitmap[i];
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    Py_SET_REFCNT((PyObject*)&sparse,1);
    return sparse;
}

static SparseMultivector unary_sparse_reverse_(SparseMultivector sparse0, PyAlgebraObject *ga){
    GradeMap gm = ga->gm;
    SparseMultivector sparse = init_sparse_empty(sparse0.size);
    if(sparse.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        int grade = gm.grade[sparse0.bitmap[i]];
        int sign = (grade & 2) ? -1 : 1;
        sparse.value[i] = sign*sparse0.value[i];
        sparse.bitmap[i] = sparse0.bitmap[i];
    }

    return sparse;
}



static BladesMultivector unary_blades_scalaradd_(BladesMultivector blades0, PyAlgebraObject *ga, ga_float value, int sign){
    BladesMultivector blades = {.size = -1};
    SparseMultivector scalar = {.size = -1};
    Py_ssize_t scalarindex = -1;
    for(Py_ssize_t i = 0; i < blades0.size; i++){
        if(blades0.grade[i] == 0){
            scalarindex = i;
            break;
        }
    }

    if((scalarindex != -1 && ga_check_value_precision(ga,sign*blades0.data[scalarindex].value[0] + value)) ||
        (scalarindex == -1 && ga_check_value_precision(ga,value))){
        Py_ssize_t size = blades0.size;
        if(scalarindex != -1)
            size -= 1;

        blades = init_blades_empty(size);
        if(blades.size == -1)
            return blades;

        Py_ssize_t k = 0;
        for(Py_ssize_t i = 0; i < blades0.size; i++){
            if(i == scalarindex) i++; // skip scalar index
            Py_ssize_t sizei = blades0.data[i].size;
            blades.data[k] = init_sparse_empty(sizei);
            blades.grade[k] = blades0.grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades.data[k].value[j] = sign*blades0.data[i].value[j];
                blades.data[k].bitmap[j] = blades0.data[i].bitmap[j];
            }
            k++;
        }
        return blades;
    }

    if(scalarindex != -1){
        blades = init_blades_empty(blades0.size);
        if(blades.size == -1)
            return blades;

        for(Py_ssize_t i = 0; i < blades0.size; i++){
            Py_ssize_t sizei = blades0.data[i].size;
            blades.data[i] = init_sparse_empty(sizei);
            blades.grade[i] = blades0.grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades.data[i].value[j] = sign*blades0.data[i].value[j];
                blades.data[i].bitmap[j] = blades0.data[i].bitmap[j];
            }
        }
        *blades.data[scalarindex].value += value;
    }else{

        scalar = init_sparse_empty(1);
        if(scalar.size == -1)
            return blades;
        blades = init_blades_empty(blades0.size + 1);
        if(blades.size == -1)
            return blades;

        *scalar.value = value;
        *scalar.bitmap = 0;
        *blades.data = scalar;

        for(Py_ssize_t i = 0; i < blades0.size; i++){
            Py_ssize_t sizei = blades0.data[i].size;
            blades.data[i+1] = init_sparse_empty(sizei);
            blades.grade[i+1] = blades0.grade[i];
            for(Py_ssize_t j = 0; j < sizei; j++){
                blades.data[i+1].value[j] = sign*blades0.data[i].value[j];
                blades.data[i+1].bitmap[j] = blades0.data[i].bitmap[j];
            }
        }
    }

    return blades;
}


static BladesMultivector unary_blades_scalarproduct_(BladesMultivector blades0, PyAlgebraObject *ga, ga_float value){
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
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
            dense.value[bitmap] += value*sub.value[j];
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga->gm);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}

static BladesMultivector binary_blades_add_(BladesMultivector blades0, BladesMultivector blades1,  PyAlgebraObject *ga, int sign){
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
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

    sparse = sparse_dense_to_blades_sparse(dense,ga->gm);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}

static BladesMultivector binary_blades_product_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    GradeMap gm = ga->gm;
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

    sparse = sparse_dense_to_blades_sparse(dense,gm);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    Py_SET_REFCNT((PyObject*)&sparse,1);

    return sparse;
}

static BladesMultivector ternary_blades_product_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    GradeMap gm = ga->gm;
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

    sparse = sparse_dense_to_blades_sparse(dense0,gm);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static BladesMultivector unary_blades_gradeproject_(BladesMultivector blades0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    GradeMap gm = ga->gm;
    BladesMultivector blades = {.size = -1};
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g)
        return blades;

    int size = 0;
    for(Py_ssize_t i = 0; i < blades0.size; i++)
        if(g[blades0.grade[i]])
            size++;

    blades = init_blades_empty(size--);
    if(blades.size == -1){
        PyMem_RawFree(g);
        return blades;
    }

    // copy the values of the selected grades
    for(Py_ssize_t i = 0; i < blades0.size; i++){
        int grade = blades0.grade[i];
        if(g[grade]){
            Py_ssize_t bsize = blades0.data[i].size;
            blades.data[size].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
            blades.data[size].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
            if(!blades.data[size].bitmap || !blades.data[size].value){
                blades_free_(blades);
                PyMem_RawFree(g);
                blades.size = -1;
                return blades;
            }
            blades.data[size].size = bsize;
            blades.grade[size] = grade;
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades.data[size].bitmap[j] = blades0.data[i].bitmap[j];
                blades.data[size].value[j] = blades0.data[i].value[j];
            }
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    return blades;
}



static BladesMultivector unary_blades_reverse_(BladesMultivector blades0, PyAlgebraObject *ga){
    BladesMultivector blades = init_blades_empty(blades0.size);
    if(blades.size == -1)
        return blades;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        int grade = blades0.grade[i];
        Py_ssize_t bsize = blades0.data[i].size;
        blades.data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades.data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades.data[i].bitmap || !blades.data[i].value){
            blades_free_(blades);
            blades.size = -1;
            return blades;
        }
        blades.data[i].size = bsize;
        blades.grade[i] = grade;

        if(grade & 2){
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades.data[i].bitmap[j] = blades0.data[i].bitmap[j];
                blades.data[i].value[j] = -blades0.data[i].value[j];
            }
        }else{
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades.data[i].bitmap[j] = blades0.data[i].bitmap[j];
                blades.data[i].value[j] = blades0.data[i].value[j];
            }
        }
    }

    return blades;
}

static DenseMultivector unary_dense_scalaradd_(DenseMultivector dense0, PyAlgebraObject *ga, ga_float value, int sign){
    DenseMultivector dense = init_dense_empty(dense0.size);
    if(dense.size == -1)
        return dense;
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        dense.value[i] = sign*dense0.value[i];

    dense.value[0] += value; // increment the scalar part

    return dense;
}


static DenseMultivector unary_dense_scalarproduct_(DenseMultivector dense0, PyAlgebraObject *ga, ga_float value){
    DenseMultivector dense = {.size = -1};
    dense = init_dense_empty(dense0.size);
    if(dense.size == -1) return dense;
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        dense.value[i] = value*dense0.value[i];

    return dense;
}

static DenseMultivector binary_dense_add_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga, int sign){
    DenseMultivector dense = {.size = -1};
    dense = init_dense_empty(dense0.size);
    if(dense.size == -1) return dense;
    for(Py_ssize_t i = 0; i < dense.size; i++)
        dense.value[i] += dense0.value[i] + sign*dense1.value[i];

    return dense;
}

static DenseMultivector binary_dense_product_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    DenseMultivector dense = {.size = -1};
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            dense.value[m.bitmap[i][j]] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_product_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            temp.value[m.bitmap[i][j]] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            dense.value[m.bitmap[i][j]] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);

    return dense;
}

static DenseMultivector unary_dense_gradeproject_(DenseMultivector dense0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    GradeMap gm = ga->gm;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    if(!g) return dense;
    dense = init_dense_empty(dense0.size);
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(g[gm.grade[i]])
            dense.value[i] = dense0.value[i];

    PyMem_RawFree(g);

    return dense;
}

static DenseMultivector unary_dense_reverse_(DenseMultivector dense0, PyAlgebraObject *ga){
    GradeMap gm = ga->gm;
    DenseMultivector dense = init_dense_empty(dense0.size);
    if(dense.size == -1)
        return dense;

    for(Py_ssize_t i = 0; i < dense0.size; i++){
        int sign = (gm.grade[i] & 2) ? -1 : 1;
        dense.value[i] = sign*dense0.value[i];
    }

    return dense;
}






   
static PyMultivectorObject *unary_sparse_scalaradd(PyMultivectorObject *data0
    , ga_float value, int sign){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_scalaradd_(
        *sparse0,
        data0->GA, value, sign);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_scalarproduct(PyMultivectorObject *data0
    , ga_float value){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_scalarproduct_(
        *sparse0,
        data0->GA, value);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_sparse_add(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , int sign){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;
    SparseMultivector *sparse1 = (SparseMultivector*)data1->data;

    *sparse_out = binary_sparse_add_(
        *sparse0,
        *sparse1,
        data0->GA, sign);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_sparse_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , ProductType type){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;
    SparseMultivector *sparse1 = (SparseMultivector*)data1->data;

    *sparse_out = binary_sparse_product_(
        *sparse0,
        *sparse1,
        data0->GA, type);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

 
static PyMultivectorObject *ternary_sparse_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , PyMultivectorObject *data2
    , ProductType type){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;
    SparseMultivector *sparse1 = (SparseMultivector*)data1->data;
    SparseMultivector *sparse2 = (SparseMultivector*)data2->data;

    *sparse_out = ternary_sparse_product_(
        *sparse0,
        *sparse1,
        *sparse2,
        data0->GA, type);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_gradeproject(PyMultivectorObject *data0
    , int *grades, Py_ssize_t size){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_gradeproject_(
        *sparse0,
        data0->GA, grades, size);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_sparse_reverse(PyMultivectorObject *data0
    ){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    SparseMultivector *sparse0 = (SparseMultivector*)data0->data;

    *sparse_out = unary_sparse_reverse_(
        *sparse0,
        data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

   
static PyMultivectorObject *unary_blades_scalaradd(PyMultivectorObject *data0
    , ga_float value, int sign){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_scalaradd_(
        *blades0,
        data0->GA, value, sign);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_blades_scalarproduct(PyMultivectorObject *data0
    , ga_float value){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_scalarproduct_(
        *blades0,
        data0->GA, value);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_blades_add(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , int sign){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;
    BladesMultivector *blades1 = (BladesMultivector*)data1->data;

    *blades_out = binary_blades_add_(
        *blades0,
        *blades1,
        data0->GA, sign);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_blades_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , ProductType type){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;
    BladesMultivector *blades1 = (BladesMultivector*)data1->data;

    *blades_out = binary_blades_product_(
        *blades0,
        *blades1,
        data0->GA, type);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

 
static PyMultivectorObject *ternary_blades_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , PyMultivectorObject *data2
    , ProductType type){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;
    BladesMultivector *blades1 = (BladesMultivector*)data1->data;
    BladesMultivector *blades2 = (BladesMultivector*)data2->data;

    *blades_out = ternary_blades_product_(
        *blades0,
        *blades1,
        *blades2,
        data0->GA, type);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_blades_gradeproject(PyMultivectorObject *data0
    , int *grades, Py_ssize_t size){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_gradeproject_(
        *blades0,
        data0->GA, grades, size);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_blades_reverse(PyMultivectorObject *data0
    ){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_blades);
    BladesMultivector *blades0 = (BladesMultivector*)data0->data;

    *blades_out = unary_blades_reverse_(
        *blades0,
        data0->GA);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        free_multivector(out);
        return NULL;
    }

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

   
static PyMultivectorObject *unary_dense_scalaradd(PyMultivectorObject *data0
    , ga_float value, int sign){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_scalaradd_(
        *dense0,
        data0->GA, value, sign);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_scalarproduct(PyMultivectorObject *data0
    , ga_float value){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_scalarproduct_(
        *dense0,
        data0->GA, value);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_dense_add(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , int sign){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;
    DenseMultivector *dense1 = (DenseMultivector*)data1->data;

    *dense_out = binary_dense_add_(
        *dense0,
        *dense1,
        data0->GA, sign);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *binary_dense_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , ProductType type){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;
    DenseMultivector *dense1 = (DenseMultivector*)data1->data;

    *dense_out = binary_dense_product_(
        *dense0,
        *dense1,
        data0->GA, type);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

 
static PyMultivectorObject *ternary_dense_product(PyMultivectorObject *data0
    , PyMultivectorObject *data1
    , PyMultivectorObject *data2
    , ProductType type){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;
    DenseMultivector *dense1 = (DenseMultivector*)data1->data;
    DenseMultivector *dense2 = (DenseMultivector*)data2->data;

    *dense_out = ternary_dense_product_(
        *dense0,
        *dense1,
        *dense2,
        data0->GA, type);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_gradeproject(PyMultivectorObject *data0
    , int *grades, Py_ssize_t size){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_gradeproject_(
        *dense0,
        data0->GA, grades, size);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

  
static PyMultivectorObject *unary_dense_reverse(PyMultivectorObject *data0
    ){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_dense);
    DenseMultivector *dense0 = (DenseMultivector*)data0->data;

    *dense_out = unary_dense_reverse_(
        *dense0,
        data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        free_multivector(out);
        return NULL;
    }

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static SparseMultivector atomic_sparse_add_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1) return dense;
    SparseMultivector sparse;
    Py_ssize_t size = 0;

    for(Py_ssize_t j = 0; j < dsize; j++){
        for(Py_ssize_t i = 0; i < data[j].size; i++){
            if(dense.bitmap[data[j].bitmap[i]] == -1){
                dense.bitmap[data[j].bitmap[i]] = data[j].bitmap[i];
                size++;
            }
            dense.value[data[j].bitmap[i]] += data[j].value[i];
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}



static SparseMultivector atomic_sparse_product_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];

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
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
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

    sparse = sparse_dense_to_blades_sparse(dense,ga->gm);
    sparse_free_(dense);
    return sparse;
}

static BladesMultivector atomic_blades_product_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga, ProductType ptype){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = ga->product[ptype];
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga->gm);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}

static DenseMultivector atomic_dense_add_(DenseMultivector *data, Py_ssize_t size,  PyAlgebraObject *ga){
    DenseMultivector dense = init_dense_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1) return dense;
    ga_float value;
    for(Py_ssize_t i = 0; i < dense.size; i++){
        value = 0;
        for(Py_ssize_t k = 0; k < size; k++)
            value += data[k].value[i];
        dense.value[i] += value;
    }

    return dense;
}


static DenseMultivector atomic_dense_product_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
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
    return temp;
}


  static PyMultivectorObject *atomic_sparse_add(PyMultivectorObject *data, Py_ssize_t size){
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        sparse_array[i] = *((SparseMultivector*)data[i].data);

    *sparse = atomic_sparse_add_(sparse_array,size,data->GA);

    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        PyMem_RawFree(sparse_array);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(sparse_array);
    return out;
}
 static PyMultivectorObject *atomic_sparse_product(PyMultivectorObject *data, Py_ssize_t size, ProductType type){
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    SparseMultivector *sparse_array = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        sparse_array[i] = *((SparseMultivector*)data[i].data);

    *sparse = atomic_sparse_product_(sparse_array,size,data->GA,type);

    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        PyMem_RawFree(sparse_array);
        free_multivector(out);
        return NULL;
    }

    out->data = sparse;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(sparse_array);
    return out;
}
  static PyMultivectorObject *atomic_blades_add(PyMultivectorObject *data, Py_ssize_t size){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_blades);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        blades_array[i] = *((BladesMultivector*)data[i].data);

    *blades = atomic_blades_add_(blades_array,size,data->GA);

    if(blades->size == -1){
        PyMem_RawFree(blades);
        PyMem_RawFree(blades_array);
        free_multivector(out);
        return NULL;
    }

    out->data = blades;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(blades_array);
    return out;
}
 static PyMultivectorObject *atomic_blades_product(PyMultivectorObject *data, Py_ssize_t size, ProductType type){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_blades);
    BladesMultivector *blades_array = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        blades_array[i] = *((BladesMultivector*)data[i].data);

    *blades = atomic_blades_product_(blades_array,size,data->GA,type);

    if(blades->size == -1){
        PyMem_RawFree(blades);
        PyMem_RawFree(blades_array);
        free_multivector(out);
        return NULL;
    }

    out->data = blades;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(blades_array);
    return out;
}
  static PyMultivectorObject *atomic_dense_add(PyMultivectorObject *data, Py_ssize_t size){
    DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_dense);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        dense_array[i] = *((DenseMultivector*)data[i].data);

    *dense = atomic_dense_add_(dense_array,size,data->GA);

    if(dense->size == -1){
        PyMem_RawFree(dense);
        PyMem_RawFree(dense_array);
        free_multivector(out);
        return NULL;
    }

    out->data = dense;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(dense_array);
    return out;
}
 static PyMultivectorObject *atomic_dense_product(PyMultivectorObject *data, Py_ssize_t size, ProductType type){
    DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_dense);
    DenseMultivector *dense_array = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        dense_array[i] = *((DenseMultivector*)data[i].data);

    *dense = atomic_dense_product_(dense_array,size,data->GA,type);

    if(dense->size == -1){
        PyMem_RawFree(dense);
        PyMem_RawFree(dense_array);
        free_multivector(out);
        return NULL;
    }

    out->data = dense;
    Py_SET_REFCNT((PyObject*)out,1);

    PyMem_RawFree(dense_array);
    return out;
}


static SparseMultivector binary_mixed_add_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, int sign){
    SparseMultivector sparse;
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1) return dense;

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

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_mixed_product_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = m.bitmap[iter0->bitmap][iter1->bitmap];
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


// mixed type binary operations

 
static PyMultivectorObject *binary_mixed_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    PyMultivectorIter *iter0 = init_multivector_iter(data0,1);
    if(!iter0) return NULL;
    PyMultivectorIter *iter1 = init_multivector_iter(data1,1);
    if(!iter1){
        free_multivector_iter(iter0,1);
        return NULL;
    }
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    if(!sparse){
        free_multivector_iter(iter1,1);
        free_multivector_iter(iter0,1);
        return NULL;
    }
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    if(!out){
        PyMem_RawFree(sparse);
        free_multivector_iter(iter1,1);
        free_multivector_iter(iter0,1);
        return NULL;
    }

    *sparse = binary_mixed_add_(iter0,iter1,data0->GA,sign);

    free_multivector_iter(iter0,1);
    free_multivector_iter(iter1,1);

    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        free_multivector(out);
        return NULL;
    }
    out->data = (void*)sparse;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

 
static PyMultivectorObject *binary_mixed_product(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType type){
    PyMultivectorIter *iter0 = init_multivector_iter(data0,1);
    if(!iter0) return NULL;
    PyMultivectorIter *iter1 = init_multivector_iter(data1,1);
    if(!iter1){
        free_multivector_iter(iter0,1);
        return NULL;
    }
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    if(!sparse){
        free_multivector_iter(iter1,1);
        free_multivector_iter(iter0,1);
        return NULL;
    }
    PyMultivectorObject *out = new_multivector(data0,MultivectorType_sparse);
    if(!out){
        PyMem_RawFree(sparse);
        free_multivector_iter(iter1,1);
        free_multivector_iter(iter0,1);
        return NULL;
    }

    *sparse = binary_mixed_product_(iter0,iter1,data0->GA,type);

    free_multivector_iter(iter0,1);
    free_multivector_iter(iter1,1);

    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        free_multivector(out);
        return NULL;
    }
    out->data = (void*)sparse;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}



static SparseMultivector atomic_mixed_add_(PyMultivectorIter *iter, Py_ssize_t dsize, PyAlgebraObject *ga){
    SparseMultivector dense = init_sparse_empty(ga->product[ProductType_geometric].size);
    if(dense.size == -1) return dense;
    SparseMultivector sparse;

    Py_ssize_t size = 0;
    for(Py_ssize_t j = 0; j < dsize; j++){
        while(iter->next(iter)){
            if(dense.bitmap[iter->bitmap] == -1){
                dense.bitmap[iter->bitmap] = iter->bitmap;
                size++;
            }
            dense.value[iter->bitmap] += iter->value;
        }iter++;
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}



static SparseMultivector atomic_mixed_product_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga, ProductType ptype){
    CliffordMap m = ga->product[ptype];
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

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}



// atomic mixed type operations

 
static PyMultivectorObject *atomic_mixed_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorIter *iter = init_multivector_iter(data,size);
    if(!iter) return NULL;
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    if(!sparse){
        free_multivector_iter(iter,size);
        return NULL;
    }
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    if(!out){
        free_multivector_iter(iter,size);
        PyMem_RawFree(sparse);
        return NULL;
    }

    *sparse = atomic_mixed_add_(iter,size,data->GA);

    free_multivector_iter(iter,size);
    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        free_multivector(out);
        return NULL;
    }
    out->data = (void*)sparse;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

 
static PyMultivectorObject *atomic_mixed_product(PyMultivectorObject *data, Py_ssize_t size, ProductType type){
    PyMultivectorIter *iter = init_multivector_iter(data,size);
    if(!iter) return NULL;
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    if(!sparse){
        free_multivector_iter(iter,size);
        return NULL;
    }
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    if(!out){
        free_multivector_iter(iter,size);
        PyMem_RawFree(sparse);
        return NULL;
    }

    *sparse = atomic_mixed_product_(iter,size,data->GA,type);

    free_multivector_iter(iter,size);
    if(sparse->size == -1){
        PyMem_RawFree(sparse);
        free_multivector(out);
        return NULL;
    }
    out->data = (void*)sparse;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}



 void sparse_free(void *sparse){
    if(!sparse)
        return;
    sparse_free_(*((SparseMultivector*)sparse));
}

static void* sparse_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    *sparse = sparse_init_(bitmap,value,size,ga);
    return (void*)sparse;
}

 void blades_free(void *blades){
    if(!blades)
        return;
    blades_free_(*((BladesMultivector*)blades));
}

static void* blades_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    *blades = blades_init_(bitmap,value,size,ga);
    return (void*)blades;
}

 void dense_free(void *dense){
    if(!dense)
        return;
    dense_free_(*((DenseMultivector*)dense));
}

static void* dense_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    *dense = dense_init_(bitmap,value,size,ga);
    return (void*)dense;
}


/*
  Method calls of PyMultivectorType:
    print(a) -> multivector_repr(self): prints the multivector
    a(2)     -> multivector_grade_project(self): grade projection
  Method calls for the number protocol:
    a*b -> multivector_geometric_product: geometric product
    a^b -> multivector_outer_product: wedge product
    a|b -> multivector_inner_product: inner product
    a+b -> multivector_add: multivector addition
    a-b -> multivector_subtract: multivector subtraction
    ~a  -> multivector_invert: reverse the order of the basis vector by changing the sign
*/



PyObject *multivector_repr(PyMultivectorObject *self){
    PyObject *mv_repr;
    PyObject *ga_repr;
    PrintTypeMV ptype = self->GA->print_type_mv;

    PyMultivectorIter *iter = init_multivector_iter(self,1);
    mv_repr = type_iter_repr(iter,ptype,get_multivector_size(self));
    free_multivector_iter(iter,1);

    if(ptype == PrintTypeMV_normal){
        ga_repr = PyObject_Repr((PyObject*)self->GA);
        return PySequence_Concat(ga_repr,mv_repr);
    }else if(ptype == PrintTypeMV_reduced){
        return mv_repr;
    }else {
        PyErr_SetString(PyExc_ValueError,"The selected print type is not valid");
        return NULL;
    }
}

static int get_scalar(PyObject *self, ga_float *value){
    if(PyFloat_Check(self)){
        *value = (ga_float)PyFloat_AsDouble(self);
        return 1;
    }
    if(PyLong_Check(self)){
        *value = (ga_float)PyLong_AsDouble(self);
        return 1;
    }
    return 0;
}

static PyMultivectorObject* generate_empty_multivector(PyMultivectorObject *self){
    PyMultivectorObject *out;
    switch(self->type){
        case MultivectorType_sparse:
            out = new_multivector(self,self->type);
            SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
            sparse->size = 0;
            sparse->value = NULL;
            sparse->bitmap = NULL;
            out->data = (void*)sparse;
            return out;
        case MultivectorType_blades:
            out = new_multivector(self,self->type);
            BladesMultivector *blades = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
            blades->size = 0;
            blades->grade = NULL;
            blades->data = NULL;
            out->data = (void*)blades;
            return out;
        case MultivectorType_dense:
            out = new_multivector(self,self->type);
            DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
            *dense = init_dense_empty(((DenseMultivector*)self->data)->size);
            Py_SET_REFCNT((PyObject*)dense,1);
            out->data = (void*)dense;
            return out;
        default:
            return NULL;
    }
    return NULL;
}



static PyObject *multivector_product(PyObject *left, PyObject *right, ProductType ptype){
    PyMultivectorObject *data0 = NULL;
    PyMultivectorObject *data1 = NULL;
    ga_float value = 0;
    gaprodfunc product;
    gascalarfunc scalar_product;

    if(get_scalar(right,&value)) // check if right is a scalar
        data0 = (PyMultivectorObject*)left;
    else if(get_scalar(left,&value)) // check if left is a scalar
        data0 = (PyMultivectorObject*)right;

    if(data0){
        // return 0 if inner product with scalar
        if(ptype == ProductType_inner)
            return (PyObject*)generate_empty_multivector(data0);
        // multiply by scalar
        scalar_product = data0->math_funcs.scalar_product[data0->type];
        if(scalar_product){
            return (PyObject*)scalar_product(data0,value);
        }else{
            PyErr_SetString(PyExc_NotImplementedError,"The scalar product for this types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(PyObject_TypeCheck(left,Py_TYPE(right))){
        data0 = (PyMultivectorObject*)left;
        data1 = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(data0->GA != data1->GA){
        PyErr_SetString(PyExc_TypeError,"operands must have been generated by the same GA class");
        return NULL;
    }

    if(data0->type == data1->type){
        product = data0->math_funcs.product[data0->type];
        if(product){
            return (PyObject*)product(data0,data1,ptype);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
            return NULL; // raise not implemented error
        }
    } else{
        product = data0->math_funcs.mixed_product;
        if(product){
            return (PyObject*)product(data0,data1,ptype);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }

    }

    return NULL;
}

static PyObject *multivector_add_subtract(PyObject *left, PyObject *right, int sign){
    PyMultivectorObject *data0 = NULL;
    PyMultivectorObject *data1 = NULL;
    ga_float value = 0;
    gaaddfunc add;
    gascalaraddfunc scalar_add;

    if(get_scalar(right,&value)){ // check if right is a scalar
        data0 = (PyMultivectorObject*)left;
        value *= sign; // multiply by sign
        sign = 1;
    }else if(get_scalar(left,&value)){ // check if left is a scalar
        data0 = (PyMultivectorObject*)right;
    }

    if(data0){
        // add a scalar
        scalar_add = data0->math_funcs.scalar_add[data0->type];
        if(scalar_add){
            return (PyObject*)scalar_add(data0,value,sign);
        }else{
            PyErr_SetString(PyExc_NotImplementedError,"The scalar product for this types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(PyObject_TypeCheck(left,Py_TYPE(right))){
        data0 = (PyMultivectorObject*)left;
        data1 = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(data0->GA != data1->GA){
        PyErr_SetString(PyExc_TypeError,"operands must have been generated by the same GA class");
        return NULL;
    }

    if(data0->type == data1->type){
        add = data0->math_funcs.add[data0->type];
        if(add){
            return (PyObject*)add(data0,data1,sign);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
            return NULL; // raise not implemented error
        }
    } else{
        add = data0->math_funcs.mixed_add;
        if(add){
            return (PyObject*)add(data0,data1,sign);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    return NULL;
}



PyObject *multivector_outer_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_outer);
}

PyObject *multivector_inner_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_inner);
}

PyObject *multivector_geometric_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_geometric);
}

PyObject *multivector_grade_project(PyMultivectorObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"grades",NULL};
    int *grades = NULL;
    PyObject *grades_obj = NULL;
    Py_ssize_t size = -1;
    PyMultivectorObject *out = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &grades_obj))
        return NULL; // raise error

    if(PyLong_Check(grades_obj)){ // check if object is an integer
        int grade = (int)PyLong_AsLong(grades_obj);
        if(grade > self->GA->gm.max_grade)
            return NULL; // raise error
        grades = (int*)PyMem_RawMalloc(sizeof(int));
        *grades = grade;
        size = 1;
    }else if(PyList_Check(grades_obj)){ // check if object is a list type
        size = PyList_Size(grades_obj);
        grades = (int*)PyMem_RawMalloc(size*sizeof(int));
        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *grade_obj = PyList_GetItem(grades_obj,i);
            if(!PyLong_Check(grade_obj))
                return NULL; // raise error
            grades[i] = (int)PyLong_AsLong(grade_obj);
            if(grades[i] > self->GA->gm.max_grade){
                PyMem_RawFree(grades);
                return NULL; // raise error
            }
        }
    }

    gaunarygradefunc grade_project = self->math_funcs.grade_project[self->type];
    if(grade_project){
        out = grade_project(self,grades,size);
    }else{
        return NULL; // raise error
    }

    PyMem_RawFree(grades);
    return (PyObject*)out;
}

PyObject *multivector_add(PyObject *left, PyObject *right){
    return (PyObject*)multivector_add_subtract(left,right,1);
}

PyObject *multivector_subtract(PyObject *left, PyObject *right){
    return (PyObject*)multivector_add_subtract(left,right,-1);
}

PyObject *multivector_invert(PyMultivectorObject *self){
    gaunaryfunc reverse = self->math_funcs.reverse[self->type];
    PyMultivectorObject *out;
    if(reverse){
        out = reverse(self);
    }else{
        return NULL; // raise error
    }
    return (PyObject*)out;
}

static PyObject *multivector_sign(PyMultivectorObject *self, ga_float value){
    gascalarfunc scalar_product = self->math_funcs.scalar_product[self->type];
    PyMultivectorObject *out;
    if(scalar_product){
        out = scalar_product(self,value);
    }else{
        return NULL; // raise error
    }
    return (PyObject*)out;
}

PyObject *multivector_negative(PyMultivectorObject *self){
    return multivector_sign(self,-1);
}

PyObject *multivector_positive(PyMultivectorObject *self){
    return multivector_sign(self,1);
}

PyObject* multivector_atomic_add(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    MultivectorType mtype = -1;
    gaatomicfunc add = NULL;
    PyMultivectorObject *data0, *data1;
    PyMultivectorObject *data_array;
    int mixed_type = 0;

    if(size <= 1){
        PyErr_SetString(PyExc_ValueError,"number of arguments must be at least two");
        return NULL;
    }
    // check if objects are multivectors
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *argi = PyTuple_GetItem(args,i);
        if(!PyObject_IsInstance(argi,cls)){
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.multivector");
            return NULL;
        }
        data0 = (PyMultivectorObject*)argi;
        // determine the type of the multivectors
        if(mtype != -1){
            if(mtype != data0->type){
                mixed_type = 1;
            }
        } else {
            mtype = data0->type;
        }
    }

    if(!mixed_type){
        if(size == 2){
            gaaddfunc binary_add;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            binary_add = data0->math_funcs.add[mtype];
            if(binary_add){
                return (PyObject*)binary_add(data0,data1,1);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The binary sum operation for these types is not implemented");
                return NULL;
            }
        }
    }

    data_array = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++)
        data_array[i] = *((PyMultivectorObject*)PyTuple_GetItem(args,i));

    if(mixed_type)
        add = data_array->math_funcs.mixed_atomic_add;
    else
        add = data_array->math_funcs.atomic_add[mtype];

    if(add){
        return (PyObject*)add(data_array,size);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The atomic sum operation for these types is not implemented");
        return NULL; // raise not implemented error
    }

    return NULL;
}




static PyObject* multivector_atomic_product(PyObject *cls, PyObject *args, ProductType ptype){
    Py_ssize_t size = PyTuple_Size(args);
    MultivectorType mtype = -1;
    gaatomicprodfunc product = NULL;
    PyMultivectorObject *data0, *data1, *data2;
    PyMultivectorObject *data_array;
    int mixed_type = 0;

    if(size <= 1){
        PyErr_SetString(PyExc_ValueError,"number of arguments must be at least two");
        return NULL;
    }
    // check if objects are multivectors
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *argi = PyTuple_GetItem(args,i);
        if(!PyObject_IsInstance(argi,cls)){
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.multivector");
            return NULL;
        }
        data0 = (PyMultivectorObject*)argi;
        // determine the type of the multivectors
        if(mtype != -1){
            if(mtype != data0->type){
                mixed_type = 1;
            }
        } else {
            mtype = data0->type;
        }
    }

    if(!mixed_type){
        if(size == 2){
            gaprodfunc binary_product;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            binary_product = data0->math_funcs.product[mtype];
            if(binary_product){
                return (PyObject*)binary_product(data0,data1,ptype);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The binary product operation for these types is not implemented");
                return NULL;
            }
        }else if(size == 3){
            gaternaryprodfunc ternary_product;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            data2 = (PyMultivectorObject*)PyTuple_GetItem(args,2);

            ternary_product = data0->math_funcs.ternary_product[mtype];
            if(ternary_product){
                return (PyObject*)ternary_product(data0,data1,data2,ptype);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The ternary product operation for these types is not implemented");
                return NULL;
            }
        }
    }

    data_array = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++)
        data_array[i] = *((PyMultivectorObject*)PyTuple_GetItem(args,i));

    if(mixed_type)
        product = data_array->math_funcs.mixed_atomic_product;
    else
        product = data_array->math_funcs.atomic_product[mtype];

    if(product){
        return (PyObject*)product(data_array,size,ptype);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The atomic product operation for these types is not implemented");
        return NULL; // raise not implemented error
    }

    return NULL;
}


PyObject* multivector_atomic_geometric_product(PyObject *cls, PyObject *args){
    return multivector_atomic_product(cls,args,ProductType_geometric);
}

PyObject* multivector_atomic_outer_product(PyObject *cls, PyObject *args){
    return multivector_atomic_product(cls,args,ProductType_outer);
}

void multivector_dealloc(PyMultivectorObject *self){
    Py_XDECREF((PyObject*)self->GA);
    gafreefunc free_type = self->data_funcs.free[self->type];
    if(free_type)
        free_type(self->data);

    if(self->data) PyMem_RawFree(self->data);
    PyMem_RawFree(self);
}


static PyMultivectorMath_Funcs multivector_math_fn = {

    .product[MultivectorType_sparse] = binary_sparse_product,
    .product[MultivectorType_blades] = binary_blades_product,
    .product[MultivectorType_dense] = binary_dense_product,

    .atomic_product[MultivectorType_sparse] = atomic_sparse_product,
    .atomic_product[MultivectorType_blades] = atomic_blades_product,
    .atomic_product[MultivectorType_dense] = atomic_dense_product,

    .ternary_product[MultivectorType_sparse] = ternary_sparse_product,
    .ternary_product[MultivectorType_blades] = ternary_blades_product,
    .ternary_product[MultivectorType_dense] = ternary_dense_product,

    .add[MultivectorType_sparse] = binary_sparse_add,
    .add[MultivectorType_blades] = binary_blades_add,
    .add[MultivectorType_dense] = binary_dense_add,

    .atomic_add[MultivectorType_sparse] = atomic_sparse_add,
    .atomic_add[MultivectorType_blades] = atomic_blades_add,
    .atomic_add[MultivectorType_dense] = atomic_dense_add,

    .grade_project[MultivectorType_sparse] = unary_sparse_gradeproject,
    .grade_project[MultivectorType_blades] = unary_blades_gradeproject,
    .grade_project[MultivectorType_dense] = unary_dense_gradeproject,

    .scalar_product[MultivectorType_sparse] = unary_sparse_scalarproduct,
    .scalar_product[MultivectorType_blades] = unary_blades_scalarproduct,
    .scalar_product[MultivectorType_dense] = unary_dense_scalarproduct,

    .scalar_add[MultivectorType_sparse] = unary_sparse_scalaradd,
    .scalar_add[MultivectorType_blades] = unary_blades_scalaradd,
    .scalar_add[MultivectorType_dense] = unary_dense_scalaradd,

    .reverse[MultivectorType_sparse] = unary_sparse_reverse,
    .reverse[MultivectorType_blades] = unary_blades_reverse,
    .reverse[MultivectorType_dense] = unary_dense_reverse,

    .mixed_add = binary_mixed_add,
    .mixed_product = binary_mixed_product,

    .mixed_atomic_add = atomic_mixed_add,
    .mixed_atomic_product = atomic_mixed_product,


};

static PyMultivectorData_Funcs multivector_data_fn = {
    .free = {sparse_free,dense_free,blades_free},
    .init = {sparse_init,dense_init,blades_init},
};

PyMultivectorObject *init_multivector(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga, PyTypeObject *obj_type, MultivectorType type){

    PyMultivectorObject *self = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    self->math_funcs = multivector_math_fn;
    self->data_funcs = multivector_data_fn;

    if(type <= MultivectorTypeMIN || type >= MultivectorTypeMAX)
        return NULL; // raise error

    gainitfunc init = self->data_funcs.init[type];
    if(init)
        self->data = init(bitmap,value,size,ga);
    else
        return NULL; // raise not implemented error

    self->type = type;
    Py_SET_TYPE(self,obj_type);
    Py_XINCREF(obj_type);
    self->GA = ga;
    Py_XINCREF((PyObject*)ga);
    Py_SET_REFCNT((PyObject*)self->data,1);
    Py_SET_REFCNT((PyObject*)self,1);
    return self;
}
#include <Python.h>
#include "gasparse.h"

// returns true if abs(v) < p
static int comp_abs(ga_float v, ga_float p){
    ga_float r = (v < 0) ? -v : v;
    return r < p;
}


static Py_ssize_t *init_grade_size(GradeMap gm){
    Py_ssize_t *gsize = (Py_ssize_t*)PyMem_RawMalloc((gm.max_grade+1)*sizeof(Py_ssize_t));
    for(Py_ssize_t i = 0; i <= gm.max_grade; i++)
        gsize[i] = 0;
    return gsize;
}


static SparseMultivector init_sparse_empty(Py_ssize_t size){
    SparseMultivector y;
    y.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    y.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    y.size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        y.bitmap[i] = -1;
        y.value[i] = 0;
    }
    return y;
}


static DenseMultivector init_dense_empty(Py_ssize_t size){
    DenseMultivector y;
    y.value = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    y.size = size;
    for(Py_ssize_t i = 0; i < size; i++)
        y.value[i] = 0;

    return y;
}


static BladesMultivector init_blades_empty(Py_ssize_t size){
    BladesMultivector y;
    y.data = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    y.grade = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    y.size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        y.data[i].bitmap = NULL;
        y.data[i].value = NULL;
        y.grade[i] = i;
    }
    return y;
}

static int sparse_iter_next(PyMultivectorIter *iter){
    SparseMultivector *sparse = (SparseMultivector*)iter->data;
    if(*iter->index >= sparse->size)
        return 0;
    iter->bitmap = sparse->bitmap[*iter->index];
    iter->value = sparse->value[(*iter->index)++];
    iter->grade = grade(iter->bitmap);
    return 1;
}

static int dense_iter_next(PyMultivectorIter *iter){
    DenseMultivector *dense = (DenseMultivector*)iter->data;
    if(*iter->index >= dense->size)
        return 0;
    iter->bitmap = *iter->index;
    iter->value = dense->value[(*iter->index)++];
    iter->grade = grade(iter->bitmap);
    return 1;
}

static int blades_iter_next(PyMultivectorIter *iter){
    BladesMultivector *blades = (BladesMultivector*)iter->data;

    if(*iter->index >= blades->size)
        return 0;
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
        case MultivectorType_scalar:;
            return 1;
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
        }else if(data[i].type == MultivectorType_scalar){
            return NULL; //raise error
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

static SparseMultivector sparse_dense_to_sparse_sparse(SparseMultivector dense_y, Py_ssize_t size){
    Py_ssize_t k = 0;
    SparseMultivector sparse_y = init_sparse_empty(size);
    for(Py_ssize_t i = 0; i < dense_y.size; i++){
        if(k < size){
            if(dense_y.bitmap[i] != -1){
                sparse_y.bitmap[k] = dense_y.bitmap[i];
                sparse_y.value[k] = dense_y.value[i];
                k++;
            }
        }
    }
    sparse_y.size = size;
    return sparse_y;
}

static BladesMultivector sparse_dense_to_blades_sparse(SparseMultivector dense, GradeMap gm){
    BladesMultivector sparse;
    Py_ssize_t ssize = 0, grade = -1;
    Py_ssize_t *gsize = init_grade_size(gm);
    Py_ssize_t *gindex = init_grade_size(gm);
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

    return sparse;
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

static SparseMultivector sparse_init_(int *bitmap, ga_float *value, Py_ssize_t size){
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

static void *sparse_init(int *bitmap, ga_float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    *sparse = sparse_init_(bitmap,value,size);
    return (void*)sparse;
}

static BladesMultivector blades_init_(int *bitmap, ga_float *value, Py_ssize_t size, GradeMap gm){
    SparseMultivector ssparse = {.bitmap = bitmap, .value = value, .size = size};
    return sparse_dense_to_blades_sparse(ssparse,gm);
}

static void* blades_init(int *bitmap, ga_float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    BladesMultivector *sparse = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    *sparse = blades_init_(bitmap,value,size,gm);
    return (void*)sparse;
}

static DenseMultivector dense_init_(int *bitmap, ga_float *value, Py_ssize_t size, Py_ssize_t algebra_size){
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

static void *dense_init(int *bitmap, ga_float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    *dense = dense_init_(bitmap,value,size,algebra_size);
    if(!dense->value)
        return NULL;
    return (void*)dense;
}

static void sparse_free_(SparseMultivector sparse){
    PyMem_RawFree(sparse.value);
    PyMem_RawFree(sparse.bitmap);
}


void sparse_free(void *sparse){
    if(!sparse)
        return;
    sparse_free_(*((SparseMultivector*)sparse));
}

static void blades_free_(BladesMultivector blades){

    if(blades.data){
        for(Py_ssize_t i = 0; i < blades.size; i++){
            if(blades.data[i].bitmap != NULL)
                PyMem_RawFree(blades.data[i].bitmap);
            if(blades.data[i].value != NULL)
                PyMem_RawFree(blades.data[i].value);
        }
        PyMem_RawFree(blades.data);
    }
    if(blades.grade)
        PyMem_RawFree(blades.grade);
}

void blades_free(void *blades){
    if(!blades)
        return;
    blades_free_(*((BladesMultivector*)blades));
}

static void dense_free_(DenseMultivector dense){
    PyMem_RawFree(dense.value);
}

void dense_free(void *dense){
    if(!dense)
        return;
    dense_free_(*((DenseMultivector*)dense));
}


static char *bitmap_to_string(int bitmap){
    Py_ssize_t size = grade((Py_ssize_t)bitmap) + 2;
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


/*
** atomic product -> multivector.geometric_product(a,b,c,...)
**                           -> multivector.outer_product(a,b,c,...)
**
**     mixed_type_atomic_product(PyMultivectorObject[size],size,ProductType)
**     sparse_atomic_product(PyMultivectorObject[size],size,ProductType)
**     blades_atomic_product(PyMultivectorObject[size],size,ProductType)
**     dense_atomic_product(PyMultivectorObject[size],size,ProductType)
**
** atomic add     -> multivector.add(a,b,c,...)
**
**     mixed_type_atomic_add(PyMultivectorObject[size],size)
**     sparse_atomic_add(PyMultivectorObject[size],size)
**     blades_atomic_add(PyMultivectorObject[size],size)
**     dense_atomic_add(PyMultivectorObject[size],size)
*/

static SparseMultivector mixed_type_atomic_product_(PyMultivectorIter *iter, Py_ssize_t size, PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    SparseMultivector dense = init_sparse_empty(m.size);
    SparseMultivector temp = init_sparse_empty(m.size);
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

    sparse_remove_small(temp,ga.precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free((void*)&dense);
    sparse_free((void*)&temp);
    Py_SET_REFCNT((PyObject*)&sparse,1);
    return sparse;
}

static PyMultivectorObject *mixed_type_atomic_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    PyMultivectorIter *iter = init_multivector_iter(data,size);
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    *sparse_out = mixed_type_atomic_product_(iter,size,*data->GA,ptype);
    free_multivector_iter(iter,size);
    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static SparseMultivector mixed_type_atomic_add_(PyMultivectorIter *iter, Py_ssize_t dsize, PyGeometricAlgebraObject ga){
    SparseMultivector dense = init_sparse_empty(ga.product[ProductType_geometric].size);
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

    sparse_remove_small(dense,ga.precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);

    sparse_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}


static PyMultivectorObject *mixed_type_atomic_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorIter *iter = init_multivector_iter(data,size);
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    *sparse_out = mixed_type_atomic_add_(iter,size,*data->GA);
    free_multivector_iter(iter,size);
    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static SparseMultivector sparse_atomic_add_(SparseMultivector *data, Py_ssize_t dsize, PyGeometricAlgebraObject ga){
    SparseMultivector dense = init_sparse_empty(ga.product[ProductType_geometric].size);
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

    sparse_remove_small(dense,ga.precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);

    sparse_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}

static PyMultivectorObject *sparse_atomic_add(PyMultivectorObject *data, Py_ssize_t size){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    SparseMultivector *sparse_data = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        sparse_data[i] = *((SparseMultivector*)data[i].data);

    *sparse_out = sparse_atomic_add_(sparse_data,size,*data->GA);
    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static SparseMultivector sparse_atomic_product_(SparseMultivector *data, Py_ssize_t dsize, PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    SparseMultivector temp = init_sparse_empty(m.size);
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

    sparse_remove_small(temp,ga.precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free((void*)&dense);
    sparse_free((void*)&temp);
    Py_SET_REFCNT((PyObject*)&sparse,1);
    return sparse;
}


static PyMultivectorObject *sparse_atomic_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_sparse);
    SparseMultivector *sparse_data = (SparseMultivector*)PyMem_RawMalloc(size*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        sparse_data[i] = *((SparseMultivector*)data[i].data);

    *sparse_out = sparse_atomic_product_(sparse_data,size,*data->GA,ptype);
    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static BladesMultivector blades_atomic_product_(BladesMultivector *data, Py_ssize_t dsize, PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    SparseMultivector temp = init_sparse_empty(m.size);
    BladesMultivector sparse;
    Py_ssize_t tsize = 1;

    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

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

    sparse_remove_small(temp,ga.precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga.gm);
    sparse_free((void*)&dense);
    sparse_free((void*)&temp);
    Py_SET_REFCNT((PyObject*)&sparse,1);
    return sparse;
}


static PyMultivectorObject *blades_atomic_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_blades);
    BladesMultivector *blades_data = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        blades_data[i] = *((BladesMultivector*)data[i].data);

    *blades_out = blades_atomic_product_(blades_data,size,*data->GA,ptype);
    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static DenseMultivector dense_atomic_product_(DenseMultivector *data, Py_ssize_t dsize, PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory and set values to zero
    DenseMultivector dense = init_dense_empty(m.size);
    DenseMultivector temp = init_dense_empty(m.size);

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

    dense_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&temp,1);
    return temp;
}


static PyMultivectorObject *dense_atomic_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_dense);
    DenseMultivector *dense_data = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        dense_data[i] = *((DenseMultivector*)data[i].data);

    *dense_out = dense_atomic_product_(dense_data,size,*data->GA,ptype);
    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static BladesMultivector blades_atomic_add_(BladesMultivector *data, Py_ssize_t size, PyGeometricAlgebraObject ga){
    SparseMultivector dense = init_sparse_empty(ga.product[ProductType_geometric].size);
    BladesMultivector sparse;
    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga.precision;

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

    sparse = sparse_dense_to_blades_sparse(dense,ga.gm);

    sparse_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}

static PyMultivectorObject *blades_atomic_add(PyMultivectorObject *data, Py_ssize_t size){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_blades);
    BladesMultivector *blades_data = (BladesMultivector*)PyMem_RawMalloc(size*sizeof(BladesMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        blades_data[i] = *((BladesMultivector*)data[i].data);

    *blades_out = blades_atomic_add_(blades_data,size,*data->GA);
    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static DenseMultivector dense_atomic_add_(DenseMultivector *data, Py_ssize_t size,  PyGeometricAlgebraObject ga){
    DenseMultivector dense = init_dense_empty(ga.product[ProductType_geometric].size);
    ga_float value;
    for(Py_ssize_t i = 0; i < dense.size; i++){
        value = 0;
        for(Py_ssize_t k = 0; k < size; k++)
            value += data[k].value[i];
        dense.value[i] += value;
    }

    Py_SET_REFCNT((PyObject*)&dense,1); // initialize reference count to 1
    return dense;
}

static PyMultivectorObject *dense_atomic_add(PyMultivectorObject *data, Py_ssize_t size){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(data,MultivectorType_dense);
    DenseMultivector *dense_data = (DenseMultivector*)PyMem_RawMalloc(size*sizeof(DenseMultivector));
    for(Py_ssize_t i = 0; i < size; i++)
        dense_data[i] = *((DenseMultivector*)data[i].data);

    *dense_out = dense_atomic_add_(dense_data,size,*data->GA);
    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}



/*
** Multiplication by scalars operations -> a*1.234 or 1.234*a
**     sparse_scalar_product(SparseMultivector,ga_float)
**     blades_scalar_product(BladesMultivector,ga_float)
**     dense_scalar_product(DenseMultivector,ga_float)
*/

static SparseMultivector sparse_scalar_product_(SparseMultivector sparse, ga_float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
    SparseMultivector sparse_out;
    // the inner product with a scalar should be 0
    if(ptype == ProductType_inner){
        // empty multivector is better than 0 valued multivector
        sparse_out.bitmap = NULL;
        sparse_out.value = NULL;
        sparse_out.size = 0;
    }else{
        sparse_out = init_sparse_empty(sparse.size);
        for(Py_ssize_t i = 0; i < sparse.size; i++){
            sparse_out.value[i] = sparse.value[i]*scalar;
            sparse_out.bitmap[i] = sparse.bitmap[i];
        }
    }

    Py_SET_REFCNT((PyObject*)&sparse_out,1);
    return sparse_out;
}


static PyMultivectorObject *sparse_scalar_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_sparse);

    ga_float *scalar = NULL;
    SparseMultivector *sparse = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        sparse = (SparseMultivector*)left_mv->data;
        scalar = (ga_float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        sparse = (SparseMultivector*)right_mv->data;
        scalar = (ga_float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *sparse_out = sparse_scalar_product_(*sparse,*scalar,*left_mv->GA,ptype);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static BladesMultivector blades_scalar_product_(BladesMultivector blades, ga_float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
    BladesMultivector blades_out;
    // the inner product with a scalar should be 0
    if(ptype == ProductType_inner){
        // empty multivector is better than 0 valued multivector
        blades_out.grade = NULL;
        blades_out.data = NULL;
        blades_out.size = 0;
    }else{
        blades_out = init_blades_empty(blades.size);
        for(Py_ssize_t i = 0; i < blades.size; i++){
            blades_out.data[i] = sparse_scalar_product_(blades.data[i],scalar,ga,ptype);
        }
    }

    Py_SET_REFCNT((PyObject*)&blades_out,1);
    return blades_out;
}



static PyMultivectorObject *blades_scalar_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_blades);

    ga_float *scalar = NULL;
    BladesMultivector *blades = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        blades = (BladesMultivector*)left_mv->data;
        scalar = (ga_float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        blades = (BladesMultivector*)right_mv->data;
        scalar = (ga_float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *blades_out = blades_scalar_product_(*blades,*scalar,*left_mv->GA,ptype);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static DenseMultivector dense_scalar_product_(DenseMultivector dense, ga_float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
    DenseMultivector dense_out;
    // the inner product with a scalar should be 0
    if(ptype == ProductType_inner){
        // empty multivector is better than 0 valued multivector
        dense_out.value = NULL;
        dense_out.size = 0;
    }else{
        dense_out = init_dense_empty(dense.size);
        for(Py_ssize_t i = 0; i < dense.size; i++){
            dense_out.value[i] = dense.value[i]*scalar;
        }
    }

    Py_SET_REFCNT((PyObject*)&dense_out,1);
    return dense_out;
}

static PyMultivectorObject *dense_scalar_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_dense);

    ga_float *scalar = NULL;
    DenseMultivector *dense = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        dense = (DenseMultivector*)left_mv->data;
        scalar = (ga_float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        dense = (DenseMultivector*)right_mv->data;
        scalar = (ga_float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *dense_out = dense_scalar_product_(*dense,*scalar,*left_mv->GA,ptype);

    out->data = (void*)dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

/*
** Same type product -> a*b and a|b and a^b
**     sparse_sparse_product(SparseMultivector,SparseMultivector,ProductType)
**     blades_blades_product(BladesMultivector,BladesMultivector,ProductType)
**     dense_dense_product(DenseMultivector,DenseMultivector,ProductType)
*/

static SparseMultivector sparse_sparse_product_(SparseMultivector left, SparseMultivector right,PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    SparseMultivector sparse;
    Py_ssize_t size = 0;

    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < left.size; i++){
        for(Py_ssize_t j = 0; j < right.size; j++){
            sign = m.sign[left.bitmap[i]][right.bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[left.bitmap[i]][right.bitmap[j]];
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += left.value[i]*right.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga.precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free((void*)&dense);

    Py_SET_REFCNT((PyObject*)&sparse,1);

    return sparse;
}

static PyMultivectorObject *sparse_sparse_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_sparse);
    SparseMultivector *left = (SparseMultivector*)left_mv->data;
    SparseMultivector *right = (SparseMultivector*)right_mv->data;

    *sparse_out = sparse_sparse_product_(*left,*right,*left_mv->GA,ptype);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static BladesMultivector blades_blades_product_(BladesMultivector left, BladesMultivector right, PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];
    GradeMap gm = ga.gm;

    SparseMultivector dense = init_sparse_empty(m.size);
    BladesMultivector sparse;
    ga_float precision = ga.precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector sleft, sright;

    for(Py_ssize_t i = 0; i < left.size; i++){
        sleft = left.data[i];
        for(Py_ssize_t j = 0; j < right.size; j++){
            sright = right.data[j];
            for(Py_ssize_t k = 0; k < sright.size; k++){
                for(Py_ssize_t l = 0; l < sleft.size; l++){
                    sign = m.sign[sleft.bitmap[l]][sright.bitmap[k]];
                    if(!sign) continue;
                    bitmap = m.bitmap[sleft.bitmap[l]][sright.bitmap[k]];
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += sleft.value[l]*sright.value[k]*sign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,gm);

    sparse_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&sparse,1);

    return sparse;
}

static PyMultivectorObject *blades_blades_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_blades);
    BladesMultivector *left = (BladesMultivector*)left_mv->data;
    BladesMultivector *right = (BladesMultivector*)right_mv->data;

    *blades_out = blades_blades_product_(*left,*right,*left_mv->GA,ptype);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static DenseMultivector dense_dense_product_(DenseMultivector left, DenseMultivector right,PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory for a dense y
    DenseMultivector dense = init_dense_empty(m.size);
    int sign;

    for(Py_ssize_t i = 0; i < left.size; i++){
        for(Py_ssize_t j = 0; j < right.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            dense.value[m.bitmap[i][j]] += left.value[i]*right.value[j]*sign;
        }
    }

    Py_SET_REFCNT((PyObject*)&dense,1);
    return dense;
}

static PyMultivectorObject *dense_dense_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_dense);
    DenseMultivector *left = (DenseMultivector*)left_mv->data;
    DenseMultivector *right = (DenseMultivector*)right_mv->data;

    *dense_out = dense_dense_product_(*left,*right,*left_mv->GA,ptype);

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

/*
** Mixed type product -> a*b and a|b and a^b:
**    mixed_mixed_product(PyMultivectorObject,PyMultivectorObject,ProductType)
**
** Mixed type add -> a+b and a-b and 1.234 + a and 1.234 - a
**    mixed_mixed_add(PyMultivectorObject,PyMultivectorObject)
**    mixed_scalar_multiply(PyMultivectorObject,ga_float,Product) -> multiply by 1 or -1
*/

static PyMultivectorObject *mixed_scalar_multiply(PyMultivectorObject *self,ga_float scalar,ProductType ptype){
    PyMultivectorObject *out = NULL;
    switch(self->type){
        case MultivectorType_sparse:;
            SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
            out = new_multivector(self,MultivectorType_sparse);
            *sparse_out = sparse_scalar_product_(*((SparseMultivector*)self->data),scalar,*self->GA,ptype);
            out->data = sparse_out;
            break;
        case MultivectorType_blades:;
            BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
            out = new_multivector(self,MultivectorType_blades);
            *blades_out = blades_scalar_product_(*((BladesMultivector*)self->data),scalar,*self->GA,ptype);
            out->data = blades_out;
            break;
        case MultivectorType_dense:;
            DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
            out = new_multivector(self,MultivectorType_dense);
            *dense_out = dense_scalar_product_(*((DenseMultivector*)self->data),scalar,*self->GA,ptype);
            out->data = dense_out;
            break;
        default:
            return NULL; // raise error
    }
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static PyMultivectorObject *mixed_mixed_add(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, int sign){
    PyMultivectorObject *out;
    PyMultivectorObject *data = (PyMultivectorObject*)PyMem_RawMalloc(2*sizeof(PyMultivectorObject));
    PyMultivectorObject *right_mv_sign;

    right_mv_sign = mixed_scalar_multiply(right_mv,(ga_float)sign,ProductType_geometric);

    data[0] = *left_mv;
    data[1] = *right_mv_sign;

    out = mixed_type_atomic_add(data,2);

    Py_SET_REFCNT((PyObject*)out,1);
    PyMem_Free(data);
    Py_XDECREF((PyObject*)right_mv_sign);
    return out;
}

static PyMultivectorObject *mixed_mixed_product(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, ProductType ptype){
    PyMultivectorObject *out;
    PyMultivectorObject *data = (PyMultivectorObject*)PyMem_RawMalloc(2*sizeof(PyMultivectorObject));

    data[0] = *left_mv;
    data[1] = *right_mv;

    out = mixed_type_atomic_product(data,2,ptype);

    Py_SET_REFCNT((PyObject*)out,1);
    PyMem_Free(data);
    return out;
}


/*
** Same type add operations -> a + b and a - b
**
** sparse_sparse_add(SparseMultivector,SparseMultivector,int)
** blades_blades_add(BladesMultivector,BladesMultivector,int)
** dense_dense_add(DenseMultivector,DenseMultivector,int)
**
*/

static SparseMultivector sparse_sparse_add_(SparseMultivector left, SparseMultivector right,  PyGeometricAlgebraObject ga, int sign){
    SparseMultivector dense = init_sparse_empty(ga.product[ProductType_geometric].size);
    SparseMultivector sparse;
    Py_ssize_t size = 0;


    for(Py_ssize_t i = 0; i < left.size; i++){
        if(dense.bitmap[left.bitmap[i]]==-1){
            dense.bitmap[left.bitmap[i]] = left.bitmap[i];
            size++;
        }
        dense.value[left.bitmap[i]] += left.value[i];
    }
    for(Py_ssize_t i = 0; i < right.size; i++){
        if(dense.bitmap[right.bitmap[i]]==-1){
            dense.bitmap[right.bitmap[i]] = right.bitmap[i];
            size++;
        }
        dense.value[right.bitmap[i]] += sign*right.value[i];
    }
    sparse_remove_small(dense,ga.precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);

    sparse_free((void*)&dense);

    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}


static PyMultivectorObject *sparse_sparse_add(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, int sign){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_sparse);
    SparseMultivector *left = (SparseMultivector*)left_mv->data;
    SparseMultivector *right = (SparseMultivector*)right_mv->data;

    *sparse_out = sparse_sparse_add_(*left,*right,*left_mv->GA,sign);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static BladesMultivector blades_blades_add_(BladesMultivector left, BladesMultivector right,  PyGeometricAlgebraObject ga, int sign){
    SparseMultivector dense = init_sparse_empty(ga.product[ProductType_geometric].size);
    BladesMultivector sparse;
    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga.precision;

    for(Py_ssize_t i = 0; i < left.size; i++){
        sub = left.data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sub.value[j];
        }
    }

    for(Py_ssize_t i = 0; i < right.size; i++){
        sub = right.data[i];
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

    sparse = sparse_dense_to_blades_sparse(dense,ga.gm);

    sparse_free((void*)&dense);
    Py_SET_REFCNT((PyObject*)&sparse,1); // initialize reference count to 1
    return sparse;
}


static PyMultivectorObject *blades_blades_add(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, int sign){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_blades);
    BladesMultivector *left = (BladesMultivector*)left_mv->data;
    BladesMultivector *right = (BladesMultivector*)right_mv->data;

    *blades_out = blades_blades_add_(*left,*right,*left_mv->GA,sign);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


static DenseMultivector dense_dense_add_(DenseMultivector left, DenseMultivector right,  PyGeometricAlgebraObject ga, int sign){
    DenseMultivector dense = init_dense_empty(ga.product[ProductType_geometric].size);

    for(Py_ssize_t i = 0; i < dense.size; i++)
        dense.value[i] += left.value[i] + sign*right.value[i];

    Py_SET_REFCNT((PyObject*)&dense,1); // initialize reference count to 1
    return dense;
}


static PyMultivectorObject *dense_dense_add(PyMultivectorObject *left_mv, PyMultivectorObject *right_mv, int sign){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(left_mv,MultivectorType_dense);
    DenseMultivector *left = (DenseMultivector*)left_mv->data;
    DenseMultivector *right = (DenseMultivector*)right_mv->data;

    *dense_out = dense_dense_add_(*left,*right,*left_mv->GA,sign);

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}


/*
** Grade projection -> a([0,5,3,...]) and a(0)
**   sparse_grade_project(SparseMultivector,grades[size],size)
**   blades_grade_project(BladesMultivector,grades[size],size)
**   dense_grade_project(BladesMultivector,grades[size],size)
**   get_grade_bool(grades[size],size,max_grade) -> gets a bool array for each grade
**
** Reversion -> ~a
**   sparse_reverse(SparseMultivector)
**   blades_reverse(SparseMultivector)
**   dense_reverse(SparseMultivector)
**
*/

static Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades){
    Py_ssize_t *g = (Py_ssize_t*)PyMem_RawMalloc(n_grades*sizeof(Py_ssize_t));
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

static SparseMultivector sparse_grade_project_(SparseMultivector sparse, GradeMap gm, int *grades, Py_ssize_t grade_size){
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    SparseMultivector sparse_out;
    int size = 0;
    for(Py_ssize_t i = 0; i < sparse.size; i++)
        if(g[gm.grade[sparse.bitmap[i]]])
            size++;

    sparse_out = init_sparse_empty(size--);

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < sparse.size; i++){
        if(g[gm.grade[sparse.bitmap[i]]]){
            sparse_out.value[size] = sparse.value[i];
            sparse_out.bitmap[size] = sparse.bitmap[i];
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);

    Py_SET_REFCNT((PyObject*)&sparse_out,1);
    return sparse_out;
}


static PyMultivectorObject *sparse_grade_project(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_sparse);
    SparseMultivector *sparse = (SparseMultivector*)self->data;

    *sparse_out = sparse_grade_project_(*sparse,self->GA->gm,grades,size);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;

}


static BladesMultivector blades_grade_project_(BladesMultivector blades, GradeMap gm, int *grades, Py_ssize_t grade_size){
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    BladesMultivector blades_out;
    int size = 0;
    for(Py_ssize_t i = 0; i < blades.size; i++)
        if(g[blades.grade[i]])
            size++;

    blades_out = init_blades_empty(size--);

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < blades.size; i++){
        int grade = blades.grade[i];
        if(g[grade]){
            Py_ssize_t bsize = blades.data[i].size;
            blades_out.data[size].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
            blades_out.data[size].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
            blades_out.data[size].size = bsize;
            blades_out.grade[size] = grade;
            for(Py_ssize_t j = 0; j < bsize; j++){
                blades_out.data[size].bitmap[j] = blades.data[i].bitmap[j];
                blades_out.data[size].value[j] = blades.data[i].value[j];
            }
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);

    Py_SET_REFCNT((PyObject*)&blades_out,1);
    return blades_out;
}


static PyMultivectorObject *blades_grade_project(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_blades);
    BladesMultivector *blades = (BladesMultivector*)self->data;

    *blades_out = blades_grade_project_(*blades,self->GA->gm,grades,size);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;

}

static DenseMultivector dense_grade_project_(DenseMultivector dense, GradeMap gm, int *grades, Py_ssize_t grade_size){
    Py_ssize_t *g = get_grade_bool(grades,grade_size,gm.max_grade + 1);
    DenseMultivector dense_out;

    dense_out = init_dense_empty(dense.size);

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(g[gm.grade[i]])
            dense_out.value[i] = dense.value[i];

    PyMem_RawFree(g);
    Py_SET_REFCNT((PyObject*)&dense_out,1);
    return dense_out;
}


static PyMultivectorObject *dense_grade_project(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_dense);
    DenseMultivector *dense = (DenseMultivector*)self->data;

    *dense_out = dense_grade_project_(*dense,self->GA->gm,grades,size);

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;

}



static SparseMultivector sparse_reverse_(SparseMultivector sparse, GradeMap gm){
    SparseMultivector sparse_out = init_sparse_empty(sparse.size);

    for(Py_ssize_t i = 0; i < sparse.size; i++){
        int grade = gm.grade[sparse.bitmap[i]];
        int sign =  (1 - (grade*(grade-1)/2 % 2)*2);
        sparse_out.value[i] = sign*sparse.value[i];
        sparse_out.bitmap[i] = sparse.bitmap[i];
    }

    Py_SET_REFCNT((PyObject*)&sparse_out,1);
    return sparse_out;
}

static PyMultivectorObject *sparse_reverse(PyMultivectorObject *self){
    SparseMultivector *sparse_out = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_sparse);
    SparseMultivector *sparse = (SparseMultivector*)self->data;

    *sparse_out = sparse_reverse_(*sparse,self->GA->gm);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static BladesMultivector blades_reverse_(BladesMultivector blades, GradeMap gm){
    BladesMultivector blades_out = init_blades_empty(blades.size);

    for(Py_ssize_t i = 0; i < blades.size; i++){
        int grade = blades.grade[i];
        int sign =  (1 - (grade*(grade-1)/2 % 2)*2);
        Py_ssize_t bsize = blades.data[i].size;
        blades_out.data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades_out.data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        blades_out.data[i].size = bsize;
        blades_out.grade[i] = grade;
        for(Py_ssize_t j = 0; j < bsize; j++){
            blades_out.data[i].bitmap[j] = blades.data[i].bitmap[j];
            blades_out.data[i].value[j] = sign*blades.data[i].value[j];
        }
    }

    Py_SET_REFCNT((PyObject*)&blades_out,1);
    return blades_out;
}

static PyMultivectorObject *blades_reverse(PyMultivectorObject *self){
    BladesMultivector *blades_out = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_blades);
    BladesMultivector *blades = (BladesMultivector*)self->data;

    *blades_out = blades_reverse_(*blades,self->GA->gm);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static DenseMultivector dense_reverse_(DenseMultivector dense, GradeMap gm){
    DenseMultivector dense_out = init_dense_empty(dense.size);

    for(Py_ssize_t i = 0; i < dense.size; i++){
        int grade = gm.grade[i];
        int sign =  (1 - (grade*(grade-1)/2 % 2)*2);
        dense_out.value[i] = sign*dense.value[i];
    }

    Py_SET_REFCNT((PyObject*)&dense_out,1);
    return dense_out;
}

static PyMultivectorObject *dense_reverse(PyMultivectorObject *self){
    DenseMultivector *dense_out = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    PyMultivectorObject *out = new_multivector(self,MultivectorType_dense);
    DenseMultivector *dense = (DenseMultivector*)self->data;

    *dense_out = dense_reverse_(*dense,self->GA->gm);

    out->data = dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
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

static PyMultivectorObject *init_scalar(PyObject *old, PyObject *scalar){

    PyMultivectorObject *scalar_mv;
    ga_float *scalar_f;
    if(PyFloat_Check(scalar) || PyLong_Check(scalar)) {
        scalar_mv = new_multivector((PyMultivectorObject*)old,MultivectorType_scalar);
        scalar_f = (ga_float*)PyMem_RawMalloc(sizeof(ga_float));
    }else {
        return NULL; // raise error
    }

    if(PyFloat_Check(scalar)){
        *scalar_f = (ga_float)PyFloat_AsDouble(scalar);
    }else if(PyLong_Check(scalar)){
        *scalar_f = (ga_float)PyLong_AsDouble(scalar);
    }
    scalar_mv->data = scalar_f;
    return scalar_mv;
}


static PyObject *multivector_product(PyObject *left, PyObject *right, ProductType ptype){
    PyMultivectorObject *left_mv = NULL;
    PyMultivectorObject *right_mv = NULL;
    MultivectorType left_type, right_type;
    // check if objects are of the same type or scalars
    int is_left_scalar = -1;
    if((left_mv = init_scalar(right,left))){ // check if left is a scalar
        right_mv = (PyMultivectorObject*)right;
        is_left_scalar = 1;
    }else if((right_mv = init_scalar(left,right))){ // check if right is a scalar
        left_mv = (PyMultivectorObject*)left;
        is_left_scalar = 0;
    }else if(PyObject_TypeCheck(left,Py_TYPE(right))){
        left_mv = (PyMultivectorObject*)left;
        right_mv = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(left_mv->GA != right_mv->GA){
        PyErr_SetString(PyExc_TypeError,"operands must have been generated by the same GA class");
        return NULL;
    }
    left_type = left_mv->type; right_type = right_mv->type;

    gaproductfunc product = left_mv->math_funcs.product[left_type][right_type];
    PyMultivectorObject *out;
    if(product){
        out = product(left_mv,right_mv,ptype);
    }
    else{
        PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
        return NULL; // raise not implemented error
    }

    // free scalar multivector
    if(is_left_scalar == 1){
        PyObject_Del(left_mv);
    }else if(is_left_scalar == 0){
        PyObject_Del(right_mv);
    }
    return (PyObject*)out;
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

static PyMultivectorObject *init_sparse_scalar(PyObject *old, PyObject *scalar){

    PyMultivectorObject *scalar_mv;
    SparseMultivector *scalar_f;
    if(PyFloat_Check(scalar) || PyLong_Check(scalar)) {
        scalar_mv = new_multivector((PyMultivectorObject*)old,MultivectorType_sparse);
        scalar_f = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
        *scalar_f = init_sparse_empty(1);
        *scalar_f->bitmap = 0;
    }else {
        return NULL; // Not a scalar
    }
    if(PyFloat_Check(scalar)){
        *scalar_f->value = (ga_float)PyFloat_AsDouble(scalar);
    }else if(PyLong_Check(scalar)){
        *scalar_f->value = (ga_float)PyLong_AsDouble(scalar);
    }
    scalar_mv->data = scalar_f;
    return scalar_mv;
}

static PyMultivectorObject *multivector_add_subtract(PyObject *left, PyObject *right, int sign){
    PyMultivectorObject *left_mv = NULL;
    PyMultivectorObject *right_mv = NULL;

    // check if objects are of the same type or scalars
    int is_left_scalar = -1;
    if((left_mv = init_sparse_scalar(right,left))){ // check if left is a scalar
        right_mv = (PyMultivectorObject*)right;
        is_left_scalar = 1;
    }else if((right_mv = init_sparse_scalar(left,right))){ // check if right is a scalar
        left_mv = (PyMultivectorObject*)left;
        is_left_scalar = 0;
    }else if(PyObject_TypeCheck(left,Py_TYPE(right))){
        left_mv = (PyMultivectorObject*)left;
        right_mv = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(left_mv->GA != right_mv->GA){
        PyErr_SetString(PyExc_TypeError,"operands must have been generated by the same GA class");
        return NULL;
    }


    gaaddfunc add = left_mv->math_funcs.add[left_mv->type][right_mv->type];
    PyMultivectorObject *out;
    if(add){
        out = add(left_mv,right_mv,sign);
    } else{
        PyErr_SetString(PyExc_NotImplementedError,"The add operation for these types is not implemented");
        return NULL; // raise not implemented error
    }

    // free scalar multivector
    if(is_left_scalar == 1)
        PyObject_Del(left_mv);
    else if(is_left_scalar == 0)
        PyObject_Del(right_mv);

    return out;
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

PyObject* multivector_atomic_add(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    MultivectorType mtype = -1;
    PyMultivectorObject *out = NULL;
    int mixed_type = 0;
    PyMultivectorObject *mvs = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *arg_i = PyTuple_GetItem(args,i);
        mvs[i] = *((PyMultivectorObject*)arg_i);
        if(!PyObject_IsInstance(arg_i,cls)){
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.multivector");
            return NULL;
        }
        if(mtype != -1){
            if(mtype != mvs[i].type)
                mixed_type = 1;
        }
        else
            mtype = mvs[i].type;
    }
    gaatomicfunc add;
    if(mixed_type)
        add = mvs->math_funcs.mixed_atomic_add;
    else
        add = mvs->math_funcs.atomic_add[mtype];

    if(add){
        out = add(mvs,size);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The atomic add operation for these types is not implemented");
        return NULL; // raise not implemented error
    }

    return (PyObject*)out;
}

static PyObject* multivector_atomic_product(PyObject *cls, PyObject *args, ProductType ptype){
    Py_ssize_t size = PyTuple_Size(args);
    MultivectorType mtype = -1;
    PyMultivectorObject *out = NULL;
    int mixed_type = 0;
    PyMultivectorObject *mvs = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *arg_i = PyTuple_GetItem(args,i);
        mvs[i] = *((PyMultivectorObject*)arg_i);
        if(!PyObject_IsInstance(arg_i,cls)){
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.multivector");
            return NULL;
        }
        if(mtype != -1){
            if(mtype != mvs[i].type)
                mixed_type = 1;
        }
        else
            mtype = mvs[i].type;
    }

    gaatomicprodfunc product = NULL;
    if(mixed_type)
        product = mvs->math_funcs.mixed_atomic_product;
    else
        product = mvs->math_funcs.atomic_product[mtype];

    if(product){
        out = product(mvs,size,ptype);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The atomic product operation for these types is not implemented");
        return NULL; // raise not implemented error
    }

    return (PyObject*)out;
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

    .mixed_atomic_add = mixed_type_atomic_add,
    .mixed_atomic_product = mixed_type_atomic_product,

    .product[MultivectorType_sparse][MultivectorType_scalar] = sparse_scalar_product,
    .product[MultivectorType_sparse][MultivectorType_sparse] = sparse_sparse_product,
    .product[MultivectorType_sparse][MultivectorType_blades] = mixed_mixed_product,
    .product[MultivectorType_sparse][MultivectorType_dense] = mixed_mixed_product,

    .product[MultivectorType_scalar][MultivectorType_sparse] = sparse_scalar_product,
    .product[MultivectorType_scalar][MultivectorType_blades] = blades_scalar_product,
    .product[MultivectorType_scalar][MultivectorType_dense] = dense_scalar_product,

    .product[MultivectorType_blades][MultivectorType_scalar] = blades_scalar_product,
    .product[MultivectorType_blades][MultivectorType_sparse] = mixed_mixed_product,
    .product[MultivectorType_blades][MultivectorType_blades] = blades_blades_product,
    .product[MultivectorType_blades][MultivectorType_dense] = mixed_mixed_product,

    .product[MultivectorType_dense][MultivectorType_scalar] = dense_scalar_product,
    .product[MultivectorType_dense][MultivectorType_sparse] = mixed_mixed_product,
    .product[MultivectorType_dense][MultivectorType_blades] = mixed_mixed_product,
    .product[MultivectorType_dense][MultivectorType_dense] = dense_dense_product,


    .add[MultivectorType_sparse][MultivectorType_sparse] = sparse_sparse_add,
    .add[MultivectorType_sparse][MultivectorType_blades] = mixed_mixed_add,
    .add[MultivectorType_sparse][MultivectorType_dense] = mixed_mixed_add,

    .add[MultivectorType_blades][MultivectorType_sparse] = mixed_mixed_add,
    .add[MultivectorType_blades][MultivectorType_blades] = blades_blades_add,
    .add[MultivectorType_blades][MultivectorType_dense] = mixed_mixed_add,

    .add[MultivectorType_dense][MultivectorType_sparse] = mixed_mixed_add,
    .add[MultivectorType_dense][MultivectorType_blades] = mixed_mixed_add,
    .add[MultivectorType_dense][MultivectorType_dense] = dense_dense_add,


    .atomic_add[MultivectorType_sparse] = sparse_atomic_add,
    .atomic_add[MultivectorType_blades] = blades_atomic_add,
    .atomic_add[MultivectorType_dense] = dense_atomic_add,

    .atomic_product[MultivectorType_sparse] = sparse_atomic_product,
    .atomic_product[MultivectorType_blades] = blades_atomic_product,
    .atomic_product[MultivectorType_dense] = dense_atomic_product,

    .grade_project[MultivectorType_sparse] = sparse_grade_project,
    .grade_project[MultivectorType_blades] = blades_grade_project,
    .grade_project[MultivectorType_dense] = dense_grade_project,

    .reverse[MultivectorType_sparse] = sparse_reverse,
    .reverse[MultivectorType_blades] = blades_reverse,
    .reverse[MultivectorType_dense] = dense_reverse,
};

static PyMultivectorData_Funcs multivector_data_fn = {
    .free = {sparse_free,dense_free,blades_free,NULL},
    .init = {sparse_init,dense_init,blades_init,NULL},
};

PyMultivectorObject *init_multivector(int *bitmap, ga_float *value, Py_ssize_t size, PyGeometricAlgebraObject *ga, PyTypeObject *obj_type, MultivectorType type){

    PyMultivectorObject *self = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    self->math_funcs = multivector_math_fn;
    self->data_funcs = multivector_data_fn;

    if(type <= MultivectorTypeMIN || type >= MultivectorTypeMAX)
        return NULL; // raise error

    gainitfunc init = self->data_funcs.init[type];
    if(init)
        self->data = init(bitmap,value,size,ga->gm,ga->product[ProductType_geometric].size);
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

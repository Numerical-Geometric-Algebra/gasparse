#include <Python.h>
#include "gasparse.h"

// returns true if abs(v) < p
static int comp_abs(float v, float p){
    float r = (v < 0) ? -v : v;
    return r < p;
}


static Py_ssize_t *init_grade_size(GradeMap gm){
    Py_ssize_t *grade_size = (Py_ssize_t *)PyMem_RawMalloc((gm.max_grade+1)*sizeof(Py_ssize_t));
    for(Py_ssize_t i = 0; i <= gm.max_grade; i++)
        grade_size[i] = 0;
    return grade_size;
}

static SparseMultivector init_sparse_empty(Py_ssize_t size){
    SparseMultivector y;
    y.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    y.value = (float*)PyMem_RawMalloc(size*sizeof(float));
    y.size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        y.bitmap[i] = -1;
        y.value[i] = 0;
    }
    return y;
}


static DenseMultivector init_dense_empty(Py_ssize_t size){
    DenseMultivector y;
    y.value = (float*)PyMem_RawMalloc(size*sizeof(float));
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

static void sparse_remove_small(SparseMultivector y, float precision, Py_ssize_t *size){
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

static void* blades_init(int *bitmap, float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    BladesMultivector *sparse = (BladesMultivector*)PyMem_RawMalloc(sizeof(BladesMultivector));
    BladesMultivector dense;
    Py_ssize_t *grade_size = init_grade_size(gm);
    Py_ssize_t n_grades = 0;
    Py_ssize_t max_grade = 0;
    // determine how many components for each grade
    for(Py_ssize_t i = 0; i < size; i++){
        Py_ssize_t grade = gm.grade[bitmap[i]];
        grade_size[grade]++;
        if(grade > max_grade)
            max_grade = grade;
    }
    // counts the number of different grades
    for(Py_ssize_t i = 0; i <= max_grade; i++)
        if(grade_size[i] > 0)
            n_grades++;

    dense = init_blades_empty(max_grade + 1);
    *sparse = init_blades_empty(n_grades);

    // convert sparse to grade dense
    for(Py_ssize_t i = 0; i < size; i++){
        Py_ssize_t grade = gm.grade[bitmap[i]];
        if(dense.data[grade].bitmap == NULL)
            dense.data[grade] = init_sparse_empty(grade_size[grade]);

        grade_size[grade]--;
        dense.data[grade].bitmap[grade_size[grade]] = bitmap[i];
        dense.data[grade].value[grade_size[grade]] = value[i];
    }

    n_grades = 0;
    for(Py_ssize_t i = 0; i <= max_grade; i++){
        if(dense.data[i].bitmap != NULL){
            sparse->data[n_grades] = dense.data[i];
            sparse->grade[n_grades] = dense.grade[i];
            Py_SET_REFCNT((PyObject*)(sparse->data + i),1);
            n_grades++;
        }
    }
    PyMem_RawFree(dense.data);
    PyMem_RawFree(dense.grade);
    PyMem_RawFree(grade_size);

    return (void*)sparse;
}


static void *dense_init(int *bitmap, float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    DenseMultivector *dense = (DenseMultivector*)PyMem_RawMalloc(sizeof(DenseMultivector));
    dense->value = (float*)PyMem_RawMalloc(algebra_size*sizeof(float));
    dense->size = algebra_size;
    // set all values to 0
    for(Py_ssize_t i = 0; i < algebra_size; i++)
        dense->value[i] = 0;
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= algebra_size)
            return NULL; // raise error
        dense->value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return (void*)dense;
}

static void *sparse_init(int *bitmap, float *value, Py_ssize_t size, GradeMap gm, Py_ssize_t algebra_size){
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    sparse->value = (float*)PyMem_RawMalloc(size*sizeof(float));
    sparse->bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    sparse->size = size;
    for(Py_ssize_t i = 0; i < size; i++){
        sparse->value[i] = value[i];
        sparse->bitmap[i] = bitmap[i];
    }
    return (void*)sparse;
}


void blades_free(void *y){
    if(y == NULL)
        return;
    BladesMultivector *blades = y;

    if(blades->data){
        for(Py_ssize_t i = 0; i < blades->size; i++){
            if(blades->data[i].bitmap != NULL)
                PyMem_RawFree(blades->data[i].bitmap);
            if(blades->data[i].value != NULL)
                PyMem_RawFree(blades->data[i].value);
        }
        PyMem_RawFree(blades->data);
    }
    PyMem_RawFree(blades->grade);
}

void sparse_free(void *y){
    if(y == NULL)
        return;
    SparseMultivector *blades = y;

    PyMem_RawFree(blades->value);
    PyMem_RawFree(blades->bitmap);
}

void dense_free(void *y){
    if(y == NULL)
        return;
    DenseMultivector *dense = y;
    PyMem_RawFree(dense->value);
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
            out_str = (char*)malloc(len*sizeof(char));

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
  Multivector math operations:
  binary operations:
      a+b -> sparse_sparse_add(2*SparseMultivector): adds two multivectors
      a*b -> sparse_sparse_product(2*SparseMultivector,ProductType): product between two multivectors
      1.234*b -> sparse_scalar_product(SparseMultivector,float,ProductType): product between a multivector and a float
  unary operations:
      ~a ->sparse_reverse(SparseMultivector): reversion of a multivector
      a([1,7]) -> sparse_grade_project(SparseMultivector,grades): Grade projection of a multivectors by grades
  atomic operations:
      sparse_atomic_add(n*SparseMultivector): adds n multivectors
      sparse_atomic_product(n*SparseMultivector,ProductType): product between n multivectors (be aware of precedence of the inner product)

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

static SparseMultivector sparse_scalar_product_(SparseMultivector sparse, float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
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

    float *scalar = NULL;
    SparseMultivector *sparse = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        sparse = (SparseMultivector*)left_mv->data;
        scalar = (float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        sparse = (SparseMultivector*)right_mv->data;
        scalar = (float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *sparse_out = sparse_scalar_product_(*sparse,*scalar,*left_mv->GA,ptype);

    out->data = sparse_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static DenseMultivector dense_scalar_product_(DenseMultivector dense, float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
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

    float *scalar = NULL;
    DenseMultivector *dense = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        dense = (DenseMultivector*)left_mv->data;
        scalar = (float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        dense = (DenseMultivector*)right_mv->data;
        scalar = (float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *dense_out = dense_scalar_product_(*dense,*scalar,*left_mv->GA,ptype);

    out->data = (void*)dense_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}




static BladesMultivector blades_scalar_product_(BladesMultivector blades, float scalar, PyGeometricAlgebraObject ga, ProductType ptype){
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

    float *scalar = NULL;
    BladesMultivector *blades = NULL;

    // check which is the scalar multivector
    if(right_mv->type == MultivectorType_scalar){
        blades = (BladesMultivector*)left_mv->data;
        scalar = (float*)right_mv->data;
    }else if(left_mv->type == MultivectorType_scalar){
        blades = (BladesMultivector*)right_mv->data;
        scalar = (float*)left_mv->data;
    }else{
        return NULL; // raise error
    }

    *blades_out = blades_scalar_product_(*blades,*scalar,*left_mv->GA,ptype);

    out->data = blades_out;
    Py_SET_REFCNT((PyObject*)out,1);
    return out;
}

static SparseMultivector sparse_sparse_product_(SparseMultivector left, SparseMultivector right,PyGeometricAlgebraObject ga, ProductType ptype){
    CliffordMap m = ga.product[ptype];

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    SparseMultivector sparse;
    Py_ssize_t size = 0;

    for(Py_ssize_t i = 0; i < left.size; i++){
        for(Py_ssize_t j = 0; j < right.size; j++){
            int sign = m.sign[left.bitmap[i]][right.bitmap[j]];
            // skip product if sign is null
            if(sign == 0)
                continue;
            Py_ssize_t bitmap = m.bitmap[left.bitmap[i]][right.bitmap[j]];
            float value = left.value[i]*right.value[j];

            // write bitmap once to memory
            if(dense.bitmap[bitmap] == -1){
                dense.bitmap[bitmap] = bitmap;
                size++;// increment size of sparse
            }
            dense.value[bitmap] += value*sign;
        }
    }

    sparse_remove_small(dense,ga.precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);

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

static PyMultivectorObject *mixed_scalar_multiply(PyMultivectorObject *self,float scalar,ProductType ptype){
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

    right_mv_sign = mixed_scalar_multiply(right_mv,(float)sign,ProductType_geometric);

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
    float *scalar_f;
    if(PyFloat_Check(scalar) || PyLong_Check(scalar)) {
        scalar_mv = new_multivector((PyMultivectorObject*)old,MultivectorType_scalar);
        scalar_f = (float*)PyMem_RawMalloc(sizeof(float));
    }else {
        return NULL; // raise error
    }

    if(PyFloat_Check(scalar)){
        *scalar_f = (float)PyFloat_AsDouble(scalar);
    }else if(PyLong_Check(scalar)){
        *scalar_f = (float)PyLong_AsDouble(scalar);
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
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or float");
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
        *scalar_f->value = (float)PyFloat_AsDouble(scalar);
    }else if(PyLong_Check(scalar)){
        *scalar_f->value = (float)PyLong_AsDouble(scalar);
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
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or float");
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
    .product[MultivectorType_sparse][MultivectorType_scalar] = sparse_scalar_product,
    .product[MultivectorType_sparse][MultivectorType_sparse] = sparse_sparse_product,
    .product[MultivectorType_sparse][MultivectorType_blades] = mixed_mixed_product,
    .product[MultivectorType_sparse][MultivectorType_dense] = mixed_mixed_product,

    .product[MultivectorType_scalar][MultivectorType_scalar] = NULL,
    .product[MultivectorType_scalar][MultivectorType_sparse] = sparse_scalar_product,
    .product[MultivectorType_scalar][MultivectorType_blades] = blades_scalar_product,
    .product[MultivectorType_scalar][MultivectorType_dense] = dense_scalar_product,

    .product[MultivectorType_blades][MultivectorType_scalar] = blades_scalar_product,
    .product[MultivectorType_blades][MultivectorType_sparse] = mixed_mixed_product,
    .product[MultivectorType_blades][MultivectorType_blades] = NULL,
    .product[MultivectorType_blades][MultivectorType_dense] = mixed_mixed_product,

    .product[MultivectorType_dense][MultivectorType_scalar] = dense_scalar_product,
    .product[MultivectorType_dense][MultivectorType_sparse] = mixed_mixed_product,
    .product[MultivectorType_dense][MultivectorType_blades] = mixed_mixed_product,
    .product[MultivectorType_dense][MultivectorType_dense] = NULL,


    .add[MultivectorType_sparse][MultivectorType_sparse] = sparse_sparse_add,
    .add[MultivectorType_sparse][MultivectorType_blades] = mixed_mixed_add,
    .add[MultivectorType_sparse][MultivectorType_dense] = mixed_mixed_add,

    .add[MultivectorType_blades][MultivectorType_sparse] = mixed_mixed_add,
    .add[MultivectorType_blades][MultivectorType_blades] = NULL,
    .add[MultivectorType_blades][MultivectorType_dense] = mixed_mixed_add,

    .add[MultivectorType_dense][MultivectorType_sparse] = mixed_mixed_add,
    .add[MultivectorType_dense][MultivectorType_blades] = mixed_mixed_add,
    .add[MultivectorType_dense][MultivectorType_dense] = NULL,


    .atomic_add[MultivectorType_sparse] = sparse_atomic_add,
    .mixed_atomic_add = mixed_type_atomic_add,
    .atomic_product[MultivectorType_sparse] = sparse_atomic_product,
    .mixed_atomic_product = mixed_type_atomic_product,
    .grade_project[MultivectorType_sparse] = sparse_grade_project,
    .reverse[MultivectorType_sparse] = sparse_reverse
};

static PyMultivectorData_Funcs multivector_data_fn = {
    .free = {sparse_free,dense_free,blades_free,NULL},
    .init = {sparse_init,dense_init,blades_init,NULL},
};

PyMultivectorObject *init_multivector(int *bitmap, float *value, Py_ssize_t size, PyGeometricAlgebraObject *ga, PyTypeObject *obj_type, MultivectorType type){

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

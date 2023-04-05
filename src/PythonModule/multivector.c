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

static SparseMultivector init_sparse_empty(size_t size){
    SparseMultivector y;
    y.bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
    y.value = (float*)PyMem_RawMalloc(size*sizeof(float));
    y.size = size;
    for(size_t i = 0; i < size; i++){
        y.bitmap[i] = -1;
        y.value[i] = 0;
    }
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

    for(Py_ssize_t i = 0; i < size; i++)
        value[i] = 0;
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
    DenseMultivector *blades = y;
    PyMem_RawFree(blades->value);
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

static PyObject *sparse_repr(void *data_ptr){
    if(data_ptr == NULL)
        return NULL; // raise error
    SparseMultivector *data = data_ptr;
    char **str_bitmap = (char**)PyMem_RawMalloc(data->size*sizeof(char*));
    char **str_value =  (char**)PyMem_RawMalloc(data->size*sizeof(char*));
    Py_ssize_t len_bitmap = 0, len_value = 0;
    PyObject *out = NULL;

    for(Py_ssize_t i = 0; i < data->size; i++){
        str_bitmap[i] = bitmap_to_string(data->bitmap[i]);
        len_bitmap += strlen(str_bitmap[i]);
        str_value[i] = PyOS_double_to_string((double)data->value[i],'f',6,0,NULL);
        len_value += strlen(str_value[i]);
    }
    len_value += data->size + 3;
    len_bitmap += data->size + 3;
    char *format_value = (char*)PyMem_RawMalloc((len_value)*sizeof(char));
    char *format_bitmap = (char*)PyMem_RawMalloc((len_bitmap)*sizeof(char));

    Py_ssize_t j = 0, k = 0, i = 0;
    for(i = 0; i < data->size - 1; i++){
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

    char *format = ".multivector([%s],blades=[%s])";
    size_t size = len_value + len_bitmap + strlen(format);
    char *format_out = (char*)PyMem_RawMalloc(size*sizeof(char));

    PyOS_snprintf(format_out,size,format,format_value,format_bitmap);
    out = Py_BuildValue("s",format_out);

    for(Py_ssize_t i = 0; i < data->size; i++){
        PyMem_RawFree(str_bitmap[i]);
        PyMem_Free(str_value[i]);
    }
    PyMem_RawFree(str_bitmap);
    PyMem_RawFree(str_value);
    PyMem_RawFree(format_value);
    PyMem_RawFree(format_bitmap);
    PyMem_RawFree(format_out);

    return out;
}

static PyObject *blades_repr(void *data){
    Py_RETURN_NONE;
}

static PyObject *dense_repr(void *data){
    Py_RETURN_NONE;
}

static PyMultivectorObject *new_multivector(PyMultivectorObject *old, MultivectorType type){ // return a multivector of type type
    PyMultivectorObject *self = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    Py_SET_REFCNT((PyObject*)self,1);
    self->funcs = old->funcs;
    self->type = type;
    Py_SET_TYPE(self,Py_TYPE(old));
    Py_XINCREF(Py_TYPE(self));
    self->GA = old->GA;
    Py_XINCREF((PyObject*)self->GA);
    self->data = NULL;
    return self;
}



static PyMultivectorObject *sparse_sparse_product(PyMultivectorObject *left_v, PyMultivectorObject *right_v, ProductType ptype){
    PyGeometricAlgebraObject *ga = left_v->GA;
    CliffordMap m = ga->product[ptype];

    SparseMultivector *left = left_v->data, *right = right_v->data;
    size_t m_size = m.size;

    int left_size = left->size;
    int right_size = right->size;

    // Allocate memory for a dense y
    SparseMultivector dense_y = init_sparse_empty(m_size);
    SparseMultivector *sparse_y = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
    PyMultivectorObject *out = new_multivector(left_v,MultivectorType_sparse);
    Py_ssize_t sparse_size = 0;

    for(Py_ssize_t i = 0; i < left_size; i++){
        for(Py_ssize_t j = 0; j < right_size; j++){
            int sign = m.sign[left->bitmap[i]][right->bitmap[j]];
            // skip product if sign is null
            if(sign == 0)
                continue;
            Py_ssize_t bitmap = m.bitmap[left->bitmap[i]][right->bitmap[j]];
            float value = left->value[i]*right->value[j];

            // write bitmap once to memory
            if(dense_y.bitmap[bitmap] == -1){
                dense_y.bitmap[bitmap] = bitmap;
                sparse_size++;// increment size of sparse
            }
            dense_y.value[bitmap] += value*sign;
        }
    }

    sparse_remove_small(dense_y,ga->precision,&sparse_size);
    *sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    sparse_free((void*)&dense_y);
    out->data = sparse_y;
    Py_SET_REFCNT((PyObject*)out->data,1);

    return out;
}



PyObject *multivector_repr(PyMultivectorObject *self){
    PyObject *mv_repr;
    PyObject *ga_repr;

    gareprfunc repr = self->funcs.repr[self->type];
    if(repr){
        mv_repr = repr(self->data);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The repr for this type is not implemented");
        return NULL; // raise error
    }
    ga_repr = PyObject_Repr((PyObject*)self->GA);
    return PySequence_Concat(ga_repr,mv_repr);
}


PyObject *multivector_geometric_product(PyObject *left, PyObject *right){
    // check if objects are of the same type
    if(!PyObject_TypeCheck(left,Py_TYPE(right))){
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type");
        return NULL;
    }
    PyMultivectorObject *left_mv = (PyMultivectorObject*)left;
    PyMultivectorObject *right_mv = (PyMultivectorObject*)right;

    if(left_mv->GA != right_mv->GA){
        PyErr_SetString(PyExc_TypeError,"operands must have been generated by the same GA class");
        return NULL;
    }
    MultivectorType left_type = left_mv->type, right_type = right_mv->type;

    gaproductfunc product = left_mv->funcs.product[left_type][right_type];
    PyMultivectorObject *out;
    if(product){
        out = product(left_mv,right_mv,ProductType_geometric);
    }
    else{
        PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
        return NULL; // raise not implemented error
    }
    return (PyObject*)out;
}

void multivector_dealloc(PyMultivectorObject *self){
    Py_XDECREF((PyObject*)self->GA);
    gafreefunc free_type = self->funcs.free[self->type];
    if(free_type)
       free_type(self->data);

    if(self->data) PyMem_RawFree(self->data);
    PyMem_RawFree(self);
}


static PyMultivector_Funcs multivector_type_fn = {
    .free = {sparse_free,dense_free,blades_free,NULL},
    .init = {sparse_init,dense_init,blades_init,NULL},
    .repr = {sparse_repr,dense_repr,blades_repr,NULL},
    .product[MultivectorType_sparse][MultivectorType_sparse] = sparse_sparse_product,// product between two sparse multivectors
};

PyMultivectorObject *init_multivector(int *bitmap, float *value, Py_ssize_t size, PyGeometricAlgebraObject *ga, PyTypeObject *obj_type, MultivectorType type){

    PyMultivectorObject *self = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    Py_SET_REFCNT((PyObject*)self,1);
    self->funcs = multivector_type_fn;

    if(type <= MultivectorTypeMIN || type >= MultivectorTypeMAX)
        return NULL; // raise error

    gainitfunc init = self->funcs.init[type];
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

/* #include "gasparse.h" */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "gasparse.h"

// returns true if abs(v) == p
static int comp_abs_eq(int v,int p){
    int r = (v < 0) ? -v : v;
    return r == p;
}


// Determines the number of one bits in an integer
// Which is equivalent to determining the grade of a bitset
Py_ssize_t grade(Py_ssize_t v){
    Py_ssize_t c = 0; // c accumulates the total bits set in v
    for (c = 0; v; c++)
        v &= v - 1; // clear the least significant bit set
    return c;
}


static void clifford_sub_algebra(Py_ssize_t k, char **s, int metric){
    Py_ssize_t m = 1 << k; // same as 2^k
    Py_ssize_t n = m << 1; // same as 2^(k+1)
    int sign;
    // This could be improved by checking if the element in the array is zero
    for(Py_ssize_t i = m; i < n; i++){// loop through the new elements
        for(Py_ssize_t j = 0; j < m; j++){// loop through old elements
            // j is indepedent of the new basis vector
            sign = ((grade(j) & 1) == 0) ? 1 : -1;// return minus one if grade is odd
            s[i][j] = sign*s[i-m][j];
            s[j][i] = s[j][i-m];
        }
        if(metric != 0 ){
            for(Py_ssize_t j = m; j < n; j++){// loop through new elements
                // These elements have the new basis vector in common
                sign = metric;
                // remove the new basis vector then determine sign
                sign *= ((grade(j-m) & 1) == 0) ? 1 : -1;
                sign *= s[i-m][j-m];// remove the new vector part
                s[i][j] = sign;
            }
        }else{//if null metric -> set all elements to zero
            for(Py_ssize_t j = m; j < n; j++)
                s[i][j] = 0;
        }
    }
}

static void map_dealloc(CliffordMap *self){
    if(self->sign){
        for(Py_ssize_t i = 0; i < self->size; i++)
            if(self->sign[i])
                PyMem_RawFree(self->sign[i]), self->sign[i] = NULL;
        PyMem_RawFree(self->sign);
        self->sign = NULL;
    }

    if(self->bitmap){
        for(Py_ssize_t i = 0; i < self->size; i++)
            if(self->bitmap[i])
                PyMem_RawFree(self->bitmap[i]), self->bitmap[i] = NULL;
        PyMem_RawFree(self->bitmap);
        self->bitmap = NULL;
    }
}

static void map_alloc(CliffordMap *m, Py_ssize_t nitems){
    char **sign;
    Py_ssize_t **bitmap;

    if(nitems <= 0){
        m->sign = NULL;
        m->bitmap = NULL;
        m->size = 0;
        return;
    }
    sign = (char**)PyMem_RawMalloc(nitems*sizeof(char*));
    bitmap = (Py_ssize_t**)PyMem_RawMalloc(nitems*sizeof(Py_ssize_t*));

    for(Py_ssize_t i = 0; i < nitems; i++){
        sign[i] = (char*)PyMem_RawMalloc(nitems*sizeof(char));
        bitmap[i] = (Py_ssize_t*)PyMem_RawMalloc(nitems*sizeof(Py_ssize_t));
    }

    m->bitmap = bitmap;
    m->sign = sign;
    m->size = nitems;
}

static void map_new(CliffordMap *map){
    map->sign = NULL;
    map->bitmap = NULL;
    map->size = -1;
}

static void map_init(CliffordMap *map, char *metric, Py_ssize_t size){
    Py_ssize_t nitems = 1 << size;
    map_alloc(map,nitems); // alloc memory for the map
    map->sign[0][0] = 1;// initialize algebra of scalars
    for(Py_ssize_t i = 0; i < size; i++)
        clifford_sub_algebra(i,map->sign,metric[i]);

    // determine each basis blade and its grade
    for(Py_ssize_t i = 0; i < nitems; i++){
        for(Py_ssize_t j = 0; j < nitems; j++){
            Py_ssize_t bitmap_ij = i ^ j;
            map->bitmap[i][j] = bitmap_ij;
        }
    }
}

static void copy_map(CliffordMap *dest, CliffordMap *origin){
    for(Py_ssize_t i = 0; i < origin->size; i++){
        for(Py_ssize_t j = 0; j < origin->size; j++){
            dest->bitmap[i][j] = origin->bitmap[i][j];
            dest->sign[i][j] = origin->sign[i][j];
        }
    }
}


static void map_add_basis(CliffordMap *map, char *metric, Py_ssize_t size, Py_ssize_t beg){
    Py_ssize_t nitems = 1 << size;
    CliffordMap temp;
    map_alloc(&temp,nitems); // alloc memory for the new map
    copy_map(&temp,map); // copy the contents of the old map
    map_dealloc(map); // PyMem_RawFree the old map memory

    for(Py_ssize_t i = beg; i < size; i++)
        clifford_sub_algebra(i,temp.sign,metric[i]);

    // determine each basis blade and its grade
    for(Py_ssize_t i = 0; i < nitems; i++){
        for(Py_ssize_t j = 0; j < nitems; j++){
            Py_ssize_t bitmap_ij = i ^ j;
            temp.bitmap[i][j] = bitmap_ij;
        }
    }

    map->bitmap = temp.bitmap;
    map->sign = temp.sign;
    map->size = nitems;
}


static void grade_map_init(GradeMap *m, Py_ssize_t size){

    Py_ssize_t max_grade = grade(size-1);
    Py_ssize_t *g_pos = (Py_ssize_t*)PyMem_RawMalloc((max_grade + 1)*sizeof(Py_ssize_t));
    m->grade = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    m->position = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    for(Py_ssize_t i = 0; i <= max_grade; i++)
        g_pos[i] = 0;

    for(Py_ssize_t i = 0; i < size; i++){
        m->grade[i] = grade(i);
        m->position[i] = g_pos[m->grade[i]]++; // assign value then increment
    }
    m->size = size;
    m->grade_size = g_pos;
    m->max_grade = max_grade;
}


static void inner_map_init(PyGeometricAlgebraObject *self){
    Py_ssize_t size = self->product[ProductType_geometric].size;
    GradeMap gm = self->gm;
    CliffordMap m = self->product[ProductType_geometric];
    map_alloc(&self->product[ProductType_inner],size);
    for(Py_ssize_t i = 0; i < size; i++){
        for(Py_ssize_t j = 0; j < size; j++){
            if(comp_abs_eq((int)gm.grade[i] - (int)gm.grade[j],(int)gm.grade[m.bitmap[i][j]])){
                self->product[ProductType_inner].bitmap[i][j] = m.bitmap[i][j];
                self->product[ProductType_inner].sign[i][j] = m.sign[i][j];
            }else{
                self->product[ProductType_inner].sign[i][j] = 0;
                self->product[ProductType_inner].bitmap[i][j] = -1;
            }
        }
    }
}

static void outer_map_init(PyGeometricAlgebraObject *self){
    Py_ssize_t size = self->product[ProductType_geometric].size;
    GradeMap gm = self->gm;
    CliffordMap m = self->product[ProductType_geometric];
    map_alloc(&self->product[ProductType_outer],size);
    for(Py_ssize_t i = 0; i < size; i++){
        for(Py_ssize_t j = 0; j < size; j++){
            if(gm.grade[i] + gm.grade[j] == gm.grade[m.bitmap[i][j]]){
                self->product[ProductType_outer].bitmap[i][j] = m.bitmap[i][j];
                self->product[ProductType_outer].sign[i][j] = m.sign[i][j];
            }else{
                self->product[ProductType_outer].sign[i][j] = 0;
                self->product[ProductType_outer].bitmap[i][j] = -1;
            }
        }
    }
}

static void inverted_map_init(CliffordMap *inv, CliffordMap *self){
    Py_ssize_t size = self->size;
    map_alloc(inv,size);
    for(Py_ssize_t i = 0; i < size; i++){
        for(Py_ssize_t j = 0; j < size; j++){
            inv->bitmap[i][self->bitmap[i][j]] = j;
            inv->sign[i][self->bitmap[i][j]] = self->sign[i][j];
        }
    }
}

static void geometric_algebra_dealloc(PyGeometricAlgebraObject *self){
    for(Py_ssize_t i = ProductTypeMIN + 1; i < ProductTypeMAX; i++)
        map_dealloc(&self->product[i]);

    if(self->gm.grade) PyMem_RawFree(self->gm.grade),self->gm.grade = NULL;
    if(self->gm.position) PyMem_RawFree(self->gm.position), self->gm.position = NULL;
    if(self->gm.grade_size) PyMem_RawFree(self->gm.grade_size), self->gm.grade_size = NULL;
    if(self->metric) PyMem_RawFree(self->metric), self->metric = NULL;
    Py_TYPE(self)->tp_free((PyObject *)self);
}


static PyObject *geometric_algebra_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    PyGeometricAlgebraObject *self;
    self = (PyGeometricAlgebraObject*)type->tp_alloc(type,0);
    if(self){
        self->p = 0, self->q = 0, self->r = 0;
        self->gm.grade = NULL;
        self->gm.position = NULL;
        self->gm.grade_size = NULL;
        self->gm.max_grade = -1;
        self->gm.size = -1;
        self->metric = NULL;
        for(Py_ssize_t i = ProductTypeMIN + 1; i < ProductTypeMAX; i++)
            map_new(&self->product[i]);

    }
    return (PyObject*) self;
}


static int geometric_algebra_init(PyGeometricAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"p","q","r","metric","print_type","print_type_mv",NULL};
    int p = 0, q = 0, r = 0; PrintType print_type = PrintTypeMIN; PrintTypeMV print_type_mv = PrintTypeMVMIN;
    PyObject *metric_obj = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "|iiiOii", kwlist, &p, &q, &r, &metric_obj, &print_type, &print_type_mv))
        return -1;
    if(!metric_obj){
        if(p < 0 || q < 0 || r < 0)
            return -1;
        else if(p == 0 && q == 0 && r == 0)
            p = 3;

        self->metric = (char*)PyMem_RawMalloc((p+q+r)*sizeof(char));
        for(Py_ssize_t i = 0; i < p; i++) self->metric[i] = 1;
        for(Py_ssize_t i = p; i < p + q; i++) self->metric[i] = -1;
        for(Py_ssize_t i = p + q; i < p + q + r; i++) self->metric[i] = 0;
    }else{
        p = 0, q = 0, r = 0;
        if(!PyList_Check(metric_obj)){
            PyErr_SetString(PyExc_TypeError,"metric must be a list");
            return -1;
        }
        Py_ssize_t size = PyList_Size(metric_obj);
        self->metric = (char*)PyMem_RawMalloc(size*sizeof(char));

        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *value = PyList_GetItem(metric_obj,i);
            if(!PyLong_Check(value)){
                PyErr_SetString(PyExc_TypeError,"items of the list metric must be of type int");
                PyMem_RawFree(self->metric);
                return -1;
            }else {
                int overflow;
                self->metric[i] = (char)PyLong_AsLongAndOverflow(value,&overflow);
                if(overflow){
                    PyMem_RawFree(self->metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return -1;
                }
                if(self->metric[i] == 1){
                    p++;
                }else if(self->metric[i] == -1){
                    q++;
                }else if(self->metric[i] == 0){
                    r++;
                }else {
                    PyMem_RawFree(self->metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return -1;
                }
            }
        }
    }
    map_init(&self->product[ProductType_geometric],self->metric,p+q+r);
    grade_map_init(&self->gm,self->product[ProductType_geometric].size);
    inner_map_init(self);
    outer_map_init(self);

    inverted_map_init(&self->product[ProductType_inverted],&self->product[ProductType_geometric]);
    // can also compute the inverse map of the other products...
    self->p = p; self->q = q; self->r = r; self->precision = 1e-6;

    // set ga print type
    self->print_type = PrintType_metric_array;
    if(print_type < PrintTypeMAX && print_type > PrintTypeMIN)
        self->print_type = print_type;

    // set multivector print type
    self->print_type_mv = PrintTypeMV_reduced;
    if(print_type_mv < PrintTypeMVMAX && print_type_mv > PrintTypeMVMIN)
        self->print_type_mv = print_type_mv;

    return 0;
}



static PyObject *geometric_algebra_add_basis(PyGeometricAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"p","q","r","metric",NULL};
    int p = 0, q = 0, r = 0;
    PyObject *metric_obj = NULL;
    char *metric;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "|iiiO", kwlist, &p, &q, &r, &metric_obj))
        return NULL;
    if(!metric_obj){
        if(p < 0 || q < 0 || r < 0){
            PyErr_SetString(PyExc_ValueError,"p,q and r can't be negative");
            return NULL;
        } else if(p == 0 && q == 0 && r == 0){
            PyErr_SetString(PyExc_ValueError,"There should be at least one non zero value for p,q or r");
            return NULL;
        }

        metric = (char*)PyMem_RawMalloc((Metric_SIZE(self) + p + q + r)*sizeof(char));
        for(Py_ssize_t i = 0; i < Metric_SIZE(self); i++)
            metric[i] = self->metric[i];
        for(Py_ssize_t i = Metric_SIZE(self); i < p + Metric_SIZE(self); i++)
            metric[i] = 1;
        for(Py_ssize_t i = p + Metric_SIZE(self); i < p + q + Metric_SIZE(self); i++)
            metric[i] = -1;
        for(Py_ssize_t i = p + q + Metric_SIZE(self); i < p + q + r + Metric_SIZE(self); i++)
            metric[i] = 0;
        PyMem_RawFree(self->metric);
        self->metric = metric;
    }else{
        if(!PyList_Check(metric_obj)){
            PyErr_SetString(PyExc_TypeError,"metric must be a list");
            return NULL;
        }
        Py_ssize_t size = PyList_Size(metric_obj);
        metric = (char*)PyMem_RawMalloc((Metric_SIZE(self) + size)*sizeof(char));

        for(Py_ssize_t i = 0; i < Metric_SIZE(self); i++)
            metric[i] = self->metric[i];

        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *value = PyList_GetItem(metric_obj,i);
            if(!PyLong_Check(value)){
                PyErr_SetString(PyExc_TypeError,"items of the list metric must be of type int");
                PyMem_RawFree(metric);
                return NULL;
            }else {
                int overflow;
                metric[i + Metric_SIZE(self)] = (char)PyLong_AsLongAndOverflow(value,&overflow);
                if(overflow){
                    PyMem_RawFree(metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return NULL;
                }
                if(metric[i] == 1){
                    p++;
                }else if(metric[i] == -1){
                    q++;
                }else if(metric[i] == 0){
                    r++;
                }else {
                    PyMem_RawFree(metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return NULL;
                }
            }
        }
    }

    map_add_basis(&self->product[ProductType_geometric],metric,Metric_SIZE(self) + p + q + r, Metric_SIZE(self));
    // PyMem_RawFree old product maps
    for(Py_ssize_t i = ProductType_geometric + 1; i < ProductTypeMAX; i++)
        map_dealloc(&self->product[i]);

    if(self->gm.grade) PyMem_RawFree(self->gm.grade);
    if(self->gm.position) PyMem_RawFree(self->gm.position);
    if(self->gm.grade_size) PyMem_RawFree(self->gm.grade_size);
    // compute all the auxiliar maps again
    grade_map_init(&self->gm,self->product[ProductType_geometric].size);
    inner_map_init(self);
    outer_map_init(self);
    inverted_map_init(&self->product[ProductType_inverted],&self->product[ProductType_geometric]);
    self->p += p; self->q += q; self->r += r;

    Py_RETURN_NONE;
}

static PyObject *geometric_algebra_repr(PyGeometricAlgebraObject *self){
    char  str[100];
    if(self->print_type == PrintType_metric){
        PyOS_snprintf(str,100,"GA(p=%lu,q=%lu,r=%lu)",self->p,self->q,self->r);
        return Py_BuildValue("s",str);
    }else if(self->print_type == PrintType_metric_array){
        PyObject *out;
        Py_ssize_t size = 3*Metric_SIZE(self) + 2;
        char *str_temp = (char*)PyMem_RawMalloc(size*sizeof(char));
        char **str_value = (char**)PyMem_RawMalloc(Metric_SIZE(self)*sizeof(char*));
        Py_ssize_t i;
        for(i = 0; i < Metric_SIZE(self)-1; i++){
            str_value[i] = (char*)PyMem_RawMalloc(4*sizeof(char));
            PyOS_snprintf(str_value[i],4,"%+d,",(int)self->metric[i]);
        }
        str_value[i] = (char*)PyMem_RawMalloc(4*sizeof(char));
        PyOS_snprintf(str_value[i],4,"%+d",(int)self->metric[i]);

        Py_ssize_t j = 0;
        for(i = 0; i < Metric_SIZE(self); i++){
            strcpy(str_temp + j,str_value[i]);
            j += strlen(str_value[i]);
        }
        str_temp[j] = '\0';
        PyOS_snprintf(str,100,"GA(metric=[%s])",str_temp);
        out = Py_BuildValue("s",str);

        for(Py_ssize_t i = 0; i < Metric_SIZE(self); i++){
            PyMem_RawFree(str_value[i]);
        }PyMem_RawFree(str_value);
        PyMem_RawFree(str_temp);

        return out;
    }else if(self->print_type == PrintType_vectors){
        PyObject *out;
        char str_temp[100];
        char **str_value = (char**)PyMem_RawMalloc(Metric_SIZE(self)*sizeof(char*));
        Py_ssize_t i;
        for(i = 0; i < Metric_SIZE(self)-1; i++){
            str_value[i] = (char*)PyMem_RawMalloc(7*sizeof(char));
            PyOS_snprintf(str_value[i],7,"e%lu:%+d,",i,(int)self->metric[i]);
        }
        str_value[i] = (char*)PyMem_RawMalloc(7*sizeof(char));
        PyOS_snprintf(str_value[i],7,"e%lu:%+d",i,(int)self->metric[i]);

        Py_ssize_t j = 0;
        for(i = 0; i < Metric_SIZE(self); i++){
            strcpy(str_temp + j, str_value[i]);
            j += strlen(str_value[i]);
        }
        str_temp[j] = '\0';
        PyOS_snprintf(str,100,"GA(vectors=[%s])",str_temp);
        out = Py_BuildValue("s",str);

        for(Py_ssize_t i = 0; i < Metric_SIZE(self); i++){
            PyMem_RawFree(str_value[i]);
        }PyMem_RawFree(str_value);

        return out;
    }
    Py_RETURN_NONE;
}

static int parse_list_as_bitmaps(PyObject *blades, int **bitmap){

    if(!PyList_Check(blades))
        return -1;

    Py_ssize_t size = PyList_Size(blades);
    *bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));

    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *blade_i = PyList_GetItem(blades,i);
        if(!PyUnicode_Check(blade_i)){
            PyMem_RawFree(*bitmap);
            return -1;
        }

        const char *blade_str = PyUnicode_AsUTF8(blade_i);
        Py_ssize_t len = strlen(blade_str);
        Py_ssize_t j = 0;
        if(*blade_str == 'e')
            j++;
        Py_ssize_t bitmap_i = 0;
        for(; j < len; j++){
            if(!Is_NonZero(blade_str[j]))
                return -1;
            bitmap_i += 1 << (int)(blade_str[j] - '1');
        }
        (*bitmap)[i] = bitmap_i;
    }
    return size;
}

static PyNumberMethods PyMultivectorNumberMethods = {
    .nb_multiply = (binaryfunc) multivector_geometric_product,
    .nb_xor = (binaryfunc) multivector_outer_product,
    .nb_or = (binaryfunc) multivector_inner_product,
    .nb_add = (binaryfunc) multivector_add,
    .nb_subtract = (binaryfunc) multivector_subtract,
    .nb_invert = (unaryfunc) multivector_invert
};

PyDoc_STRVAR(add_doc, "adds a bunch of multivectors.");
PyDoc_STRVAR(product_doc, "multiplies a bunch of multivectors.");

PyMethodDef multivector_methods[] = {
    {"add", (PyCFunction) multivector_atomic_add, METH_VARARGS|METH_CLASS, add_doc},
    {"geometric_product", (PyCFunction) multivector_atomic_geometric_product, METH_VARARGS|METH_CLASS, product_doc},
    {"outer_product", (PyCFunction) multivector_atomic_outer_product, METH_VARARGS|METH_CLASS, product_doc},
    {NULL},
};


static PyTypeObject PyMultivectorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gasparse.multivector",
    .tp_doc = PyDoc_STR("Builds a multivector in different types (sparse,dense,blades)"),
    .tp_basicsize = sizeof(PyMultivectorObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_dealloc = (destructor) multivector_dealloc,
    .tp_repr = (reprfunc) multivector_repr,
    .tp_str = (reprfunc) multivector_repr,
    .tp_call = (ternaryfunc) multivector_grade_project,
    .tp_new = NULL,
    .tp_as_number = &PyMultivectorNumberMethods,
    .tp_methods = multivector_methods
};


static PyObject *geometric_algebra_multivector(PyGeometricAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"values","blades","precision","dtype","addmethod",NULL};
    PyObject *values = NULL, *blades = NULL, *dtype = NULL, *addmethod = NULL, *out = NULL;
    int precision;
    int *bitmaps_int = NULL;
    float *values_float;
    Py_ssize_t size;
    MultivectorType type = MultivectorType_sparse; // default type

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iOO", kwlist, &values,&blades,&precision,&dtype,&addmethod))
        return NULL;
    if(!values)
        return NULL;
    if(!PyList_Check(values)){
        PyErr_SetString(PyExc_TypeError,"values must be a list");
        return NULL;
    }
    size = PyList_Size(values);
    values_float = (float*)PyMem_RawMalloc(size*sizeof(float));
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *value_i = PyList_GetItem(values,i);
        if(PyFloat_Check(value_i))
            values_float[i] = (float)PyFloat_AsDouble(value_i);
        else if(PyLong_Check(value_i))
            values_float[i] = (float)PyLong_AsLong(value_i);
        else{
            PyErr_SetString(PyExc_TypeError,"Elements of the list of values must be float");
            PyMem_RawFree(values_float);
            return NULL;
        }
    }

    if(blades){
        int size_bitmap = parse_list_as_bitmaps(blades,&bitmaps_int);
        if(size_bitmap != (int)size){
            PyErr_SetString(PyExc_TypeError,"blades must be of the same size as values");
            return NULL;
        }
    }

    if(dtype && PyUnicode_Check(dtype)){
        const char *dtype_str = PyUnicode_AsUTF8(dtype);
        if(!strcmp(dtype_str,"sparse")){
            type = MultivectorType_sparse;
        }else if(!strcmp(dtype_str,"dense")){
            type = MultivectorType_dense;
        }else if(!strcmp(dtype_str,"blades")){
            type = MultivectorType_blades;
        }else if(!strcmp(dtype_str,"scalar")){
            type = MultivectorType_scalar;
        }else {
            // raise warning saying the default type
        }
    }
    Py_INCREF(self);
    Py_INCREF((PyObject*)&PyMultivectorType);
    out = (PyObject*)init_multivector(bitmaps_int,values_float,size,self,&PyMultivectorType,type);

    PyMem_RawFree(values_float);
    PyMem_RawFree(bitmaps_int);
    return out;
}


static PyMethodDef geometric_algebra_methods[] = {
    {"add_basis", (PyCFunction) geometric_algebra_add_basis, METH_VARARGS | METH_KEYWORDS,
     "adds basis vectors to the algebra"
    },
    {"multivector",(PyCFunction) geometric_algebra_multivector, METH_VARARGS | METH_KEYWORDS,
     "generate a multivector"
    },
    {NULL}  /* Sentinel */
};


static PyTypeObject PyGeometricAlgebraType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gasparse.GA",
    .tp_doc = PyDoc_STR("Construction of a geometric algebra given a metric"),
    .tp_basicsize = sizeof(PyGeometricAlgebraObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = (newfunc) geometric_algebra_new,
    .tp_init = (initproc) geometric_algebra_init,
    .tp_dealloc = (destructor) geometric_algebra_dealloc,
    .tp_repr = (reprfunc) geometric_algebra_repr,
    .tp_str = (reprfunc) geometric_algebra_repr,
    .tp_methods = geometric_algebra_methods,
};

static PyModuleDef gasparse_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "gasparse",
    .m_doc = "Implementation of sparse multivectors.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_gasparse(void){
    PyObject *m;
    if (PyType_Ready(&PyGeometricAlgebraType) < 0)
        return NULL;

    if (PyType_Ready(&PyMultivectorType) < 0)
        return NULL;

    m = PyModule_Create(&gasparse_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&PyGeometricAlgebraType);
    if (PyModule_AddObject(m, "GA", (PyObject *) &PyGeometricAlgebraType) < 0) {
        Py_DECREF(&PyGeometricAlgebraType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyMultivectorType);
    if (PyModule_AddObject(m, "multivector", (PyObject *)&PyMultivectorType) < 0) {
        Py_DECREF(&PyMultivectorType);
        Py_DECREF(m);
        return NULL;
    }


    return m;
}

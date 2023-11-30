#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "gasparse.h"
/* #include "multivector.h" */
#ifdef INCLUDE_GENCODE
#include "multivector_gen.h"
#endif


// replace this with a builtin
// returns true if abs(v) == p
static int comp_abs_eq(int v,int p){
    int r = (v < 0) ? -v : v;
    return r == p;
}

static void clifford_sub_algebra(Py_ssize_t k, char **s, int metric){
    Py_ssize_t m = 1 << k; // same as 2^k
    Py_ssize_t n = m << 1; // same as 2^(k+1)
    int sign;
    // This could be improved by checking if the element in the array is zero
    for(Py_ssize_t i = m; i < n; i++){// loop through the new elements
        for(Py_ssize_t j = 0; j < m; j++){// loop through old elements
            // j is indepedent of the new basis vector
            sign = ((GRADE(j) & 1) == 0) ? 1 : -1;// return minus one if grade is odd
            s[i][j] = sign*s[i-m][j];
            s[j][i] = s[j][i-m];
        }
        if(metric != 0 ){
            for(Py_ssize_t j = m; j < n; j++){// loop through new elements
                // These elements have the new basis vector in common
                sign = metric;
                // remove the new basis vector then determine sign
                sign *= ((GRADE(j-m) & 1) == 0) ? 1 : -1;
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
            PyMem_RawFree(self->sign[i]), self->sign[i] = NULL;
        PyMem_RawFree(self->sign);
        self->sign = NULL;
    }

    if(self->bitmap){
        for(Py_ssize_t i = 0; i < self->size; i++)
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
    if(!bitmap || !sign){
        m->size = -1;
        PyMem_RawFree(sign);
        PyMem_RawFree(bitmap);
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for the map");
        return;
    }
    m->bitmap = bitmap;
    m->sign = sign;
    for(Py_ssize_t i = 0; i < nitems; i++){
        sign[i] = (char*)PyMem_RawMalloc(nitems*sizeof(char));
        bitmap[i] = (Py_ssize_t*)PyMem_RawMalloc(nitems*sizeof(Py_ssize_t));
        if(!sign[i] || !bitmap[i]){
            m->size = i + 1;
            map_dealloc(m);
            m->size = -1;
            PyErr_SetString(PyExc_MemoryError,"Error allocating memory for the sign or bitmap of the map");
            return;
        }
    }

    m->size = nitems;
}

static void map_sign_alloc(CliffordMap *m, Py_ssize_t nitems){
    char **sign;
    if(nitems <= 0){
        m->sign = NULL;
        m->bitmap = NULL;
        m->size = 0;
        return;
    }
    sign = (char**)PyMem_RawMalloc(nitems*sizeof(char*));
    if(!sign){
        m->size = -1;
        PyMem_RawFree(sign);
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for the map");
        return;
    }
    m->sign = sign;
    for(Py_ssize_t i = 0; i < nitems; i++){
        sign[i] = (char*)PyMem_RawMalloc(nitems*sizeof(char));

        if(!sign[i]){
            m->size = i + 1;
            map_dealloc(m);
            m->size = -1;
            PyErr_SetString(PyExc_MemoryError,"Error allocating memory for the sign of the map");
            return;
        }
    }

    m->bitmap = NULL;
    m->size = nitems;
}

static void map_new(CliffordMap *map){
    map->sign = NULL;
    map->bitmap = NULL;
    map->size = -1;
}

static void map_init(CliffordMap *map, char *metric, Py_ssize_t size){
    if(size == -1) return;
    Py_ssize_t nitems = 1 << size;
    map_alloc(map,nitems); // alloc memory for the map
    if(map->size == -1) return;
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


static void map_sign_init(CliffordMap *map, char *metric, Py_ssize_t size){
    if(size == -1) return;
    Py_ssize_t nitems = 1 << size;
    map_sign_alloc(map,nitems); // alloc memory only for the signs
    if(map->size == -1) return;
    map->sign[0][0] = 1;// initialize algebra of scalars
    for(Py_ssize_t i = 0; i < size; i++)
        clifford_sub_algebra(i,map->sign,metric[i]);
}

static void copy_map(CliffordMap *dest, CliffordMap *origin){
    for(Py_ssize_t i = 0; i < origin->size; i++){
        for(Py_ssize_t j = 0; j < origin->size; j++){
            dest->bitmap[i][j] = origin->bitmap[i][j];
            dest->sign[i][j] = origin->sign[i][j];
        }
    }
}

static void map_reset(CliffordMap *m){
    for(Py_ssize_t i = 0; i < m->size; i++){
        for(Py_ssize_t j = 0; j < m->size; j++){
            m->bitmap[i][j] = -1;
            m->sign[i][j] = 0;
        }
    }
}


static void map_add_basis(CliffordMap *map, char *metric, Py_ssize_t size, Py_ssize_t beg){
    if(size == -1) return;
    Py_ssize_t nitems = 1 << size;
    CliffordMap temp;
    map_alloc(&temp,nitems); // alloc memory for the new map
    if(temp.size == -1) {
        map->size = -1;
        return;
    }
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
    if(size == -1) return;
    Py_ssize_t max_grade = GRADE(size-1);
    Py_ssize_t *g_pos = (Py_ssize_t*)PyMem_RawMalloc((max_grade + 1)*sizeof(Py_ssize_t));
    m->grade = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    m->position = (Py_ssize_t*)PyMem_RawMalloc(size*sizeof(Py_ssize_t));
    if(!m->grade || !m->position || !g_pos){
        m->size = -1;
        return;
    }
    for(Py_ssize_t i = 0; i <= max_grade; i++)
        g_pos[i] = 0;

    for(Py_ssize_t i = 0; i < size; i++){
        m->grade[i] = GRADE(i);
        m->position[i] = g_pos[m->grade[i]]++; // assign value then increment
    }
    m->size = size;
    m->grade_size = g_pos;
    m->max_grade = max_grade;
}


static void inner_map_init(PyAlgebraObject *self){
    Py_ssize_t size = self->product[ProductType_geometric].size;
    if(size == -1) return;
    GradeMap gm = self->gm;
    CliffordMap m = self->product[ProductType_geometric];
    map_alloc(&self->product[ProductType_inner],size);
    if(self->product[ProductType_inner].size == -1) return;
    for(Py_ssize_t i = 0; i < size; i++){
        for(Py_ssize_t j = 0; j < size; j++){
            if(gm.grade[i] == 0 || gm.grade[j] == 0){
                // inner product with grade 0 elements is 0
                self->product[ProductType_inner].sign[i][j] = 0;
                self->product[ProductType_inner].bitmap[i][j] = -1;
            }else if(comp_abs_eq((int)gm.grade[i] - (int)gm.grade[j],(int)gm.grade[m.bitmap[i][j]])){
                self->product[ProductType_inner].bitmap[i][j] = m.bitmap[i][j];
                self->product[ProductType_inner].sign[i][j] = m.sign[i][j];
            }else{
                self->product[ProductType_inner].sign[i][j] = 0;
                self->product[ProductType_inner].bitmap[i][j] = -1;
            }
        }
    }
}

static void outer_map_init(PyAlgebraObject *self){
    Py_ssize_t size = self->product[ProductType_geometric].size;
    GradeMap gm = self->gm;
    CliffordMap m = self->product[ProductType_geometric];
    map_alloc(&self->product[ProductType_outer],size);
    if(self->product[ProductType_outer].size == -1) return;
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

static void regressive_map_init(PyAlgebraObject *self){
    CliffordMap m = self->product[ProductType_outer];
    Py_ssize_t size = m.size;
    DualMap dm = self->dm;
    Py_ssize_t l,r; int lsign,rsign;

    int sign = METRIC_SIZE(self) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    if(dm.size <= 0 || m.size <= 0) return;
    map_alloc(&self->product[ProductType_regressive],size);
    for(Py_ssize_t i = 0; i < size; i++){ //iterate left multivectors
        // compute dual
        l = dm.bitmap[i];
        lsign = dm.sign[i];
        for(Py_ssize_t j = 0; j < size; j++){
            r = dm.bitmap[j];
            rsign = dm.sign[j];
            if(m.bitmap[l][r] != -1){
                self->product[ProductType_regressive].bitmap[i][j] = dm.bitmap[m.bitmap[l][r]];
                self->product[ProductType_regressive].sign[i][j] = sign*rsign*lsign*m.sign[l][r];
            }else{
                self->product[ProductType_regressive].sign[i][j] = 0;
                self->product[ProductType_regressive].bitmap[i][j] = -1;
            }
        }
    }
}


static void inverted_map_init(CliffordMap *inv, CliffordMap *self){
    Py_ssize_t size = self->size;
    map_alloc(inv,size);
    map_reset(inv);
    if(inv->size == -1) return;
    for(Py_ssize_t i = 0; i < size; i++){
        for(Py_ssize_t j = 0; j < size; j++){
            if(self->bitmap[i][j] != -1){
                inv->bitmap[i][self->bitmap[i][j]] = j;
                inv->sign[i][self->bitmap[i][j]] = self->sign[i][j];
            }
        }
    }
}

static DualMap dual_map_init(Py_ssize_t n){
    DualMap dm;
    Py_ssize_t nitems = 1 << n;
    Py_ssize_t pss = nitems - 1;
    int psssignreverse = n & 2 ? -1 : 1;
    dm.bitmap = (Py_ssize_t*)PyMem_RawMalloc(nitems*sizeof(Py_ssize_t));
    dm.sign = (char*)PyMem_RawMalloc(nitems*sizeof(char));
    dm.size = nitems;
    int sum, sign; Py_ssize_t j;
    for(Py_ssize_t i = 0; i < nitems; i++){
        j = i; sum = 0;
        while(j){
            int tz = __builtin_ctzll((unsigned long long)j); // count the number of trailing zeros
            sum += tz;
            j -= 1 << tz; // remove set bit
        }
        sign = sum & 1 ? -1 : 1;
        dm.bitmap[i] = pss ^ i;
        dm.sign[i] = psssignreverse*sign;
    }
    return dm;
}

static DualMap dual_map_sign_init(Py_ssize_t n){
    DualMap dm;
    Py_ssize_t nitems = 1 << n;
    int psssignreverse = n & 2 ? -1 : 1;
    dm.bitmap = NULL;
    dm.sign = (char*)PyMem_RawMalloc(nitems*sizeof(char));
    dm.size = nitems;
    int sum, sign; Py_ssize_t j;
    for(Py_ssize_t i = 0; i < nitems; i++){
        j = i; sum = 0;
        while(j){
            int tz = __builtin_ctzll((unsigned long long)j); // count the number of trailing zeros
            sum += tz;
            j -= 1 << tz; // remove set bit
        }
        sign = sum & 1 ? -1 : 1;
        dm.sign[i] = psssignreverse*sign;
    }
    return dm;
}

static void algebra_dealloc(PyAlgebraObject *self){
    for(Py_ssize_t i = ProductTypeMIN + 1; i < ProductTypeMAX; i++)
        map_dealloc(&self->product[i]);

    PyMem_RawFree(self->gm.grade); self->gm.grade = NULL;
    PyMem_RawFree(self->gm.position); self->gm.position = NULL;
    PyMem_RawFree(self->gm.grade_size); self->gm.grade_size = NULL;
    PyMem_RawFree(self->metric); self->metric = NULL;
    PyMem_RawFree(self->dm.bitmap); self->dm.bitmap = NULL;
    PyMem_RawFree(self->dm.sign); self->dm.sign = NULL;
    PyMem_RawFree(self->mdefault.type_name);
    PyMem_RawFree(self->types);
    Py_TYPE(self)->tp_free((PyObject *)self);
}


static PyObject *algebra_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    PyAlgebraObject *self;
    self = (PyAlgebraObject*)type->tp_alloc(type,0);
    if(self){
        self->p = 0, self->q = 0, self->r = 0;
        self->gm.grade = NULL;
        self->gm.position = NULL;
        self->gm.grade_size = NULL;
        self->gm.max_grade = -1;
        self->gm.size = -1;
        self->dm.size = -1;
        self->dm.sign = NULL;
        self->dm.bitmap = NULL;
        self->metric = NULL;
        self->mdefault.type_name = NULL;
        for(Py_ssize_t i = ProductTypeMIN + 1; i < ProductTypeMAX; i++)
            map_new(&self->product[i]);

    }
    return (PyObject*)self;
}

static ComputationMode get_computation_mode_value(char *name){
    if(!name) return ComputationMode_generic; // default mode
    if(!strcmp(name,"generated")) return ComputationMode_generated;
    if(!strcmp(name,"large")) return ComputationMode_large;
    if(!strcmp(name,"generic")) return ComputationMode_generic;
    if(!strcmp(name,"devgeneration")) return ComputationMode_devgeneration;
    return -1;
}


static int compute_metric(PyAlgebraObject *self, PyObject *metric, int p, int q, int r){
    if(!metric){
        if(p <= 0 && q <= 0 && r <= 0)
            return 0;
        if(p < 0 || q < 0 || r < 0)
            return 0;

        self->metric = (char*)PyMem_RawMalloc((p+q+r)*sizeof(char));
        for(Py_ssize_t i = 0; i < p; i++) self->metric[i] = 1;
        for(Py_ssize_t i = p; i < p + q; i++) self->metric[i] = -1;
        for(Py_ssize_t i = p + q; i < p + q + r; i++) self->metric[i] = 0;
        self->p = p; self->q = q; self->r  = r;
    }else{
        self->p = 0, self->q = 0, self->r = 0;
        if(!PyList_Check(metric)){
            PyErr_SetString(PyExc_TypeError,"metric must be a list");
            return 0;
        }
        Py_ssize_t size = PyList_Size(metric);
        self->metric = (char*)PyMem_RawMalloc(size*sizeof(char));

        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *value = PyList_GetItem(metric,i);
            if(!PyLong_Check(value)){
                PyErr_SetString(PyExc_TypeError,"items of the list metric must be of type int");
                PyMem_RawFree(self->metric);
                return 0;
            }else {
                int overflow;
                self->metric[i] = (char)PyLong_AsLongAndOverflow(value,&overflow);
                if(overflow){
                    PyMem_RawFree(self->metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return 0;
                }
                if(self->metric[i] == 1){
                    self->p++;
                }else if(self->metric[i] == -1){
                    self->q++;
                }else if(self->metric[i] == 0){
                    self->r++;
                }else {
                    PyMem_RawFree(self->metric);
                    PyErr_SetString(PyExc_ValueError,"metric can only have values -1,1,0.");
                    return 0;
                }
            }
        }
    }
    return 1;
}

#ifdef INCLUDE_GENCODE
// looks for a name or a metric in the generated types array
static int get_generated_types_indices(PyAlgebraObject *self, char *name, Py_ssize_t *index){
    Py_ssize_t k = 0;
    if(name && !self->metric){
        for(Py_ssize_t i = 0; i < N_GEN_SUBTYPES; i++)
            if(!strcmp(gen_subtypes_array[i].name,name))
                index[k++] = i;
        if(k){// copy the metric from the type
            Py_ssize_t l = *index;
            Py_ssize_t size = gen_subtypes_array[l].msize;
            self->p = 0; self->q = 0; self->r = 0;
            self->metric = (char*)PyMem_RawMalloc(size*sizeof(char));
            for(Py_ssize_t i = 0; i < size; i++){
                self->metric[i] = gen_subtypes_array[l].metric[i];
                if(self->metric[i] == 1) self->p++;
                if(self->metric[i] == -1) self->q++;
                if(self->metric[i] == 0) self->r++;
            }
            self->asize = gen_subtypes_array[l].asize;
        }else return -1;
    }else if(self->metric){
        for(Py_ssize_t i = 0; i < N_GEN_SUBTYPES; i++){
            if(METRIC_SIZE(self) == gen_subtypes_array[i].msize){
                char *metric = gen_subtypes_array[i].metric;
                int check = 1;
                Py_ssize_t size = gen_subtypes_array[i].msize;
                for(Py_ssize_t j = 0; j < size; j++){
                    // compare metric values
                    if(metric[j] != self->metric[j]){
                        check = 0;
                        break;
                    }
                }
                if(check) index[k++] = i;
            }
        }
        if(k) self->asize = gen_subtypes_array[*index].asize;
    }
    return k;
}
#endif

static int algebra_init(PyAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"p","q","r","metric","print_type","print_type_mv","compute_mode","name",NULL};
    int p = 0, q = 0, r = 0; PrintType print_type = PrintTypeMIN; PrintTypeMV print_type_mv = PrintTypeMVMIN;
    PyObject *metric = NULL;
    char *mode_name = NULL, *algebra_name = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "|iiiOiiss", kwlist,
                                    &p, &q, &r,
                                    &metric,
                                    &print_type,
                                    &print_type_mv,
                                    &mode_name,
                                    &algebra_name)) return -1;

    compute_metric(self,metric,p,q,r);

    ComputationMode mode = get_computation_mode_value(mode_name);
    switch(mode){
        case ComputationMode_generated:;
#ifndef INCLUDE_GENCODE
            PyErr_SetString(PyExc_ValueError,
                            "Generated computation mode not available in this module use gasparsegen!");
            return -1;
#else
            Py_ssize_t k = -1;
            Py_ssize_t index[10] = {-1};
            if((k = get_generated_types_indices(self,algebra_name,index)) <= 0){
                PyErr_SetString(PyExc_ValueError,
                                "The asked algebra is not available generated try the other computation modes");
                return -1; // algebra not found error
            }
            self->types = (PyMultivectorSubType*)PyMem_RawMalloc(k*sizeof(PyMultivectorSubType));
            self->tsize = k;
            for(Py_ssize_t i = 0; i < k; i++)
                self->types[i] = gen_subtypes_array[index[i]];
            self->mixed = &cast_multivector_mixed_fn;
            // set default type
            self->mdefault.type_name = (char*)PyMem_RawMalloc((strlen("dense")+1)*sizeof(char));
            strcpy(self->mdefault.type_name,"dense");
#endif
            break;
        case ComputationMode_large:;
            if(!self->metric) return -1;
            fill_missing_funcs();
            map_sign_init(self->product,self->metric,METRIC_SIZE(self));
            self->dm = dual_map_sign_init(METRIC_SIZE(self));
            self->asize = self->product->size;
            self->tsize = 3;
            self->types = (PyMultivectorSubType*)PyMem_RawMalloc(self->tsize*sizeof(PyMultivectorSubType));
            for(Py_ssize_t i = 0; i < self->tsize; i++)
                self->types[i] = largemultivector_subtypes_array[i];
            self->mixed = &largemultivector_mixed_fn;
            // set default type
            self->mdefault.type_name = (char*)PyMem_RawMalloc((strlen("sparse")+1)*sizeof(char));
            strcpy(self->mdefault.type_name,"sparse");
            break;
        case ComputationMode_generic:;
            if(!self->metric) return -1;
            map_init(&self->product[ProductType_geometric],self->metric,METRIC_SIZE(self));
            self->dm = dual_map_init(METRIC_SIZE(self));
            self->asize = self->product->size;

            grade_map_init(&self->gm,self->asize);
            inner_map_init(self);
            outer_map_init(self);
            regressive_map_init(self);

            self->tsize = 3;
            self->types = (PyMultivectorSubType*)PyMem_RawMalloc(self->tsize*sizeof(PyMultivectorSubType));
            for(Py_ssize_t i = 0; i < self->tsize; i++)
                self->types[i] = multivector_subtypes_array[i];
            self->mixed = &multivector_mixed_fn;
            // set default type
            self->mdefault.type_name = (char*)PyMem_RawMalloc((strlen("sparse")+1)*sizeof(char));
            strcpy(self->mdefault.type_name,"sparse");

            break;
        case ComputationMode_devgeneration:;
            if(!self->metric) return -1;
            map_init(&self->product[ProductType_geometric],self->metric,METRIC_SIZE(self));
            self->dm = dual_map_init(METRIC_SIZE(self));
            self->asize = self->product->size;

            grade_map_init(&self->gm,self->asize);
            inner_map_init(self);
            outer_map_init(self);
            regressive_map_init(self);

            inverted_map_init(&self->product[ProductType_geometricinverted],&self->product[ProductType_geometric]);
            inverted_map_init(&self->product[ProductType_innerinverted],&self->product[ProductType_inner]);
            inverted_map_init(&self->product[ProductType_outerinverted],&self->product[ProductType_outer]);
            inverted_map_init(&self->product[ProductType_regressiveinverted],&self->product[ProductType_regressive]);

            self->tsize = 0;
            self->types = NULL;
            self->mdefault.type_name = NULL;

            break;
        default:
            return -1;
    }
    // set ga print type
    self->print_type = PrintType_metric_array;
    if(print_type < PrintTypeMAX && print_type > PrintTypeMIN)
        self->print_type = print_type;

    // set multivector print type
    self->print_type_mv = PrintTypeMV_reduced;
    if(print_type_mv < PrintTypeMVMAX && print_type_mv > PrintTypeMVMIN)
        self->print_type_mv = print_type_mv;

    self->precision = 1e-6;

    return 0;
}



static PyObject *algebra_add_basis(PyAlgebraObject *self, PyObject *args, PyObject *kwds){
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

        metric = (char*)PyMem_RawMalloc((METRIC_SIZE(self) + p + q + r)*sizeof(char));
        for(Py_ssize_t i = 0; i < METRIC_SIZE(self); i++)
            metric[i] = self->metric[i];
        for(Py_ssize_t i = METRIC_SIZE(self); i < p + METRIC_SIZE(self); i++)
            metric[i] = 1;
        for(Py_ssize_t i = p + METRIC_SIZE(self); i < p + q + METRIC_SIZE(self); i++)
            metric[i] = -1;
        for(Py_ssize_t i = p + q + METRIC_SIZE(self); i < p + q + r + METRIC_SIZE(self); i++)
            metric[i] = 0;
        PyMem_RawFree(self->metric);
        self->metric = metric;
    }else{
        if(!PyList_Check(metric_obj)){
            PyErr_SetString(PyExc_TypeError,"metric must be a list");
            return NULL;
        }
        Py_ssize_t size = PyList_Size(metric_obj);
        metric = (char*)PyMem_RawMalloc((METRIC_SIZE(self) + size)*sizeof(char));

        for(Py_ssize_t i = 0; i < METRIC_SIZE(self); i++)
            metric[i] = self->metric[i];

        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *value = PyList_GetItem(metric_obj,i);
            if(!PyLong_Check(value)){
                PyErr_SetString(PyExc_TypeError,"items of the list metric must be of type int");
                PyMem_RawFree(metric);
                return NULL;
            }else {
                int overflow;
                metric[i + METRIC_SIZE(self)] = (char)PyLong_AsLongAndOverflow(value,&overflow);
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

    map_add_basis(&self->product[ProductType_geometric],metric,METRIC_SIZE(self) + p + q + r, METRIC_SIZE(self));
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

    inverted_map_init(&self->product[ProductType_geometricinverted],&self->product[ProductType_geometric]);
    inverted_map_init(&self->product[ProductType_innerinverted],&self->product[ProductType_inner]);
    inverted_map_init(&self->product[ProductType_outerinverted],&self->product[ProductType_outer]);

    self->p += p; self->q += q; self->r += r;

    Py_RETURN_NONE;
}

static PyObject *algebra_repr(PyAlgebraObject *self){
    char  str[100];
    if(self->print_type == PrintType_metric){
        PyOS_snprintf(str,100,"GA(p=%lu,q=%lu,r=%lu)",self->p,self->q,self->r);
        return Py_BuildValue("s",str);
    }else if(self->print_type == PrintType_metric_array){
        PyObject *out;
        Py_ssize_t size = 3*METRIC_SIZE(self) + 2;
        char *str_temp = (char*)PyMem_RawMalloc(size*sizeof(char));
        char **str_value = (char**)PyMem_RawMalloc(METRIC_SIZE(self)*sizeof(char*));
        Py_ssize_t i;
        for(i = 0; i < METRIC_SIZE(self)-1; i++){
            str_value[i] = (char*)PyMem_RawMalloc(4*sizeof(char));
            PyOS_snprintf(str_value[i],4,"%+d,",(int)self->metric[i]);
        }
        str_value[i] = (char*)PyMem_RawMalloc(4*sizeof(char));
        PyOS_snprintf(str_value[i],4,"%+d",(int)self->metric[i]);

        Py_ssize_t j = 0;
        for(i = 0; i < METRIC_SIZE(self); i++){
            strcpy(str_temp + j,str_value[i]);
            j += strlen(str_value[i]);
        }
        str_temp[j] = '\0';
        PyOS_snprintf(str,100,"GA(metric=[%s])",str_temp);
        out = Py_BuildValue("s",str);

        for(Py_ssize_t i = 0; i < METRIC_SIZE(self); i++){
            PyMem_RawFree(str_value[i]);
        }PyMem_RawFree(str_value);
        PyMem_RawFree(str_temp);

        return out;
    }else if(self->print_type == PrintType_vectors){
        PyObject *out;
        char str_temp[100];
        char **str_value = (char**)PyMem_RawMalloc(METRIC_SIZE(self)*sizeof(char*));
        Py_ssize_t i;
        for(i = 0; i < METRIC_SIZE(self)-1; i++){
            str_value[i] = (char*)PyMem_RawMalloc(7*sizeof(char));
            PyOS_snprintf(str_value[i],7,"e%lu:%+d,",i,(int)self->metric[i]);
        }
        str_value[i] = (char*)PyMem_RawMalloc(7*sizeof(char));
        PyOS_snprintf(str_value[i],7,"e%lu:%+d",i,(int)self->metric[i]);

        Py_ssize_t j = 0;
        for(i = 0; i < METRIC_SIZE(self); i++){
            strcpy(str_temp + j, str_value[i]);
            j += strlen(str_value[i]);
        }
        str_temp[j] = '\0';
        PyOS_snprintf(str,100,"GA(vectors=[%s])",str_temp);
        out = Py_BuildValue("s",str);

        for(Py_ssize_t i = 0; i < METRIC_SIZE(self); i++){
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
            if(!IS_NONZERO(blade_str[j]))
                return -1;
            bitmap_i += 1 << (int)(blade_str[j] - '1');
        }
        (*bitmap)[i] = bitmap_i;
    }
    return size;
}

static PyObject* algebra_metric(PyAlgebraObject *self, PyObject *Py_UNUSED(ignored)){
    Py_ssize_t size = METRIC_SIZE(self);
    PyObject *metric_list = PyList_New(size);

    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *metrici = PyLong_FromLong((int)self->metric[i]);
        PyList_SetItem(metric_list,i,metrici);
    }

    return metric_list;
}


static PyObject *algebra_dualmap(PyAlgebraObject *self, PyObject *Py_UNUSED(ignored)){
    Py_ssize_t size = self->dm.size;
    PyObject *sign_list = PyList_New(size);
    PyObject *bitmap_list = PyList_New(size);
    PyObject *tuple = PyTuple_New(2);
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *signi = PyLong_FromLong(self->dm.sign[i]);
        PyObject *bitmapi = PyLong_FromLong(self->dm.bitmap[i]);
        PyList_SetItem(sign_list,i,signi);
        PyList_SetItem(bitmap_list,i,bitmapi);
    }
    PyTuple_SetItem(tuple,0,bitmap_list);
    PyTuple_SetItem(tuple,1,sign_list);

    return tuple;
}


static PyObject* algebra_grademap(PyAlgebraObject *self, PyObject *Py_UNUSED(ignored)){
    Py_ssize_t size = self->gm.size;
    PyObject *grade_list = PyList_New(size);
    PyObject *position_list = PyList_New(size);
    PyObject *gradesize_list = PyList_New(self->gm.max_grade+1);
    PyObject *tuple = PyTuple_New(3);


    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *gradei = PyLong_FromLong(self->gm.grade[i]);
        PyObject *positioni = PyLong_FromLong(self->gm.position[i]);
        PyList_SetItem(grade_list,i,gradei);
        PyList_SetItem(position_list,i,positioni);
    }
    for(Py_ssize_t i = 0; i < self->gm.max_grade+1; i++){
        PyObject *gradesizei = PyLong_FromLong(self->gm.grade_size[i]);
        PyList_SetItem(gradesize_list,i,gradesizei);
    }

    PyTuple_SetItem(tuple,0,grade_list);
    PyTuple_SetItem(tuple,1,position_list);
    PyTuple_SetItem(tuple,2,gradesize_list);
    return tuple;
}

static PyObject* algebra_cayley_table(PyAlgebraObject *self, PyObject *args){
    ProductType type = ProductType_geometricinverted;
    CliffordMap m = self->product[type];
    Py_ssize_t algebra_size = m.size;
    PyObject *sign_list = PyList_New(algebra_size);
    PyObject *bitmap_list = PyList_New(algebra_size);
    PyObject *tuple = PyTuple_New(2);

    char *type_str = NULL;

    PyArg_ParseTuple(args,"|s",&type_str);

    if(type_str){
        if(!strcmp("geometric",type_str)){
            type = ProductType_geometric;
        }else if(!strcmp("inner",type_str)){
            type = ProductType_inner;
        }else if(!strcmp("outer",type_str)){
            type = ProductType_outer;
        }else if(!strcmp("regressive",type_str)){
            type = ProductType_regressive;
        }else if(!strcmp("geometricinverted",type_str)){
            type = ProductType_geometricinverted;
        }else if(!strcmp("innerinverted",type_str)){
            type = ProductType_innerinverted;
        }else if(!strcmp("outerinverted",type_str)){
            type = ProductType_outerinverted;
        }else if(!strcmp("regressiveinverted",type_str)){
            type = ProductType_regressiveinverted;
        }
    }
    m = self->product[type];

    for(Py_ssize_t i = 0; i < algebra_size; i++){
        PyObject *bitmap_sublist = PyList_New(algebra_size);
        PyObject *sign_sublist = PyList_New(algebra_size);
        for(Py_ssize_t j = 0; j < algebra_size; j++){
            PyObject *bitmapij = PyLong_FromLong((long)m.bitmap[i][j]);
            PyObject *signij = PyLong_FromLong((long)m.sign[i][j]);
            PyList_SetItem(bitmap_sublist,j,bitmapij);
            PyList_SetItem(sign_sublist,j,signij);
        }
        PyList_SetItem(bitmap_list,i,bitmap_sublist);
        PyList_SetItem(sign_list,i,sign_sublist);
    }

    PyTuple_SetItem(tuple,0,bitmap_list);
    PyTuple_SetItem(tuple,1,sign_list);
    return tuple;
}


static PyNumberMethods PyMultivectorNumberMethods = {
    .nb_multiply = (binaryfunc) multivector_geometric_product,
    .nb_xor = (binaryfunc) multivector_outer_product,
    .nb_and = (binaryfunc) multivector_regressive_product,
    .nb_or = (binaryfunc) multivector_inner_product,
    .nb_add = (binaryfunc) multivector_add,
    .nb_subtract = (binaryfunc) multivector_subtract,
    .nb_invert = (unaryfunc) multivector_invert,
    .nb_negative = (unaryfunc) multivector_negative,
    .nb_positive = (unaryfunc) multivector_positive,

};


PyDoc_STRVAR(add_doc, "adds a bunch of multivectors.");
PyDoc_STRVAR(dual_doc, "dualizes the multivector.");
PyDoc_STRVAR(undual_doc, "undualizes the multivector.");
PyDoc_STRVAR(product_doc, "multiplies a bunch of multivectors.");
PyDoc_STRVAR(exponential_doc, "takes the exponential of multivectors.");
PyDoc_STRVAR(list_doc, "Returns a list with each coefficient of the multivector.");

PyMethodDef multivector_methods[] = {
    {"dual", (PyCFunction)multivector_dual, METH_NOARGS, dual_doc},
    {"undual", (PyCFunction)multivector_undual, METH_NOARGS, undual_doc},
    {"add", (PyCFunction) multivector_atomic_add, METH_VARARGS|METH_CLASS, add_doc},
    {"geometric_product", (PyCFunction) multivector_atomic_geometric_product, METH_VARARGS|METH_CLASS, product_doc},
    {"outer_product", (PyCFunction) multivector_atomic_outer_product, METH_VARARGS|METH_CLASS, product_doc},
    {"exp", (PyCFunction) multivector_exponential, METH_VARARGS|METH_CLASS, exponential_doc},
    {"list", (PyCFunction)multivector_list, METH_NOARGS, list_doc},
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

static PyMultivectorObject *populate_multivector_types(PyAlgebraObject *self){
    Py_ssize_t type = -1;
    PyMultivectorObject *multivector = NULL;
    char *dtype_str = self->mdefault.type_name;
    if(self->types == NULL){
        PyErr_SetString(PyExc_TypeError,"the operation table is empty");
        return NULL;
    }
    if(dtype_str){
        for(Py_ssize_t i = 0; i < self->tsize; i++){
            if(!strncmp(self->types[i].type_name,dtype_str,strlen(dtype_str))){
                type = i; break;
            }
        }
    }else type = 0;

    if(type == -1){
        PyErr_SetString(PyExc_TypeError,"type not found");
        return NULL; // raise error couldn't find the type
    }

    multivector = (PyMultivectorObject*)PyMem_RawMalloc(sizeof(PyMultivectorObject));
    multivector->type = self->types[type];
    multivector->mixed = self->mixed;

    Py_SET_TYPE(multivector,&PyMultivectorType);
    Py_XINCREF(&PyMultivectorType);
    multivector->GA = self;
    Py_XINCREF((PyObject*)self);
    Py_SET_REFCNT((PyObject*)multivector,1);

    return multivector;
}

static Py_ssize_t parse_list_as_values(PyObject *values, ga_float **values_float){
    if(!PyList_Check(values)){
        PyErr_SetString(PyExc_TypeError,"values must be a list");
        return -1;
    }
    Py_ssize_t size = PyList_Size(values);
    if(size <= 0) return -1;
    *values_float = (ga_float*)PyMem_RawMalloc(size*sizeof(ga_float));
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *value_i = PyList_GetItem(values,i);
        if(PyFloat_Check(value_i))
            (*values_float)[i] = (ga_float)PyFloat_AsDouble(value_i);
        else if(PyLong_Check(value_i))
            (*values_float)[i] = (ga_float)PyLong_AsLong(value_i);
        else{
            PyErr_SetString(PyExc_TypeError,"Elements of the list of values must be ga_float");
            PyMem_RawFree(*values_float);
            return -1;
        }
    }
    return size;
}

static PyObject *algebra_set_multivector_defaults(PyAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"dtype","precision",NULL};
    char *dtype = NULL;
    double precision = 1e-12;

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|f", kwlist,&dtype,&precision))
        return NULL;

    if(!dtype)
        return NULL;

    PyMem_RawFree(self->mdefault.type_name); // free previous
    self->mdefault.type_name = (char*)PyMem_RawMalloc((strlen(dtype)+1)*sizeof(char));
    strcpy(self->mdefault.type_name,dtype);
    self->precision = (ga_float)precision;
    Py_RETURN_NONE;
}


static PyObject *algebra_multivector(PyAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"values","blades",NULL};
    PyObject *values = NULL, *blades = NULL;
    int *bitmaps_int = NULL;
    ga_float *values_float = NULL;
    Py_ssize_t size,bsize;
    PyMultivectorObject *multivector;

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &values,&blades))
        return NULL;
    if(!values || !blades)
        return NULL;

    size = parse_list_as_values(values,&values_float);
    if(size <= 0){
        PyErr_SetString(PyExc_TypeError,"values must be a non empty list of integers or floats");
        PyMem_RawFree(values_float);
        return NULL;
    }
    bsize = parse_list_as_bitmaps(blades,&bitmaps_int);
    if(bsize != size){
        PyMem_RawFree(values_float);
        PyMem_RawFree(bitmaps_int);
        PyErr_SetString(PyExc_TypeError,"blades must be of the same size as values");
        return NULL;
    }

    multivector = populate_multivector_types(self);
    if(!multivector){
        PyMem_RawFree(values_float);
        PyMem_RawFree(bitmaps_int);
        return NULL;
    }

    gainitfunc init = multivector->type.data_funcs->init;
    if(init)
        multivector->data = init(bitmaps_int,values_float,size,self);
    else{
        PyMem_RawFree(values_float);
        PyMem_RawFree(bitmaps_int);
        return NULL; // raise not implemented error
    }

    PyMem_RawFree(values_float);
    PyMem_RawFree(bitmaps_int);

    return (PyObject*)multivector;
}

static PyObject *algebra_blades(PyAlgebraObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"blades","grades",NULL};
    PyObject *grades = NULL, *blades = NULL;
    int *bitmap = NULL;
    int *grade_array = NULL;
    Py_ssize_t *grade_bool = NULL;
    int **bitmap_array = NULL;
    char **bitmap_char = NULL;
    ga_float **value_array = NULL;
    Py_ssize_t size,gsize;
    PyMultivectorObject **multivectors = NULL;
    PyObject *dict_blades = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,&blades,&grades))
        return NULL;

    if(blades && grades)
        return NULL; // raise error

    if(blades){
        size = parse_list_as_bitmaps(blades,&bitmap);
        value_array = (ga_float**)PyMem_RawMalloc(size*sizeof(ga_float*));
        bitmap_array = (int**)PyMem_RawMalloc(size*sizeof(int*));
        for(Py_ssize_t i = 0; i < size; i++){
            value_array[i] = (ga_float*)PyMem_RawMalloc(sizeof(ga_float));
            bitmap_array[i] = (int*)PyMem_RawMalloc(sizeof(int));
            *(value_array[i]) = 1;
            *(bitmap_array[i]) = bitmap[i];
        }
    }else if(grades){
        gsize = parse_list_as_grades(self,grades,&grade_array);
        if(gsize <= 0) return NULL;
        grade_bool = get_grade_bool(grade_array,gsize,MAX_GRADE(self)+1);
        size = self->asize;
        Py_ssize_t psize = 0;
        for(Py_ssize_t i = 0; i < size; i++){
            if(grade_bool[GRADE(i)])
                psize++;
        }
        value_array = (ga_float**)PyMem_RawMalloc(psize*sizeof(ga_float*));
        bitmap_array = (int**)PyMem_RawMalloc(psize*sizeof(int*));
        Py_ssize_t j = 0;
        for(Py_ssize_t i = 0; i < size; i++){
            if(grade_bool[GRADE(i)] && j < psize){
                value_array[j] = (ga_float*)PyMem_RawMalloc(sizeof(ga_float));
                bitmap_array[j] = (int*)PyMem_RawMalloc(sizeof(int));
                *(value_array[j]) = 1;
                *(bitmap_array[j]) = i;
                j++;
            }else if(j>psize){
                break;
            }
        }
        size = psize;
    }else{
        size = self->asize;
        value_array = (ga_float**)PyMem_RawMalloc(size*sizeof(ga_float*));
        bitmap_array = (int**)PyMem_RawMalloc(size*sizeof(int*));
        for(Py_ssize_t i = 0; i < size; i++){
            value_array[i] = (ga_float*)PyMem_RawMalloc(sizeof(ga_float));
            bitmap_array[i] = (int*)PyMem_RawMalloc(sizeof(int));
            *(value_array[i]) = 1;
            *(bitmap_array[i]) = i;
        }
    }

    multivectors = (PyMultivectorObject**)PyMem_RawMalloc(size*sizeof(PyMultivectorObject*));
    if(!multivectors) goto fail;
    for(Py_ssize_t i = 0; i < size; i++)
        multivectors[i] = NULL;
    for(Py_ssize_t i = 0; i < size; i++){
        multivectors[i] = populate_multivector_types(self);
        if(!multivectors[i])
            goto fail;

    }

    bitmap_char = (char**)PyMem_RawMalloc(size*sizeof(char*));
    if(!bitmap_char) goto fail;
    for(Py_ssize_t i = 0; i < size; i++){
        bitmap_char[i] = bitmap_to_string(*(bitmap_array[i]));
        if(!bitmap_char[i])
            goto fail;
    }

    gainitfunc init = (*multivectors)->type.data_funcs->init;
    if(init)
        for(Py_ssize_t i = 0; i < size; i++)
            multivectors[i]->data = init(bitmap_array[i],value_array[i],1,self);
    else
        goto fail;


    dict_blades = PyDict_New();
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *key = Py_BuildValue("s",bitmap_char[i]);
        PyDict_SetItem(dict_blades,key,(PyObject*)multivectors[i]);
        Py_XDECREF(key);
        Py_XDECREF((PyObject*)multivectors[i]);
    }

    goto success;

fail:
    if(multivectors)
        for(Py_ssize_t i = 0; i < size; i++)
            Py_XDECREF(multivectors[i]);

success:
    if(bitmap_array)
        for(Py_ssize_t i = 0; i < size; i++)
            PyMem_RawFree(bitmap_array[i]);
    if(value_array)
        for(Py_ssize_t i = 0; i < size; i++)
            PyMem_RawFree(value_array[i]);
    if(bitmap_char)
        for(Py_ssize_t i = 0; i < size; i++)
            PyMem_RawFree(bitmap_char[i]);

    PyMem_RawFree(bitmap_array);
    PyMem_RawFree(bitmap_char);
    PyMem_RawFree(value_array);
    PyMem_RawFree(grade_bool);
    PyMem_RawFree(bitmap);
    PyMem_RawFree(multivectors);

    return dict_blades;
}


static PyMethodDef algebra_methods[] = {
    {"metric", (PyCFunction)algebra_metric, METH_NOARGS,
     "returns the metric array of the algebra" },
    {"dualmap", (PyCFunction)algebra_dualmap, METH_NOARGS,
     "returns the signs, and bitmaps of the dual algebra" },
    {"grademap", (PyCFunction)algebra_grademap, METH_NOARGS,
     "returns the grades, positions and grade sizes of the algebra" },
    {"cayley", (PyCFunction)algebra_cayley_table, METH_VARARGS,
     "returns the signs and bitmaps of the cayley table" },
    {"add_basis", (PyCFunction) algebra_add_basis, METH_VARARGS | METH_KEYWORDS,
     "adds basis vectors to the algebra" },
    {"multivector",(PyCFunction) algebra_multivector, METH_VARARGS | METH_KEYWORDS,
     "generate a multivector" },
    {"blades",(PyCFunction) algebra_blades, METH_VARARGS | METH_KEYWORDS,
     "generate blades for the algebra" },
    {"default",(PyCFunction) algebra_set_multivector_defaults, METH_VARARGS | METH_KEYWORDS,
     "set the default types of the multivector" },
    {NULL}  /* Sentinel */
};


static PyTypeObject PyGeometricAlgebraType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gasparse.GA",
    .tp_doc = PyDoc_STR("Construction of a geometric algebra given a metric"),
    .tp_basicsize = sizeof(PyAlgebraObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = (newfunc) algebra_new,
    .tp_init = (initproc) algebra_init,
    .tp_dealloc = (destructor) algebra_dealloc,
    .tp_repr = (reprfunc) algebra_repr,
    .tp_str = (reprfunc) algebra_repr,
    .tp_methods = algebra_methods,
};


PyDoc_STRVAR(gasparse_doc, "Implementation of sparse multivectors.");


static PyModuleDef gasparse_module = {
    PyModuleDef_HEAD_INIT,
#ifdef INCLUDE_GENCODE
    .m_name = "gasparsegen",
#else
    .m_name = "gasparse",
#endif
    .m_doc = gasparse_doc,
    .m_size = -1,
};

#ifdef INCLUDE_GENCODE
PyMODINIT_FUNC PyInit_gasparsegen(void){
#else
PyMODINIT_FUNC PyInit_gasparse(void){
#endif
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

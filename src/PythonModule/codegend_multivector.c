
/*
static *int get_grades_from_tuple(PyObject *tuple, Py_ssize_t *size, int max_grade){
    *size = PyTuple_Size(tuple);
    int *grades = (int*)PyMem_RawMalloc((*size)*sizeof(int));
    for(Py_ssize_t i = 0; i < *size; i++){
        PyObject *grade_obj = PyTuple_GetItem(tuple,i);
        if(!PyLong_Check(grade_obj))
            return NULL; // free grades and raise type error
        grades[i] = (int)PyLong_AsLong(grade_obj);
        if(grades[i] > max_grade)
            return NULL; // free grades and raise value error
    }
    return grades;
} */



typedef struct _gen0_DenseMultivector{
    ga_float value[8];
}gen0_DenseMultivector;

typedef struct _gen0_BladesMultivector{
    ga_float value0[1];
    ga_float value1[3];
    ga_float value2[3];
    ga_float value3[1];
}gen0_BladesMultivector;

static void gen0_dense_grade0project(gen0_DenseMultivector *dense0, gen0_DenseMultivector *dense){

    dense->value[0] = dense0->value[0];
}

static void gen0_dense_grade1project(gen0_DenseMultivector *dense0, gen0_DenseMultivector *dense){

    dense->value[1] = dense0->value[1];
    dense->value[2] = dense0->value[2];
    dense->value[4] = dense0->value[4];
}

static void gen0_dense_grade2project(gen0_DenseMultivector *dense0, gen0_DenseMultivector *dense){

    dense->value[3] = dense0->value[3];
    dense->value[5] = dense0->value[5];
    dense->value[6] = dense0->value[6];
}

static void gen0_dense_grade3project(gen0_DenseMultivector *dense0, gen0_DenseMultivector *dense){

    dense->value[7] = dense0->value[7];
}


static void gen0_blades_grade0project(gen0_BladesMultivector *blades0, gen0_BladesMultivector *blades){
    memcpy(blades->value0,blades0->value0,1);
}

static void gen0_blades_grade1project(gen0_BladesMultivector *blades0, gen0_BladesMultivector *blades){
    memcpy(blades->value1,blades0->value1,3);
}

static void gen0_blades_grade2project(gen0_BladesMultivector *blades0, gen0_BladesMultivector *blades){
    memcpy(blades->value2,blades0->value2,3);
}

static void gen0_blades_grade3project(gen0_BladesMultivector *blades0, gen0_BladesMultivector *blades){
    memcpy(blades->value3,blades0->value3,1);
}


#define GEN0_BLADES_GRADE0PROJECT(blades0,blades)\
    (memcpy(blades->value0,blades0->value0,1))


#define GEN0_BLADES_GRADE1PROJECT(blades0,blades)\
    (memcpy(blades->value1,blades0->value1,3))


#define GEN0_BLADES_GRADE2PROJECT(blades0,blades)\
    (memcpy(blades->value2,blades0->value2,3))


#define GEN0_BLADES_GRADE3PROJECT(blades0,blades)\
    (memcpy(blades->value3,blades0->value3,1))




typedef void (*gen0densegradeprojectfunc)(gen0_DenseMultivector*,gen0_DenseMultivector*);
typedef struct gen0_DenseGradeProject_funcs{
    gen0densegradeprojectfunc gradeproject[4];
}gen0_DenseGradeProject_func;

static gen0_DenseGradeProject_func gen0denseproject = {
    .gradeproject[0] = gen0_dense_grade0project,
    .gradeproject[1] = gen0_dense_grade1project,
    .gradeproject[2] = gen0_dense_grade2project,
    .gradeproject[3] = gen0_dense_grade3project,
};

typedef void (*gen0bladesgradeprojectfunc)(gen0_BladesMultivector*,gen0_BladesMultivector*);
typedef struct gen0_BladesGradeProject_funcs{
    gen0bladesgradeprojectfunc gradeproject[4];
}gen0_BladesGradeProject_func;

static gen0_BladesGradeProject_func gen0bladesproject = {
    .gradeproject[0] = gen0_blades_grade0project,
    .gradeproject[1] = gen0_blades_grade1project,
    .gradeproject[2] = gen0_blades_grade2project,
    .gradeproject[3] = gen0_blades_grade3project,
};


typedef struct _gen0_GradesBitmap{
        ga_float bitmap0[1];
        ga_float bitmap1[3];
        ga_float bitmap2[3];
        ga_float bitmap3[1];
    }gen0_GradesBitmap;

gen0_GradesBitmap gen0_gradesbitmap = {
    {0},
    {1,2,4},
    {3,5,6},
    {7},
}


static int gen0_blades_iter_next(PyMultivectorIter *iter){
    gen0_BladesMultivector *blades = (gen0_BladesMultivector*)iter->data;


    switch(*iter->index){
        case 0:
            iter->value = blades.value0[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap0[iter->index[1]];
            iter->grade = 0;
            if(++iter->index[1] >= 1){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 1:
            iter->value = blades.value1[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap1[iter->index[1]];
            iter->grade = 1;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 2:
            iter->value = blades.value2[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap2[iter->index[1]];
            iter->grade = 2;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 3:
            iter->value = blades.value3[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap3[iter->index[1]];
            iter->grade = 3;
            if(++iter->index[1] >= 1){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        default: // reset indices
            iter->index[1] = 0;
            iter->index[0] = 0;
            return 0; // stop loop
    }

    return 0;
}

static int gen0_dense_iter_next(PyMultivectorIter *iter){
    gen0_DenseMultivector *dense = (gen0_DenseMultivector*)iter->data;
    if(*iter->index >= 8){
        *iter->index = 0;
        return 0;
    }
    iter->bitmap = *iter->index;
    iter->value = dense->value[(*iter->index)++];
    iter->grade = GRADE(iter->bitmap);
    return 1;
}


static PyMultivectorObject *gen0_dense_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;
    gen0_DenseMultivector dense = {{0}};
    gen0_DenseMultivector *pdense;
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;

    for(Py_ssize_t i = 0; i < size; i++){
        gen0densegradeprojectfunc gradeproject =
                    gen0bladesproject.gradeproject[grades[i]];
        if(gradeproject != NULL)
            gradeproject(pdense0,&dense);
        else
            return NULL; // raise not implemented error
    }
    // should allocate all the necessary memory and also set the reference count to 1
    out = init_multivector(self);
    pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    *pdense = dense;
    out->data = (void*)pdense;
    PyMem_RawFree(grades);
    return (PyObject*)out;
}

static PyMultivectorObject *gen0_blades_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;

    gen0_BladesMultivector blades = {{0},{0},{0},{0},}
    gen0_BladesMultivector *pblades;
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;

    for(Py_ssize_t i = 0; i < size; i++){
        gen0bladesgradeprojectfunc gradeproject =
                    gen0bladesproject.gradeproject[grades[i]];
        if(gradeproject != NULL)
            gradeproject(pblades0,&blades);
        else
            return NULL; // raise not implemented error
    }
    // should allocate all the necessary memory and also set the reference count to 1
    out = init_multivector(self);
    pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    *pblades = blades;
    out->data = (void*)pblades;
    return (PyObject*)out;
}




static gen0_DenseMultivector gen0_dense_geometricproduct(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1){
    gen0_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[1]*dense1.value[1]
    +dense0.value[2]*dense1.value[2]
    -dense0.value[3]*dense1.value[3]
    +dense0.value[4]*dense1.value[4]
    -dense0.value[5]*dense1.value[5]
    -dense0.value[6]*dense1.value[6]
    -dense0.value[7]*dense1.value[7]
;
    dense.value[1] =
    -dense0.value[2]*dense1.value[3]
    +dense0.value[3]*dense1.value[2]
    -dense0.value[4]*dense1.value[5]
    +dense0.value[5]*dense1.value[4]
    -dense0.value[6]*dense1.value[7]
    -dense0.value[7]*dense1.value[6]
;
    dense.value[2] =
    +dense0.value[1]*dense1.value[3]
    -dense0.value[3]*dense1.value[1]
    -dense0.value[4]*dense1.value[6]
    +dense0.value[5]*dense1.value[7]
    +dense0.value[6]*dense1.value[4]
    +dense0.value[7]*dense1.value[5]
;
    dense.value[3] =
    +dense0.value[4]*dense1.value[7]
    +dense0.value[7]*dense1.value[4]
;
    dense.value[4] =
    +dense0.value[1]*dense1.value[5]
    +dense0.value[2]*dense1.value[6]
    -dense0.value[3]*dense1.value[7]
    -dense0.value[5]*dense1.value[1]
    -dense0.value[6]*dense1.value[2]
    -dense0.value[7]*dense1.value[3]
;
    dense.value[5] =
    -dense0.value[2]*dense1.value[7]
    -dense0.value[7]*dense1.value[2]
;
    dense.value[6] =
    +dense0.value[1]*dense1.value[7]
    +dense0.value[7]*dense1.value[1]
;
    return dense;
}

static gen0_DenseMultivector gen0_dense_add(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1, int sign){
    gen0_DenseMultivector dense = {{0}};
    if(sign == -1){
        dense.value[0] = dense0.value[0] - dense1.value[0];
        dense.value[1] = dense0.value[1] - dense1.value[1];
        dense.value[2] = dense0.value[2] - dense1.value[2];
        dense.value[3] = dense0.value[3] - dense1.value[3];
        dense.value[4] = dense0.value[4] - dense1.value[4];
        dense.value[5] = dense0.value[5] - dense1.value[5];
        dense.value[6] = dense0.value[6] - dense1.value[6];
        dense.value[7] = dense0.value[7] - dense1.value[7];
    }else if(sign == 1){
        dense.value[0] = dense0.value[0] + dense1.value[0];
        dense.value[1] = dense0.value[1] + dense1.value[1];
        dense.value[2] = dense0.value[2] + dense1.value[2];
        dense.value[3] = dense0.value[3] + dense1.value[3];
        dense.value[4] = dense0.value[4] + dense1.value[4];
        dense.value[5] = dense0.value[5] + dense1.value[5];
        dense.value[6] = dense0.value[6] + dense1.value[6];
        dense.value[7] = dense0.value[7] + dense1.value[7];
    } else{
        dense.value[0] = dense0.value[0] + sign*dense1.value[0];
        dense.value[1] = dense0.value[1] + sign*dense1.value[1];
        dense.value[2] = dense0.value[2] + sign*dense1.value[2];
        dense.value[3] = dense0.value[3] + sign*dense1.value[3];
        dense.value[4] = dense0.value[4] + sign*dense1.value[4];
        dense.value[5] = dense0.value[5] + sign*dense1.value[5];
        dense.value[6] = dense0.value[6] + sign*dense1.value[6];
        dense.value[7] = dense0.value[7] + sign*dense1.value[7];
    }
    return dense;
}


static gen0_DenseMultivector gen0_dense_scalaradd(gen1_DenseMultivector dense0, ga_float value, int sign){
    gen0_DenseMultivector dense = {{0}};
    if(sign == -1){
        dense.value[0] = -dense0.value[0];
        dense.value[1] = -dense0.value[1];
        dense.value[2] = -dense0.value[2];
        dense.value[3] = -dense0.value[3];
        dense.value[4] = -dense0.value[4];
        dense.value[5] = -dense0.value[5];
        dense.value[6] = -dense0.value[6];
        dense.value[7] = -dense0.value[7];
    }else if(sign == 1){
        dense.value[0] = dense0.value[0];
        dense.value[1] = dense0.value[1];
        dense.value[2] = dense0.value[2];
        dense.value[3] = dense0.value[3];
        dense.value[4] = dense0.value[4];
        dense.value[5] = dense0.value[5];
        dense.value[6] = dense0.value[6];
        dense.value[7] = dense0.value[7];
    } else{
        dense.value[0] = sign*dense0.value[0];
        dense.value[1] = sign*dense0.value[1];
        dense.value[2] = sign*dense0.value[2];
        dense.value[3] = sign*dense0.value[3];
        dense.value[4] = sign*dense0.value[4];
        dense.value[5] = sign*dense0.value[5];
        dense.value[6] = sign*dense0.value[6];
        dense.value[7] = sign*dense0.value[7];
    }
    dense.value[0] += value;
    return dense;
}

static gen0_DenseMultivector gen0_dense_scalarproducr(gen1_DenseMultivector dense0, ga_float value){
    gen0_DenseMultivector dense = {{0}};

    dense.value[0] = value*dense0.value[0];
    dense.value[1] = value*dense0.value[1];
    dense.value[2] = value*dense0.value[2];
    dense.value[3] = value*dense0.value[3];
    dense.value[4] = value*dense0.value[4];
    dense.value[5] = value*dense0.value[5];
    dense.value[6] = value*dense0.value[6];
    dense.value[7] = value*dense0.value[7];

    return dense;
}


PyMultivectorObject *gen0_dense_multivector_geometricproduct(PyMultivectorObject *data0, PyMultivectorObject *data1){
    gen0_DenseMultivector *dense0 = (gen0_DenseMultivector*)data0->data;
    gen0_DenseMultivector *dense1 = (gen0_DenseMultivector*)data1->data;
    gen0_DenseMultivector *dense  = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!pdense0 || !pdense1 || !pdense){
        PyMem_RawFree(pdense);
        return NULL; // raise error
    }

    PyMultivectorObject *out = new_multivector(data0);
    *dense = gen0_dense_geometricproduct(dense0,dense1);
    out->data = (void*)dense;
    return out;
}

static gen0_BladesMultivector gen0_blades_geometricproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},}

    blades.value0[0] =
    +blades0.value1[0]*blades1.value1[0]
    +blades0.value1[1]*blades1.value1[1]
    -blades0.value2[0]*blades1.value2[0]
    +blades0.value1[2]*blades1.value1[2]
    -blades0.value2[1]*blades1.value2[1]
    -blades0.value2[2]*blades1.value2[2]
    -blades0.value3[0]*blades1.value3[0]
;
    blades.value1[0] =
    -blades0.value1[1]*blades1.value2[0]
    +blades0.value2[0]*blades1.value1[1]
    -blades0.value1[2]*blades1.value2[1]
    +blades0.value2[1]*blades1.value1[2]
    -blades0.value2[2]*blades1.value3[0]
    -blades0.value3[0]*blades1.value2[2]
;
    blades.value1[1] =
    +blades0.value1[0]*blades1.value2[0]
    -blades0.value2[0]*blades1.value1[0]
    -blades0.value1[2]*blades1.value2[2]
    +blades0.value2[1]*blades1.value3[0]
    +blades0.value2[2]*blades1.value1[2]
    +blades0.value3[0]*blades1.value2[1]
;
    blades.value2[0] =
    +blades0.value1[2]*blades1.value3[0]
    +blades0.value3[0]*blades1.value1[2]
;
    blades.value1[2] =
    +blades0.value1[0]*blades1.value2[1]
    +blades0.value1[1]*blades1.value2[2]
    -blades0.value2[0]*blades1.value3[0]
    -blades0.value2[1]*blades1.value1[0]
    -blades0.value2[2]*blades1.value1[1]
    -blades0.value3[0]*blades1.value2[0]
;
    blades.value2[1] =
    -blades0.value1[1]*blades1.value3[0]
    -blades0.value3[0]*blades1.value1[1]
;
    blades.value2[2] =
    +blades0.value1[0]*blades1.value3[0]
    +blades0.value3[0]*blades1.value1[0]
;
    return blades;
}

static gen0_BladesMultivector gen0_blades_add(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1, int sign){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},}

    if(sign == -1){
        blades.value0[0] = blades0.value0[0] - blades1.value0[0];
        blades.value1[0] = blades0.value1[0] - blades1.value1[0];
        blades.value1[1] = blades0.value1[1] - blades1.value1[1];
        blades.value1[2] = blades0.value1[2] - blades1.value1[2];
        blades.value2[0] = blades0.value2[0] - blades1.value2[0];
        blades.value2[1] = blades0.value2[1] - blades1.value2[1];
        blades.value2[2] = blades0.value2[2] - blades1.value2[2];
        blades.value3[0] = blades0.value3[0] - blades1.value3[0];
    }else if(sign == 1){
        blades.value0[0] = blades0.value0[0] + blades1.value0[0];
        blades.value1[0] = blades0.value1[0] + blades1.value1[0];
        blades.value1[1] = blades0.value1[1] + blades1.value1[1];
        blades.value1[2] = blades0.value1[2] + blades1.value1[2];
        blades.value2[0] = blades0.value2[0] + blades1.value2[0];
        blades.value2[1] = blades0.value2[1] + blades1.value2[1];
        blades.value2[2] = blades0.value2[2] + blades1.value2[2];
        blades.value3[0] = blades0.value3[0] + blades1.value3[0];
    }else{
        blades.value0[0] = blades0.value0[0] + sign*blades1.value0[0];
        blades.value1[0] = blades0.value1[0] + sign*blades1.value1[0];
        blades.value1[1] = blades0.value1[1] + sign*blades1.value1[1];
        blades.value1[2] = blades0.value1[2] + sign*blades1.value1[2];
        blades.value2[0] = blades0.value2[0] + sign*blades1.value2[0];
        blades.value2[1] = blades0.value2[1] + sign*blades1.value2[1];
        blades.value2[2] = blades0.value2[2] + sign*blades1.value2[2];
        blades.value3[0] = blades0.value3[0] + sign*blades1.value3[0];
    }
    return blades;
}


static gen0_BladesMultivector gen0_blades_scalaradd(gen0_BladesMultivector blades0, ga_float value, int sign){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},}

    if(sign == -1){
        blades.value0[0] = -blades0.value0[0];
        blades.value1[0] = -blades0.value1[0];
        blades.value1[1] = -blades0.value1[1];
        blades.value1[2] = -blades0.value1[2];
        blades.value2[0] = -blades0.value2[0];
        blades.value2[1] = -blades0.value2[1];
        blades.value2[2] = -blades0.value2[2];
        blades.value3[0] = -blades0.value3[0];
    }else if(sign == 1){
        blades.value0[0] = blades0.value0[0];
        blades.value1[0] = blades0.value1[0];
        blades.value1[1] = blades0.value1[1];
        blades.value1[2] = blades0.value1[2];
        blades.value2[0] = blades0.value2[0];
        blades.value2[1] = blades0.value2[1];
        blades.value2[2] = blades0.value2[2];
        blades.value3[0] = blades0.value3[0];
    }else{
        blades.value0[0] = sign*blades1.value0[0];
        blades.value1[0] = sign*blades1.value1[0];
        blades.value1[1] = sign*blades1.value1[1];
        blades.value1[2] = sign*blades1.value1[2];
        blades.value2[0] = sign*blades1.value2[0];
        blades.value2[1] = sign*blades1.value2[1];
        blades.value2[2] = sign*blades1.value2[2];
        blades.value3[0] = sign*blades1.value3[0];
    }
    blades.value0[0] += value;
    return blades;
}


static gen0_BladesMultivector gen0_blades_scalarproduct(gen0_BladesMultivector blades0, ga_float value){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},}

    blades.value0[0] = value*blades0.value0[0];
    blades.value1[0] = value*blades0.value1[0];
    blades.value1[1] = value*blades0.value1[1];
    blades.value1[2] = value*blades0.value1[2];
    blades.value2[0] = value*blades0.value2[0];
    blades.value2[1] = value*blades0.value2[1];
    blades.value2[2] = value*blades0.value2[2];
    blades.value3[0] = value*blades0.value3[0];
    return blades;
}



PyMultivectorObject *gen0_blades_multivector_geometricproduct(PyMultivectorObject *data0, PyMultivectorObject *data1){
    gen0_BladesMultivector *blades0 = (gen0_BladesMultivector*)data0->data;
    gen0_BladesMultivector *blades1 = (gen0_BladesMultivector*)data1->data;
    gen0_BladesMultivector *blades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!pblades0 || !pblades1 || !pblades){
        PyMem_RawFree(pblades);
        return NULL; // raise error
    }

    PyMultivectorObject *out = new_multivector(data0);

    *blades = gen0_blades_geometricproduct(blades0,blades1);
    out->data = (void*)blades;
    return out;
}




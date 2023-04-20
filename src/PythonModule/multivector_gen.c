#include <Python.h>
#include "gasparse.h"
#include "multivector.h"
#include "multivector_gen.h"

typedef struct _gen0_DenseMultivector{
    ga_float value[8];
}gen0_DenseMultivector;

typedef struct _gen0_BladesMultivector{
    ga_float value0[1];
    ga_float value1[3];
    ga_float value2[3];
    ga_float value3[1];
}gen0_BladesMultivector;


typedef struct _gen0_GradesBitmap{
    int bitmap0[1];
    int bitmap1[3];
    int bitmap2[3];
    int bitmap3[1];
}gen0_GradesBitmap;

static gen0_GradesBitmap gen0_gradesbitmap = {
    {0},
    {1,2,4},
    {3,5,6},
    {7},
};

static int gen0_grades_position[8] = {0,0,1,0,2,1,2,0};

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

static gen0_DenseMultivector dense0_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen0_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= 8){
            return dense; // raise error
        }
        dense.value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return dense;
}

static gen0_BladesMultivector blades0_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    for(Py_ssize_t i = 0; i < size; i++){
        switch(GRADE(bitmap[i])){
            case 0:
                blades.value0[gen0_grades_position[bitmap[i]]] += value[i];
                break;
            case 1:
                blades.value1[gen0_grades_position[bitmap[i]]] += value[i];
                break;
            case 2:
                blades.value2[gen0_grades_position[bitmap[i]]] += value[i];
                break;
            case 3:
                blades.value3[gen0_grades_position[bitmap[i]]] += value[i];
                break;
            default:
                return blades; // raise error
        }
    }
    return blades;
}

static PyMultivectorIter dense0_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data = data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type; // don't know what is the type of this
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = 8;
    iter.next = data->data_funcs.iter_next[data->type];
    return iter;
}

static PyMultivectorIter blades0_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data= data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(2*sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.index[1] = 0;
    iter.size = 2;
    iter.niters = 8;
    iter.next = data->data_funcs.iter_next[data->type];
    return iter;
}


static int blades0_iternext(PyMultivectorIter *iter){
    gen0_BladesMultivector *blades = (gen0_BladesMultivector*)iter->data;


    switch(*iter->index){
        case 0:
            iter->value = blades->value0[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap0[iter->index[1]];
            iter->grade = 0;
            if(++iter->index[1] >= 1){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 1:
            iter->value = blades->value1[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap1[iter->index[1]];
            iter->grade = 1;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 2:
            iter->value = blades->value2[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap2[iter->index[1]];
            iter->grade = 2;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 3:
            iter->value = blades->value3[iter->index[1]];
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

static int dense0_iternext(PyMultivectorIter *iter){
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

typedef struct _gen1_DenseMultivector{
    ga_float value[32];
}gen1_DenseMultivector;

typedef struct _gen1_BladesMultivector{
    ga_float value0[1];
    ga_float value1[5];
    ga_float value2[10];
    ga_float value3[10];
    ga_float value4[5];
    ga_float value5[1];
}gen1_BladesMultivector;


typedef struct _gen1_GradesBitmap{
    int bitmap0[1];
    int bitmap1[5];
    int bitmap2[10];
    int bitmap3[10];
    int bitmap4[5];
    int bitmap5[1];
}gen1_GradesBitmap;

static gen1_GradesBitmap gen1_gradesbitmap = {
    {0},
    {1,2,4,8,16},
    {3,5,6,9,10,12,17,18,20,24},
    {7,11,13,14,19,21,22,25,26,28},
    {15,23,27,29,30},
    {31},
};

static int gen1_grades_position[32] = {0,0,1,0,2,1,2,0,3,3,4,1,5,2,3,0,4,6,7,4,8,5,6,1,9,7,8,2,9,3,4,0};

static void gen1_dense_grade0project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[0] = dense0->value[0];
}

static void gen1_dense_grade1project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[1] = dense0->value[1];
    dense->value[2] = dense0->value[2];
    dense->value[4] = dense0->value[4];
    dense->value[8] = dense0->value[8];
    dense->value[16] = dense0->value[16];
}

static void gen1_dense_grade2project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[3] = dense0->value[3];
    dense->value[5] = dense0->value[5];
    dense->value[6] = dense0->value[6];
    dense->value[9] = dense0->value[9];
    dense->value[10] = dense0->value[10];
    dense->value[12] = dense0->value[12];
    dense->value[17] = dense0->value[17];
    dense->value[18] = dense0->value[18];
    dense->value[20] = dense0->value[20];
    dense->value[24] = dense0->value[24];
}

static void gen1_dense_grade3project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[7] = dense0->value[7];
    dense->value[11] = dense0->value[11];
    dense->value[13] = dense0->value[13];
    dense->value[14] = dense0->value[14];
    dense->value[19] = dense0->value[19];
    dense->value[21] = dense0->value[21];
    dense->value[22] = dense0->value[22];
    dense->value[25] = dense0->value[25];
    dense->value[26] = dense0->value[26];
    dense->value[28] = dense0->value[28];
}

static void gen1_dense_grade4project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[15] = dense0->value[15];
    dense->value[23] = dense0->value[23];
    dense->value[27] = dense0->value[27];
    dense->value[29] = dense0->value[29];
    dense->value[30] = dense0->value[30];
}

static void gen1_dense_grade5project(gen1_DenseMultivector *dense0, gen1_DenseMultivector *dense){

    dense->value[31] = dense0->value[31];
}


static void gen1_blades_grade0project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value0,blades0->value0,1);
}

static void gen1_blades_grade1project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value1,blades0->value1,5);
}

static void gen1_blades_grade2project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value2,blades0->value2,10);
}

static void gen1_blades_grade3project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value3,blades0->value3,10);
}

static void gen1_blades_grade4project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value4,blades0->value4,5);
}

static void gen1_blades_grade5project(gen1_BladesMultivector *blades0, gen1_BladesMultivector *blades){
    memcpy(blades->value5,blades0->value5,1);
}


#define GEN1_BLADES_GRADE0PROJECT(blades0,blades)\
    (memcpy(blades->value0,blades0->value0,1))


#define GEN1_BLADES_GRADE1PROJECT(blades0,blades)\
    (memcpy(blades->value1,blades0->value1,5))


#define GEN1_BLADES_GRADE2PROJECT(blades0,blades)\
    (memcpy(blades->value2,blades0->value2,10))


#define GEN1_BLADES_GRADE3PROJECT(blades0,blades)\
    (memcpy(blades->value3,blades0->value3,10))


#define GEN1_BLADES_GRADE4PROJECT(blades0,blades)\
    (memcpy(blades->value4,blades0->value4,5))


#define GEN1_BLADES_GRADE5PROJECT(blades0,blades)\
    (memcpy(blades->value5,blades0->value5,1))




typedef void (*gen1densegradeprojectfunc)(gen1_DenseMultivector*,gen1_DenseMultivector*);

typedef struct gen1_DenseGradeProject_funcs{
    gen1densegradeprojectfunc gradeproject[6];
}gen1_DenseGradeProject_func;

static gen1_DenseGradeProject_func gen1denseproject = {
    .gradeproject[0] = gen1_dense_grade0project,
    .gradeproject[1] = gen1_dense_grade1project,
    .gradeproject[2] = gen1_dense_grade2project,
    .gradeproject[3] = gen1_dense_grade3project,
    .gradeproject[4] = gen1_dense_grade4project,
    .gradeproject[5] = gen1_dense_grade5project,
};

typedef void (*gen1bladesgradeprojectfunc)(gen1_BladesMultivector*,gen1_BladesMultivector*);
typedef struct gen1_BladesGradeProject_funcs{
    gen1bladesgradeprojectfunc gradeproject[6];
}gen1_BladesGradeProject_func;

static gen1_BladesGradeProject_func gen1bladesproject = {
    .gradeproject[0] = gen1_blades_grade0project,
    .gradeproject[1] = gen1_blades_grade1project,
    .gradeproject[2] = gen1_blades_grade2project,
    .gradeproject[3] = gen1_blades_grade3project,
    .gradeproject[4] = gen1_blades_grade4project,
    .gradeproject[5] = gen1_blades_grade5project,
};

static gen1_DenseMultivector dense1_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen1_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= 32){
            return dense; // raise error
        }
        dense.value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return dense;
}

static gen1_BladesMultivector blades1_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    for(Py_ssize_t i = 0; i < size; i++){
        switch(GRADE(bitmap[i])){
            case 0:
                blades.value0[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            case 1:
                blades.value1[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            case 2:
                blades.value2[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            case 3:
                blades.value3[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            case 4:
                blades.value4[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            case 5:
                blades.value5[gen1_grades_position[bitmap[i]]] += value[i];
                break;
            default:
                return blades; // raise error
        }
    }
    return blades;
}

static PyMultivectorIter dense1_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data = data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type; // don't know what is the type of this
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = 32;
    iter.next = data->data_funcs.iter_next[data->type];
    return iter;
}

static PyMultivectorIter blades1_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data= data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(2*sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.index[1] = 0;
    iter.size = 2;
    iter.niters = 32;
    iter.next = data->data_funcs.iter_next[data->type];
    return iter;
}


static int blades1_iternext(PyMultivectorIter *iter){
    gen1_BladesMultivector *blades = (gen1_BladesMultivector*)iter->data;


    switch(*iter->index){
        case 0:
            iter->value = blades->value0[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap0[iter->index[1]];
            iter->grade = 0;
            if(++iter->index[1] >= 1){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 1:
            iter->value = blades->value1[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap1[iter->index[1]];
            iter->grade = 1;
            if(++iter->index[1] >= 5){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 2:
            iter->value = blades->value2[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap2[iter->index[1]];
            iter->grade = 2;
            if(++iter->index[1] >= 10){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 3:
            iter->value = blades->value3[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap3[iter->index[1]];
            iter->grade = 3;
            if(++iter->index[1] >= 10){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 4:
            iter->value = blades->value4[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap4[iter->index[1]];
            iter->grade = 4;
            if(++iter->index[1] >= 5){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            break;
        case 5:
            iter->value = blades->value5[iter->index[1]];
            iter->bitmap = gen1_gradesbitmap.bitmap5[iter->index[1]];
            iter->grade = 5;
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

static int dense1_iternext(PyMultivectorIter *iter){
    gen1_DenseMultivector *dense = (gen1_DenseMultivector*)iter->data;
    if(*iter->index >= 32){
        *iter->index = 0;
        return 0;
    }
    iter->bitmap = *iter->index;
    iter->value = dense->value[(*iter->index)++];
    iter->grade = GRADE(iter->bitmap);
    return 1;
}




static gen0_DenseMultivector gen0_dense_geometricproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1){
    gen0_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[0]*dense1.value[0]
    +dense0.value[1]*dense1.value[1]
    +dense0.value[2]*dense1.value[2]
    -dense0.value[3]*dense1.value[3]
    +dense0.value[4]*dense1.value[4]
    -dense0.value[5]*dense1.value[5]
    -dense0.value[6]*dense1.value[6]
    -dense0.value[7]*dense1.value[7]
;
    dense.value[1] =
    +dense0.value[0]*dense1.value[1]
    +dense0.value[1]*dense1.value[0]
    -dense0.value[2]*dense1.value[3]
    +dense0.value[3]*dense1.value[2]
    -dense0.value[4]*dense1.value[5]
    +dense0.value[5]*dense1.value[4]
    -dense0.value[6]*dense1.value[7]
    -dense0.value[7]*dense1.value[6]
;
    dense.value[2] =
    +dense0.value[0]*dense1.value[2]
    +dense0.value[1]*dense1.value[3]
    +dense0.value[2]*dense1.value[0]
    -dense0.value[3]*dense1.value[1]
    -dense0.value[4]*dense1.value[6]
    +dense0.value[5]*dense1.value[7]
    +dense0.value[6]*dense1.value[4]
    +dense0.value[7]*dense1.value[5]
;
    dense.value[3] =
    +dense0.value[0]*dense1.value[3]
    +dense0.value[1]*dense1.value[2]
    -dense0.value[2]*dense1.value[1]
    +dense0.value[3]*dense1.value[0]
    +dense0.value[4]*dense1.value[7]
    -dense0.value[5]*dense1.value[6]
    +dense0.value[6]*dense1.value[5]
    +dense0.value[7]*dense1.value[4]
;
    dense.value[4] =
    +dense0.value[0]*dense1.value[4]
    +dense0.value[1]*dense1.value[5]
    +dense0.value[2]*dense1.value[6]
    -dense0.value[3]*dense1.value[7]
    +dense0.value[4]*dense1.value[0]
    -dense0.value[5]*dense1.value[1]
    -dense0.value[6]*dense1.value[2]
    -dense0.value[7]*dense1.value[3]
;
    dense.value[5] =
    +dense0.value[0]*dense1.value[5]
    +dense0.value[1]*dense1.value[4]
    -dense0.value[2]*dense1.value[7]
    +dense0.value[3]*dense1.value[6]
    -dense0.value[4]*dense1.value[1]
    +dense0.value[5]*dense1.value[0]
    -dense0.value[6]*dense1.value[3]
    -dense0.value[7]*dense1.value[2]
;
    dense.value[6] =
    +dense0.value[0]*dense1.value[6]
    +dense0.value[1]*dense1.value[7]
    +dense0.value[2]*dense1.value[4]
    -dense0.value[3]*dense1.value[5]
    -dense0.value[4]*dense1.value[2]
    +dense0.value[5]*dense1.value[3]
    +dense0.value[6]*dense1.value[0]
    +dense0.value[7]*dense1.value[1]
;
    dense.value[7] =
    +dense0.value[0]*dense1.value[7]
    +dense0.value[1]*dense1.value[6]
    -dense0.value[2]*dense1.value[5]
    +dense0.value[3]*dense1.value[4]
    +dense0.value[4]*dense1.value[3]
    -dense0.value[5]*dense1.value[2]
    +dense0.value[6]*dense1.value[1]
    +dense0.value[7]*dense1.value[0]
;
    return dense;
}


static gen0_DenseMultivector gen0_dense_innerproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1){
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


static gen0_DenseMultivector gen0_dense_outerproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1){
    gen0_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[0]*dense1.value[0]
;
    dense.value[1] =
    +dense0.value[0]*dense1.value[1]
    +dense0.value[1]*dense1.value[0]
;
    dense.value[2] =
    +dense0.value[0]*dense1.value[2]
    +dense0.value[2]*dense1.value[0]
;
    dense.value[3] =
    +dense0.value[0]*dense1.value[3]
    +dense0.value[1]*dense1.value[2]
    -dense0.value[2]*dense1.value[1]
    +dense0.value[3]*dense1.value[0]
;
    dense.value[4] =
    +dense0.value[0]*dense1.value[4]
    +dense0.value[4]*dense1.value[0]
;
    dense.value[5] =
    +dense0.value[0]*dense1.value[5]
    +dense0.value[1]*dense1.value[4]
    -dense0.value[4]*dense1.value[1]
    +dense0.value[5]*dense1.value[0]
;
    dense.value[6] =
    +dense0.value[0]*dense1.value[6]
    +dense0.value[2]*dense1.value[4]
    -dense0.value[4]*dense1.value[2]
    +dense0.value[6]*dense1.value[0]
;
    dense.value[7] =
    +dense0.value[0]*dense1.value[7]
    +dense0.value[1]*dense1.value[6]
    -dense0.value[2]*dense1.value[5]
    +dense0.value[3]*dense1.value[4]
    +dense0.value[4]*dense1.value[3]
    -dense0.value[5]*dense1.value[2]
    +dense0.value[6]*dense1.value[1]
    +dense0.value[7]*dense1.value[0]
;
    return dense;
}


static gen0_DenseMultivector gen0_dense_atomicadd(gen0_DenseMultivector *dense_array, Py_ssize_t size){
    gen0_DenseMultivector dense = {{0}};

    for(Py_ssize_t i = 0; i < size; i++){
        dense.value[0] += dense_array[i].value[0];
        dense.value[1] += dense_array[i].value[1];
        dense.value[2] += dense_array[i].value[2];
        dense.value[3] += dense_array[i].value[3];
        dense.value[4] += dense_array[i].value[4];
        dense.value[5] += dense_array[i].value[5];
        dense.value[6] += dense_array[i].value[6];
        dense.value[7] += dense_array[i].value[7];
    }

    return dense;
}




static gen0_DenseMultivector gen0_dense_add(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1, int sign){
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


static gen0_DenseMultivector gen0_dense_scalaradd(gen0_DenseMultivector dense0, ga_float value, int sign){
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

static gen0_DenseMultivector gen0_dense_scalarproduct(gen0_DenseMultivector dense0, ga_float value){
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

static gen0_DenseMultivector gen0_dense_reverse(gen0_DenseMultivector dense0){
    gen0_DenseMultivector dense = {{0}};

    dense.value[0] = dense0.value[0];
    dense.value[1] = dense0.value[1];
    dense.value[2] = dense0.value[2];
    dense.value[3] = -dense0.value[3];
    dense.value[4] = dense0.value[4];
    dense.value[5] = -dense0.value[5];
    dense.value[6] = -dense0.value[6];
    dense.value[7] = -dense0.value[7];

    return dense;
}



static gen0_BladesMultivector gen0_blades_geometricproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    blades.value0[0] =
    +blades0.value0[0]*blades1.value0[0]
    +blades0.value1[0]*blades1.value1[0]
    +blades0.value1[1]*blades1.value1[1]
    -blades0.value2[0]*blades1.value2[0]
    +blades0.value1[2]*blades1.value1[2]
    -blades0.value2[1]*blades1.value2[1]
    -blades0.value2[2]*blades1.value2[2]
    -blades0.value3[0]*blades1.value3[0]
;
    blades.value1[0] =
    +blades0.value0[0]*blades1.value1[0]
    +blades0.value1[0]*blades1.value0[0]
    -blades0.value1[1]*blades1.value2[0]
    +blades0.value2[0]*blades1.value1[1]
    -blades0.value1[2]*blades1.value2[1]
    +blades0.value2[1]*blades1.value1[2]
    -blades0.value2[2]*blades1.value3[0]
    -blades0.value3[0]*blades1.value2[2]
;
    blades.value1[1] =
    +blades0.value0[0]*blades1.value1[1]
    +blades0.value1[0]*blades1.value2[0]
    +blades0.value1[1]*blades1.value0[0]
    -blades0.value2[0]*blades1.value1[0]
    -blades0.value1[2]*blades1.value2[2]
    +blades0.value2[1]*blades1.value3[0]
    +blades0.value2[2]*blades1.value1[2]
    +blades0.value3[0]*blades1.value2[1]
;
    blades.value2[0] =
    +blades0.value0[0]*blades1.value2[0]
    +blades0.value1[0]*blades1.value1[1]
    -blades0.value1[1]*blades1.value1[0]
    +blades0.value2[0]*blades1.value0[0]
    +blades0.value1[2]*blades1.value3[0]
    -blades0.value2[1]*blades1.value2[2]
    +blades0.value2[2]*blades1.value2[1]
    +blades0.value3[0]*blades1.value1[2]
;
    blades.value1[2] =
    +blades0.value0[0]*blades1.value1[2]
    +blades0.value1[0]*blades1.value2[1]
    +blades0.value1[1]*blades1.value2[2]
    -blades0.value2[0]*blades1.value3[0]
    +blades0.value1[2]*blades1.value0[0]
    -blades0.value2[1]*blades1.value1[0]
    -blades0.value2[2]*blades1.value1[1]
    -blades0.value3[0]*blades1.value2[0]
;
    blades.value2[1] =
    +blades0.value0[0]*blades1.value2[1]
    +blades0.value1[0]*blades1.value1[2]
    -blades0.value1[1]*blades1.value3[0]
    +blades0.value2[0]*blades1.value2[2]
    -blades0.value1[2]*blades1.value1[0]
    +blades0.value2[1]*blades1.value0[0]
    -blades0.value2[2]*blades1.value2[0]
    -blades0.value3[0]*blades1.value1[1]
;
    blades.value2[2] =
    +blades0.value0[0]*blades1.value2[2]
    +blades0.value1[0]*blades1.value3[0]
    +blades0.value1[1]*blades1.value1[2]
    -blades0.value2[0]*blades1.value2[1]
    -blades0.value1[2]*blades1.value1[1]
    +blades0.value2[1]*blades1.value2[0]
    +blades0.value2[2]*blades1.value0[0]
    +blades0.value3[0]*blades1.value1[0]
;
    blades.value3[0] =
    +blades0.value0[0]*blades1.value3[0]
    +blades0.value1[0]*blades1.value2[2]
    -blades0.value1[1]*blades1.value2[1]
    +blades0.value2[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value2[0]
    -blades0.value2[1]*blades1.value1[1]
    +blades0.value2[2]*blades1.value1[0]
    +blades0.value3[0]*blades1.value0[0]
;
    return blades;
}


static gen0_BladesMultivector gen0_blades_innerproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

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


static gen0_BladesMultivector gen0_blades_outerproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    blades.value0[0] =
    +blades0.value0[0]*blades1.value0[0]
;
    blades.value1[0] =
    +blades0.value0[0]*blades1.value1[0]
    +blades0.value1[0]*blades1.value0[0]
;
    blades.value1[1] =
    +blades0.value0[0]*blades1.value1[1]
    +blades0.value1[1]*blades1.value0[0]
;
    blades.value2[0] =
    +blades0.value0[0]*blades1.value2[0]
    +blades0.value1[0]*blades1.value1[1]
    -blades0.value1[1]*blades1.value1[0]
    +blades0.value2[0]*blades1.value0[0]
;
    blades.value1[2] =
    +blades0.value0[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value0[0]
;
    blades.value2[1] =
    +blades0.value0[0]*blades1.value2[1]
    +blades0.value1[0]*blades1.value1[2]
    -blades0.value1[2]*blades1.value1[0]
    +blades0.value2[1]*blades1.value0[0]
;
    blades.value2[2] =
    +blades0.value0[0]*blades1.value2[2]
    +blades0.value1[1]*blades1.value1[2]
    -blades0.value1[2]*blades1.value1[1]
    +blades0.value2[2]*blades1.value0[0]
;
    blades.value3[0] =
    +blades0.value0[0]*blades1.value3[0]
    +blades0.value1[0]*blades1.value2[2]
    -blades0.value1[1]*blades1.value2[1]
    +blades0.value2[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value2[0]
    -blades0.value2[1]*blades1.value1[1]
    +blades0.value2[2]*blades1.value1[0]
    +blades0.value3[0]*blades1.value0[0]
;
    return blades;
}


static gen0_BladesMultivector gen0_blades_atomicadd(gen0_BladesMultivector *blades_array, Py_ssize_t size){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    for(Py_ssize_t i = 0; i < size; i++){
       blades.value0[0] += blades_array[i].value0[0];
       blades.value1[0] += blades_array[i].value1[0];
       blades.value1[1] += blades_array[i].value1[1];
       blades.value1[2] += blades_array[i].value1[2];
       blades.value2[0] += blades_array[i].value2[0];
       blades.value2[1] += blades_array[i].value2[1];
       blades.value2[2] += blades_array[i].value2[2];
       blades.value3[0] += blades_array[i].value3[0];
    }
    return blades;
}

static gen0_BladesMultivector gen0_blades_add(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1, int sign){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

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
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

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
        blades.value0[0] = sign*blades0.value0[0];
        blades.value1[0] = sign*blades0.value1[0];
        blades.value1[1] = sign*blades0.value1[1];
        blades.value1[2] = sign*blades0.value1[2];
        blades.value2[0] = sign*blades0.value2[0];
        blades.value2[1] = sign*blades0.value2[1];
        blades.value2[2] = sign*blades0.value2[2];
        blades.value3[0] = sign*blades0.value3[0];
    }
    blades.value0[0] += value;
    return blades;
}


static gen0_BladesMultivector gen0_blades_scalarproduct(gen0_BladesMultivector blades0, ga_float value){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

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

static gen0_BladesMultivector gen0_blades_reverse(gen0_BladesMultivector blades0){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    blades.value0[0] = blades0.value0[0];
    blades.value1[0] = blades0.value1[0];
    blades.value1[1] = blades0.value1[1];
    blades.value1[2] = blades0.value1[2];
    blades.value2[0] = -blades0.value2[0];
    blades.value2[1] = -blades0.value2[1];
    blades.value2[2] = -blades0.value2[2];
    blades.value3[0] = -blades0.value3[0];
    return blades;
}


static PyMultivectorObject *binary_dense0_product(PyMultivectorObject *data0, PyMultivectorObject *data1,ProductType ptype){
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)data0->data;
    gen0_DenseMultivector *pdense1 = (gen0_DenseMultivector*)data1->data;
    gen0_DenseMultivector *pdense  = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pdense0 || !pdense1 || !pdense || !out){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = gen0_dense_geometricproduct(*pdense0,*pdense1);
            break;
        case ProductType_inner:
            *pdense = gen0_dense_innerproduct(*pdense0,*pdense1);
            break;
        case ProductType_outer:
            *pdense = gen0_dense_outerproduct(*pdense0,*pdense1);
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades0_product(PyMultivectorObject *data0, PyMultivectorObject *data1,ProductType ptype){
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)data0->data;
    gen0_BladesMultivector *pblades1 = (gen0_BladesMultivector*)data1->data;
    gen0_BladesMultivector *pblades  = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pblades0 || !pblades1 || !pblades || !out){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = gen0_blades_geometricproduct(*pblades0,*pblades1);
            break;
        case ProductType_inner:
            *pblades = gen0_blades_innerproduct(*pblades0,*pblades1);
            break;
        case ProductType_outer:
            *pblades = gen0_blades_outerproduct(*pblades0,*pblades1);
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *ternary_dense0_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2,ProductType ptype){
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)data0->data;
    gen0_DenseMultivector *pdense1 = (gen0_DenseMultivector*)data1->data;
    gen0_DenseMultivector *pdense2 = (gen0_DenseMultivector*)data2->data;
    gen0_DenseMultivector *pdense  = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pdense0 || !pdense1 || !pdense2 || !pdense || !out){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = gen0_dense_geometricproduct(*pdense0,*pdense1);
            *pdense = gen0_dense_geometricproduct(*pdense,*pdense2);
            break;
        case ProductType_inner:
            *pdense = gen0_dense_innerproduct(*pdense0,*pdense1);
            *pdense = gen0_dense_innerproduct(*pdense,*pdense2);
            break;
        case ProductType_outer:
            *pdense = gen0_dense_outerproduct(*pdense0,*pdense1);
            *pdense = gen0_dense_outerproduct(*pdense,*pdense2);
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *ternary_blades0_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2,ProductType ptype){
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)data0->data;
    gen0_BladesMultivector *pblades1 = (gen0_BladesMultivector*)data1->data;
    gen0_BladesMultivector *pblades2 = (gen0_BladesMultivector*)data2->data;
    gen0_BladesMultivector *pblades  = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pblades0 || !pblades1 || !pblades2 || !pblades || !out){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = gen0_blades_geometricproduct(*pblades0,*pblades1);
            *pblades = gen0_blades_geometricproduct(*pblades,*pblades2);
            break;
        case ProductType_inner:
            *pblades = gen0_blades_innerproduct(*pblades0,*pblades1);
            *pblades = gen0_blades_innerproduct(*pblades,*pblades2);
            break;
        case ProductType_outer:
            *pblades = gen0_blades_outerproduct(*pblades0,*pblades1);
            *pblades = gen0_blades_outerproduct(*pblades,*pblades2);
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *unary_dense0_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;
    gen0_DenseMultivector dense = {{0}};
    gen0_DenseMultivector *pdense;
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;

    for(Py_ssize_t i = 0; i < size; i++){
        gen0densegradeprojectfunc gradeproject =
                    gen0denseproject.gradeproject[grades[i]];
        if(gradeproject != NULL)
            gradeproject(pdense0,&dense);
        else
            return NULL; // raise not implemented error
    }
    // should allocate all the necessary memory and also set the reference count to 1
    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    *pdense = dense;
    out->data = (void*)pdense;
    PyMem_RawFree(grades);
    return out;
}

static PyMultivectorObject *unary_blades0_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;

    gen0_BladesMultivector blades = {{0},{0},{0},{0},};
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
    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    *pblades = blades;
    out->data = (void*)pblades;
    return out;
}

static PyMultivectorObject* atomic_dense0_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorObject *out = new_multivector(data,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)PyMem_RawMalloc(size*sizeof(gen0_DenseMultivector));
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense0);
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    for(Py_ssize_t i = 0; i < size; i++)
        pdense0[i] = *((gen0_DenseMultivector*)data[i].data);

    *pdense = gen0_dense_atomicadd(pdense0,size);
    out->data = (void*)pdense;
    return out;
}
static PyMultivectorObject* atomic_blades0_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorObject *out = new_multivector(data,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)PyMem_RawMalloc(size*sizeof(gen0_BladesMultivector));
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades0);
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    for(Py_ssize_t i = 0; i < size; i++)
        pblades0[i] = *((gen0_BladesMultivector*)data[i].data);

    *pblades = gen0_blades_atomicadd(pblades0,size);
    out->data = (void*)pblades;
    return out;
}

static PyMultivectorObject* binary_dense0_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    PyMultivectorObject *out = new_multivector(data0,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)data0->data;
    gen0_DenseMultivector *pdense1 = (gen0_DenseMultivector*)data1->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense1 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_add(*pdense0,*pdense1,sign);
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject* binary_blades0_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    PyMultivectorObject *out = new_multivector(data0,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)data0->data;
    gen0_BladesMultivector *pblades1 = (gen0_BladesMultivector*)data1->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades1 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_add(*pblades0,*pblades1,sign);
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}


static PyMultivectorObject* atomic_dense0_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    if(size < 2) return NULL;
    PyMultivectorObject *out = new_multivector(data,-1);
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    gen0_DenseMultivector dense;
    if(!out  || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    switch(ptype){
        case ProductType_geometric:
            dense = gen0_dense_geometricproduct(
                      *((gen0_DenseMultivector*)data[0].data),
                      *((gen0_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen0_dense_geometricproduct(
                          dense,
                          *((gen0_DenseMultivector*)data[i].data));
            }
            break;
        case ProductType_inner:
            dense = gen0_dense_innerproduct(
                      *((gen0_DenseMultivector*)data[0].data),
                      *((gen0_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen0_dense_innerproduct(
                          dense,
                          *((gen0_DenseMultivector*)data[i].data));
            }
            break;
        case ProductType_outer:
            dense = gen0_dense_outerproduct(
                      *((gen0_DenseMultivector*)data[0].data),
                      *((gen0_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen0_dense_outerproduct(
                          dense,
                          *((gen0_DenseMultivector*)data[i].data));
            }
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }
    *pdense = dense;
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject* atomic_blades0_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    if(size < 2) return NULL;
    PyMultivectorObject *out = new_multivector(data,-1);
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    gen0_BladesMultivector blades;
    if(!out  || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    switch(ptype){
        case ProductType_geometric:
            blades = gen0_blades_geometricproduct(
                      *((gen0_BladesMultivector*)data[0].data),
                      *((gen0_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen0_blades_geometricproduct(
                          blades,
                          *((gen0_BladesMultivector*)data[i].data));
            }
            break;
        case ProductType_inner:
            blades = gen0_blades_innerproduct(
                      *((gen0_BladesMultivector*)data[0].data),
                      *((gen0_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen0_blades_innerproduct(
                          blades,
                          *((gen0_BladesMultivector*)data[i].data));
            }
            break;
        case ProductType_outer:
            blades = gen0_blades_outerproduct(
                      *((gen0_BladesMultivector*)data[0].data),
                      *((gen0_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen0_blades_outerproduct(
                          blades,
                          *((gen0_BladesMultivector*)data[i].data));
            }
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }
    *pblades = blades;
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_dense0_scalarproduct(PyMultivectorObject *self, ga_float value){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_scalarproduct(*pdense0,value); // multiply by a scalar
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades0_scalarproduct(PyMultivectorObject *self, ga_float value){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_scalarproduct(*pblades0,value); // multiply by a scalar
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_dense0_scalaradd(PyMultivectorObject *self, ga_float value, int sign){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_scalaradd(*pdense0,value,sign); // add a scalar
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades0_scalaradd(PyMultivectorObject *self, ga_float value, int sign){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_scalaradd(*pblades0,value,sign); // add a scalar
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *unary_dense0_reverse(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_reverse(*pdense0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *unary_blades0_reverse(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_reverse(*pblades0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static gen1_DenseMultivector gen1_dense_geometricproduct(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1){
    gen1_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[0]*dense1.value[0]
    +dense0.value[1]*dense1.value[1]
    +dense0.value[2]*dense1.value[2]
    -dense0.value[3]*dense1.value[3]
    +dense0.value[4]*dense1.value[4]
    -dense0.value[5]*dense1.value[5]
    -dense0.value[6]*dense1.value[6]
    -dense0.value[7]*dense1.value[7]
    +dense0.value[8]*dense1.value[8]
    -dense0.value[9]*dense1.value[9]
    -dense0.value[10]*dense1.value[10]
    -dense0.value[11]*dense1.value[11]
    -dense0.value[12]*dense1.value[12]
    -dense0.value[13]*dense1.value[13]
    -dense0.value[14]*dense1.value[14]
    +dense0.value[15]*dense1.value[15]
    -dense0.value[16]*dense1.value[16]
    +dense0.value[17]*dense1.value[17]
    +dense0.value[18]*dense1.value[18]
    +dense0.value[19]*dense1.value[19]
    +dense0.value[20]*dense1.value[20]
    +dense0.value[21]*dense1.value[21]
    +dense0.value[22]*dense1.value[22]
    -dense0.value[23]*dense1.value[23]
    +dense0.value[24]*dense1.value[24]
    +dense0.value[25]*dense1.value[25]
    +dense0.value[26]*dense1.value[26]
    -dense0.value[27]*dense1.value[27]
    +dense0.value[28]*dense1.value[28]
    -dense0.value[29]*dense1.value[29]
    -dense0.value[30]*dense1.value[30]
    -dense0.value[31]*dense1.value[31]
;
    dense.value[1] =
    +dense0.value[0]*dense1.value[1]
    +dense0.value[1]*dense1.value[0]
    -dense0.value[2]*dense1.value[3]
    +dense0.value[3]*dense1.value[2]
    -dense0.value[4]*dense1.value[5]
    +dense0.value[5]*dense1.value[4]
    -dense0.value[6]*dense1.value[7]
    -dense0.value[7]*dense1.value[6]
    -dense0.value[8]*dense1.value[9]
    +dense0.value[9]*dense1.value[8]
    -dense0.value[10]*dense1.value[11]
    -dense0.value[11]*dense1.value[10]
    -dense0.value[12]*dense1.value[13]
    -dense0.value[13]*dense1.value[12]
    +dense0.value[14]*dense1.value[15]
    -dense0.value[15]*dense1.value[14]
    +dense0.value[16]*dense1.value[17]
    -dense0.value[17]*dense1.value[16]
    +dense0.value[18]*dense1.value[19]
    +dense0.value[19]*dense1.value[18]
    +dense0.value[20]*dense1.value[21]
    +dense0.value[21]*dense1.value[20]
    -dense0.value[22]*dense1.value[23]
    +dense0.value[23]*dense1.value[22]
    +dense0.value[24]*dense1.value[25]
    +dense0.value[25]*dense1.value[24]
    -dense0.value[26]*dense1.value[27]
    +dense0.value[27]*dense1.value[26]
    -dense0.value[28]*dense1.value[29]
    +dense0.value[29]*dense1.value[28]
    -dense0.value[30]*dense1.value[31]
    -dense0.value[31]*dense1.value[30]
;
    dense.value[2] =
    +dense0.value[0]*dense1.value[2]
    +dense0.value[1]*dense1.value[3]
    +dense0.value[2]*dense1.value[0]
    -dense0.value[3]*dense1.value[1]
    -dense0.value[4]*dense1.value[6]
    +dense0.value[5]*dense1.value[7]
    +dense0.value[6]*dense1.value[4]
    +dense0.value[7]*dense1.value[5]
    -dense0.value[8]*dense1.value[10]
    +dense0.value[9]*dense1.value[11]
    +dense0.value[10]*dense1.value[8]
    +dense0.value[11]*dense1.value[9]
    -dense0.value[12]*dense1.value[14]
    -dense0.value[13]*dense1.value[15]
    -dense0.value[14]*dense1.value[12]
    +dense0.value[15]*dense1.value[13]
    +dense0.value[16]*dense1.value[18]
    -dense0.value[17]*dense1.value[19]
    -dense0.value[18]*dense1.value[16]
    -dense0.value[19]*dense1.value[17]
    +dense0.value[20]*dense1.value[22]
    +dense0.value[21]*dense1.value[23]
    +dense0.value[22]*dense1.value[20]
    -dense0.value[23]*dense1.value[21]
    +dense0.value[24]*dense1.value[26]
    +dense0.value[25]*dense1.value[27]
    +dense0.value[26]*dense1.value[24]
    -dense0.value[27]*dense1.value[25]
    -dense0.value[28]*dense1.value[30]
    +dense0.value[29]*dense1.value[31]
    +dense0.value[30]*dense1.value[28]
    +dense0.value[31]*dense1.value[29]
;
    dense.value[3] =
    +dense0.value[0]*dense1.value[3]
    +dense0.value[1]*dense1.value[2]
    -dense0.value[2]*dense1.value[1]
    +dense0.value[3]*dense1.value[0]
    +dense0.value[4]*dense1.value[7]
    -dense0.value[5]*dense1.value[6]
    +dense0.value[6]*dense1.value[5]
    +dense0.value[7]*dense1.value[4]
    +dense0.value[8]*dense1.value[11]
    -dense0.value[9]*dense1.value[10]
    +dense0.value[10]*dense1.value[9]
    +dense0.value[11]*dense1.value[8]
    -dense0.value[12]*dense1.value[15]
    -dense0.value[13]*dense1.value[14]
    +dense0.value[14]*dense1.value[13]
    -dense0.value[15]*dense1.value[12]
    -dense0.value[16]*dense1.value[19]
    +dense0.value[17]*dense1.value[18]
    -dense0.value[18]*dense1.value[17]
    -dense0.value[19]*dense1.value[16]
    +dense0.value[20]*dense1.value[23]
    +dense0.value[21]*dense1.value[22]
    -dense0.value[22]*dense1.value[21]
    +dense0.value[23]*dense1.value[20]
    +dense0.value[24]*dense1.value[27]
    +dense0.value[25]*dense1.value[26]
    -dense0.value[26]*dense1.value[25]
    +dense0.value[27]*dense1.value[24]
    +dense0.value[28]*dense1.value[31]
    -dense0.value[29]*dense1.value[30]
    +dense0.value[30]*dense1.value[29]
    +dense0.value[31]*dense1.value[28]
;
    dense.value[4] =
    +dense0.value[0]*dense1.value[4]
    +dense0.value[1]*dense1.value[5]
    +dense0.value[2]*dense1.value[6]
    -dense0.value[3]*dense1.value[7]
    +dense0.value[4]*dense1.value[0]
    -dense0.value[5]*dense1.value[1]
    -dense0.value[6]*dense1.value[2]
    -dense0.value[7]*dense1.value[3]
    -dense0.value[8]*dense1.value[12]
    +dense0.value[9]*dense1.value[13]
    +dense0.value[10]*dense1.value[14]
    +dense0.value[11]*dense1.value[15]
    +dense0.value[12]*dense1.value[8]
    +dense0.value[13]*dense1.value[9]
    +dense0.value[14]*dense1.value[10]
    -dense0.value[15]*dense1.value[11]
    +dense0.value[16]*dense1.value[20]
    -dense0.value[17]*dense1.value[21]
    -dense0.value[18]*dense1.value[22]
    -dense0.value[19]*dense1.value[23]
    -dense0.value[20]*dense1.value[16]
    -dense0.value[21]*dense1.value[17]
    -dense0.value[22]*dense1.value[18]
    +dense0.value[23]*dense1.value[19]
    +dense0.value[24]*dense1.value[28]
    +dense0.value[25]*dense1.value[29]
    +dense0.value[26]*dense1.value[30]
    -dense0.value[27]*dense1.value[31]
    +dense0.value[28]*dense1.value[24]
    -dense0.value[29]*dense1.value[25]
    -dense0.value[30]*dense1.value[26]
    -dense0.value[31]*dense1.value[27]
;
    dense.value[5] =
    +dense0.value[0]*dense1.value[5]
    +dense0.value[1]*dense1.value[4]
    -dense0.value[2]*dense1.value[7]
    +dense0.value[3]*dense1.value[6]
    -dense0.value[4]*dense1.value[1]
    +dense0.value[5]*dense1.value[0]
    -dense0.value[6]*dense1.value[3]
    -dense0.value[7]*dense1.value[2]
    +dense0.value[8]*dense1.value[13]
    -dense0.value[9]*dense1.value[12]
    +dense0.value[10]*dense1.value[15]
    +dense0.value[11]*dense1.value[14]
    +dense0.value[12]*dense1.value[9]
    +dense0.value[13]*dense1.value[8]
    -dense0.value[14]*dense1.value[11]
    +dense0.value[15]*dense1.value[10]
    -dense0.value[16]*dense1.value[21]
    +dense0.value[17]*dense1.value[20]
    -dense0.value[18]*dense1.value[23]
    -dense0.value[19]*dense1.value[22]
    -dense0.value[20]*dense1.value[17]
    -dense0.value[21]*dense1.value[16]
    +dense0.value[22]*dense1.value[19]
    -dense0.value[23]*dense1.value[18]
    +dense0.value[24]*dense1.value[29]
    +dense0.value[25]*dense1.value[28]
    -dense0.value[26]*dense1.value[31]
    +dense0.value[27]*dense1.value[30]
    -dense0.value[28]*dense1.value[25]
    +dense0.value[29]*dense1.value[24]
    -dense0.value[30]*dense1.value[27]
    -dense0.value[31]*dense1.value[26]
;
    dense.value[6] =
    +dense0.value[0]*dense1.value[6]
    +dense0.value[1]*dense1.value[7]
    +dense0.value[2]*dense1.value[4]
    -dense0.value[3]*dense1.value[5]
    -dense0.value[4]*dense1.value[2]
    +dense0.value[5]*dense1.value[3]
    +dense0.value[6]*dense1.value[0]
    +dense0.value[7]*dense1.value[1]
    +dense0.value[8]*dense1.value[14]
    -dense0.value[9]*dense1.value[15]
    -dense0.value[10]*dense1.value[12]
    -dense0.value[11]*dense1.value[13]
    +dense0.value[12]*dense1.value[10]
    +dense0.value[13]*dense1.value[11]
    +dense0.value[14]*dense1.value[8]
    -dense0.value[15]*dense1.value[9]
    -dense0.value[16]*dense1.value[22]
    +dense0.value[17]*dense1.value[23]
    +dense0.value[18]*dense1.value[20]
    +dense0.value[19]*dense1.value[21]
    -dense0.value[20]*dense1.value[18]
    -dense0.value[21]*dense1.value[19]
    -dense0.value[22]*dense1.value[16]
    +dense0.value[23]*dense1.value[17]
    +dense0.value[24]*dense1.value[30]
    +dense0.value[25]*dense1.value[31]
    +dense0.value[26]*dense1.value[28]
    -dense0.value[27]*dense1.value[29]
    -dense0.value[28]*dense1.value[26]
    +dense0.value[29]*dense1.value[27]
    +dense0.value[30]*dense1.value[24]
    +dense0.value[31]*dense1.value[25]
;
    dense.value[7] =
    +dense0.value[0]*dense1.value[7]
    +dense0.value[1]*dense1.value[6]
    -dense0.value[2]*dense1.value[5]
    +dense0.value[3]*dense1.value[4]
    +dense0.value[4]*dense1.value[3]
    -dense0.value[5]*dense1.value[2]
    +dense0.value[6]*dense1.value[1]
    +dense0.value[7]*dense1.value[0]
    -dense0.value[8]*dense1.value[15]
    +dense0.value[9]*dense1.value[14]
    -dense0.value[10]*dense1.value[13]
    -dense0.value[11]*dense1.value[12]
    +dense0.value[12]*dense1.value[11]
    +dense0.value[13]*dense1.value[10]
    -dense0.value[14]*dense1.value[9]
    +dense0.value[15]*dense1.value[8]
    +dense0.value[16]*dense1.value[23]
    -dense0.value[17]*dense1.value[22]
    +dense0.value[18]*dense1.value[21]
    +dense0.value[19]*dense1.value[20]
    -dense0.value[20]*dense1.value[19]
    -dense0.value[21]*dense1.value[18]
    +dense0.value[22]*dense1.value[17]
    -dense0.value[23]*dense1.value[16]
    +dense0.value[24]*dense1.value[31]
    +dense0.value[25]*dense1.value[30]
    -dense0.value[26]*dense1.value[29]
    +dense0.value[27]*dense1.value[28]
    +dense0.value[28]*dense1.value[27]
    -dense0.value[29]*dense1.value[26]
    +dense0.value[30]*dense1.value[25]
    +dense0.value[31]*dense1.value[24]
;
    dense.value[8] =
    +dense0.value[0]*dense1.value[8]
    +dense0.value[1]*dense1.value[9]
    +dense0.value[2]*dense1.value[10]
    -dense0.value[3]*dense1.value[11]
    +dense0.value[4]*dense1.value[12]
    -dense0.value[5]*dense1.value[13]
    -dense0.value[6]*dense1.value[14]
    -dense0.value[7]*dense1.value[15]
    +dense0.value[8]*dense1.value[0]
    -dense0.value[9]*dense1.value[1]
    -dense0.value[10]*dense1.value[2]
    -dense0.value[11]*dense1.value[3]
    -dense0.value[12]*dense1.value[4]
    -dense0.value[13]*dense1.value[5]
    -dense0.value[14]*dense1.value[6]
    +dense0.value[15]*dense1.value[7]
    +dense0.value[16]*dense1.value[24]
    -dense0.value[17]*dense1.value[25]
    -dense0.value[18]*dense1.value[26]
    -dense0.value[19]*dense1.value[27]
    -dense0.value[20]*dense1.value[28]
    -dense0.value[21]*dense1.value[29]
    -dense0.value[22]*dense1.value[30]
    +dense0.value[23]*dense1.value[31]
    -dense0.value[24]*dense1.value[16]
    -dense0.value[25]*dense1.value[17]
    -dense0.value[26]*dense1.value[18]
    +dense0.value[27]*dense1.value[19]
    -dense0.value[28]*dense1.value[20]
    +dense0.value[29]*dense1.value[21]
    +dense0.value[30]*dense1.value[22]
    +dense0.value[31]*dense1.value[23]
;
    dense.value[9] =
    +dense0.value[0]*dense1.value[9]
    +dense0.value[1]*dense1.value[8]
    -dense0.value[2]*dense1.value[11]
    +dense0.value[3]*dense1.value[10]
    -dense0.value[4]*dense1.value[13]
    +dense0.value[5]*dense1.value[12]
    -dense0.value[6]*dense1.value[15]
    -dense0.value[7]*dense1.value[14]
    -dense0.value[8]*dense1.value[1]
    +dense0.value[9]*dense1.value[0]
    -dense0.value[10]*dense1.value[3]
    -dense0.value[11]*dense1.value[2]
    -dense0.value[12]*dense1.value[5]
    -dense0.value[13]*dense1.value[4]
    +dense0.value[14]*dense1.value[7]
    -dense0.value[15]*dense1.value[6]
    -dense0.value[16]*dense1.value[25]
    +dense0.value[17]*dense1.value[24]
    -dense0.value[18]*dense1.value[27]
    -dense0.value[19]*dense1.value[26]
    -dense0.value[20]*dense1.value[29]
    -dense0.value[21]*dense1.value[28]
    +dense0.value[22]*dense1.value[31]
    -dense0.value[23]*dense1.value[30]
    -dense0.value[24]*dense1.value[17]
    -dense0.value[25]*dense1.value[16]
    +dense0.value[26]*dense1.value[19]
    -dense0.value[27]*dense1.value[18]
    +dense0.value[28]*dense1.value[21]
    -dense0.value[29]*dense1.value[20]
    +dense0.value[30]*dense1.value[23]
    +dense0.value[31]*dense1.value[22]
;
    dense.value[10] =
    +dense0.value[0]*dense1.value[10]
    +dense0.value[1]*dense1.value[11]
    +dense0.value[2]*dense1.value[8]
    -dense0.value[3]*dense1.value[9]
    -dense0.value[4]*dense1.value[14]
    +dense0.value[5]*dense1.value[15]
    +dense0.value[6]*dense1.value[12]
    +dense0.value[7]*dense1.value[13]
    -dense0.value[8]*dense1.value[2]
    +dense0.value[9]*dense1.value[3]
    +dense0.value[10]*dense1.value[0]
    +dense0.value[11]*dense1.value[1]
    -dense0.value[12]*dense1.value[6]
    -dense0.value[13]*dense1.value[7]
    -dense0.value[14]*dense1.value[4]
    +dense0.value[15]*dense1.value[5]
    -dense0.value[16]*dense1.value[26]
    +dense0.value[17]*dense1.value[27]
    +dense0.value[18]*dense1.value[24]
    +dense0.value[19]*dense1.value[25]
    -dense0.value[20]*dense1.value[30]
    -dense0.value[21]*dense1.value[31]
    -dense0.value[22]*dense1.value[28]
    +dense0.value[23]*dense1.value[29]
    -dense0.value[24]*dense1.value[18]
    -dense0.value[25]*dense1.value[19]
    -dense0.value[26]*dense1.value[16]
    +dense0.value[27]*dense1.value[17]
    +dense0.value[28]*dense1.value[22]
    -dense0.value[29]*dense1.value[23]
    -dense0.value[30]*dense1.value[20]
    -dense0.value[31]*dense1.value[21]
;
    dense.value[11] =
    +dense0.value[0]*dense1.value[11]
    +dense0.value[1]*dense1.value[10]
    -dense0.value[2]*dense1.value[9]
    +dense0.value[3]*dense1.value[8]
    +dense0.value[4]*dense1.value[15]
    -dense0.value[5]*dense1.value[14]
    +dense0.value[6]*dense1.value[13]
    +dense0.value[7]*dense1.value[12]
    +dense0.value[8]*dense1.value[3]
    -dense0.value[9]*dense1.value[2]
    +dense0.value[10]*dense1.value[1]
    +dense0.value[11]*dense1.value[0]
    -dense0.value[12]*dense1.value[7]
    -dense0.value[13]*dense1.value[6]
    +dense0.value[14]*dense1.value[5]
    -dense0.value[15]*dense1.value[4]
    +dense0.value[16]*dense1.value[27]
    -dense0.value[17]*dense1.value[26]
    +dense0.value[18]*dense1.value[25]
    +dense0.value[19]*dense1.value[24]
    -dense0.value[20]*dense1.value[31]
    -dense0.value[21]*dense1.value[30]
    +dense0.value[22]*dense1.value[29]
    -dense0.value[23]*dense1.value[28]
    -dense0.value[24]*dense1.value[19]
    -dense0.value[25]*dense1.value[18]
    +dense0.value[26]*dense1.value[17]
    -dense0.value[27]*dense1.value[16]
    -dense0.value[28]*dense1.value[23]
    +dense0.value[29]*dense1.value[22]
    -dense0.value[30]*dense1.value[21]
    -dense0.value[31]*dense1.value[20]
;
    dense.value[12] =
    +dense0.value[0]*dense1.value[12]
    +dense0.value[1]*dense1.value[13]
    +dense0.value[2]*dense1.value[14]
    -dense0.value[3]*dense1.value[15]
    +dense0.value[4]*dense1.value[8]
    -dense0.value[5]*dense1.value[9]
    -dense0.value[6]*dense1.value[10]
    -dense0.value[7]*dense1.value[11]
    -dense0.value[8]*dense1.value[4]
    +dense0.value[9]*dense1.value[5]
    +dense0.value[10]*dense1.value[6]
    +dense0.value[11]*dense1.value[7]
    +dense0.value[12]*dense1.value[0]
    +dense0.value[13]*dense1.value[1]
    +dense0.value[14]*dense1.value[2]
    -dense0.value[15]*dense1.value[3]
    -dense0.value[16]*dense1.value[28]
    +dense0.value[17]*dense1.value[29]
    +dense0.value[18]*dense1.value[30]
    +dense0.value[19]*dense1.value[31]
    +dense0.value[20]*dense1.value[24]
    +dense0.value[21]*dense1.value[25]
    +dense0.value[22]*dense1.value[26]
    -dense0.value[23]*dense1.value[27]
    -dense0.value[24]*dense1.value[20]
    -dense0.value[25]*dense1.value[21]
    -dense0.value[26]*dense1.value[22]
    +dense0.value[27]*dense1.value[23]
    -dense0.value[28]*dense1.value[16]
    +dense0.value[29]*dense1.value[17]
    +dense0.value[30]*dense1.value[18]
    +dense0.value[31]*dense1.value[19]
;
    dense.value[13] =
    +dense0.value[0]*dense1.value[13]
    +dense0.value[1]*dense1.value[12]
    -dense0.value[2]*dense1.value[15]
    +dense0.value[3]*dense1.value[14]
    -dense0.value[4]*dense1.value[9]
    +dense0.value[5]*dense1.value[8]
    -dense0.value[6]*dense1.value[11]
    -dense0.value[7]*dense1.value[10]
    +dense0.value[8]*dense1.value[5]
    -dense0.value[9]*dense1.value[4]
    +dense0.value[10]*dense1.value[7]
    +dense0.value[11]*dense1.value[6]
    +dense0.value[12]*dense1.value[1]
    +dense0.value[13]*dense1.value[0]
    -dense0.value[14]*dense1.value[3]
    +dense0.value[15]*dense1.value[2]
    +dense0.value[16]*dense1.value[29]
    -dense0.value[17]*dense1.value[28]
    +dense0.value[18]*dense1.value[31]
    +dense0.value[19]*dense1.value[30]
    +dense0.value[20]*dense1.value[25]
    +dense0.value[21]*dense1.value[24]
    -dense0.value[22]*dense1.value[27]
    +dense0.value[23]*dense1.value[26]
    -dense0.value[24]*dense1.value[21]
    -dense0.value[25]*dense1.value[20]
    +dense0.value[26]*dense1.value[23]
    -dense0.value[27]*dense1.value[22]
    +dense0.value[28]*dense1.value[17]
    -dense0.value[29]*dense1.value[16]
    +dense0.value[30]*dense1.value[19]
    +dense0.value[31]*dense1.value[18]
;
    dense.value[14] =
    +dense0.value[0]*dense1.value[14]
    +dense0.value[1]*dense1.value[15]
    +dense0.value[2]*dense1.value[12]
    -dense0.value[3]*dense1.value[13]
    -dense0.value[4]*dense1.value[10]
    +dense0.value[5]*dense1.value[11]
    +dense0.value[6]*dense1.value[8]
    +dense0.value[7]*dense1.value[9]
    +dense0.value[8]*dense1.value[6]
    -dense0.value[9]*dense1.value[7]
    -dense0.value[10]*dense1.value[4]
    -dense0.value[11]*dense1.value[5]
    +dense0.value[12]*dense1.value[2]
    +dense0.value[13]*dense1.value[3]
    +dense0.value[14]*dense1.value[0]
    -dense0.value[15]*dense1.value[1]
    +dense0.value[16]*dense1.value[30]
    -dense0.value[17]*dense1.value[31]
    -dense0.value[18]*dense1.value[28]
    -dense0.value[19]*dense1.value[29]
    +dense0.value[20]*dense1.value[26]
    +dense0.value[21]*dense1.value[27]
    +dense0.value[22]*dense1.value[24]
    -dense0.value[23]*dense1.value[25]
    -dense0.value[24]*dense1.value[22]
    -dense0.value[25]*dense1.value[23]
    -dense0.value[26]*dense1.value[20]
    +dense0.value[27]*dense1.value[21]
    +dense0.value[28]*dense1.value[18]
    -dense0.value[29]*dense1.value[19]
    -dense0.value[30]*dense1.value[16]
    -dense0.value[31]*dense1.value[17]
;
    dense.value[15] =
    +dense0.value[0]*dense1.value[15]
    +dense0.value[1]*dense1.value[14]
    -dense0.value[2]*dense1.value[13]
    +dense0.value[3]*dense1.value[12]
    +dense0.value[4]*dense1.value[11]
    -dense0.value[5]*dense1.value[10]
    +dense0.value[6]*dense1.value[9]
    +dense0.value[7]*dense1.value[8]
    -dense0.value[8]*dense1.value[7]
    +dense0.value[9]*dense1.value[6]
    -dense0.value[10]*dense1.value[5]
    -dense0.value[11]*dense1.value[4]
    +dense0.value[12]*dense1.value[3]
    +dense0.value[13]*dense1.value[2]
    -dense0.value[14]*dense1.value[1]
    +dense0.value[15]*dense1.value[0]
    -dense0.value[16]*dense1.value[31]
    +dense0.value[17]*dense1.value[30]
    -dense0.value[18]*dense1.value[29]
    -dense0.value[19]*dense1.value[28]
    +dense0.value[20]*dense1.value[27]
    +dense0.value[21]*dense1.value[26]
    -dense0.value[22]*dense1.value[25]
    +dense0.value[23]*dense1.value[24]
    -dense0.value[24]*dense1.value[23]
    -dense0.value[25]*dense1.value[22]
    +dense0.value[26]*dense1.value[21]
    -dense0.value[27]*dense1.value[20]
    -dense0.value[28]*dense1.value[19]
    +dense0.value[29]*dense1.value[18]
    -dense0.value[30]*dense1.value[17]
    -dense0.value[31]*dense1.value[16]
;
    dense.value[16] =
    +dense0.value[0]*dense1.value[16]
    +dense0.value[1]*dense1.value[17]
    +dense0.value[2]*dense1.value[18]
    -dense0.value[3]*dense1.value[19]
    +dense0.value[4]*dense1.value[20]
    -dense0.value[5]*dense1.value[21]
    -dense0.value[6]*dense1.value[22]
    -dense0.value[7]*dense1.value[23]
    +dense0.value[8]*dense1.value[24]
    -dense0.value[9]*dense1.value[25]
    -dense0.value[10]*dense1.value[26]
    -dense0.value[11]*dense1.value[27]
    -dense0.value[12]*dense1.value[28]
    -dense0.value[13]*dense1.value[29]
    -dense0.value[14]*dense1.value[30]
    +dense0.value[15]*dense1.value[31]
    +dense0.value[16]*dense1.value[0]
    -dense0.value[17]*dense1.value[1]
    -dense0.value[18]*dense1.value[2]
    -dense0.value[19]*dense1.value[3]
    -dense0.value[20]*dense1.value[4]
    -dense0.value[21]*dense1.value[5]
    -dense0.value[22]*dense1.value[6]
    +dense0.value[23]*dense1.value[7]
    -dense0.value[24]*dense1.value[8]
    -dense0.value[25]*dense1.value[9]
    -dense0.value[26]*dense1.value[10]
    +dense0.value[27]*dense1.value[11]
    -dense0.value[28]*dense1.value[12]
    +dense0.value[29]*dense1.value[13]
    +dense0.value[30]*dense1.value[14]
    +dense0.value[31]*dense1.value[15]
;
    dense.value[17] =
    +dense0.value[0]*dense1.value[17]
    +dense0.value[1]*dense1.value[16]
    -dense0.value[2]*dense1.value[19]
    +dense0.value[3]*dense1.value[18]
    -dense0.value[4]*dense1.value[21]
    +dense0.value[5]*dense1.value[20]
    -dense0.value[6]*dense1.value[23]
    -dense0.value[7]*dense1.value[22]
    -dense0.value[8]*dense1.value[25]
    +dense0.value[9]*dense1.value[24]
    -dense0.value[10]*dense1.value[27]
    -dense0.value[11]*dense1.value[26]
    -dense0.value[12]*dense1.value[29]
    -dense0.value[13]*dense1.value[28]
    +dense0.value[14]*dense1.value[31]
    -dense0.value[15]*dense1.value[30]
    -dense0.value[16]*dense1.value[1]
    +dense0.value[17]*dense1.value[0]
    -dense0.value[18]*dense1.value[3]
    -dense0.value[19]*dense1.value[2]
    -dense0.value[20]*dense1.value[5]
    -dense0.value[21]*dense1.value[4]
    +dense0.value[22]*dense1.value[7]
    -dense0.value[23]*dense1.value[6]
    -dense0.value[24]*dense1.value[9]
    -dense0.value[25]*dense1.value[8]
    +dense0.value[26]*dense1.value[11]
    -dense0.value[27]*dense1.value[10]
    +dense0.value[28]*dense1.value[13]
    -dense0.value[29]*dense1.value[12]
    +dense0.value[30]*dense1.value[15]
    +dense0.value[31]*dense1.value[14]
;
    dense.value[18] =
    +dense0.value[0]*dense1.value[18]
    +dense0.value[1]*dense1.value[19]
    +dense0.value[2]*dense1.value[16]
    -dense0.value[3]*dense1.value[17]
    -dense0.value[4]*dense1.value[22]
    +dense0.value[5]*dense1.value[23]
    +dense0.value[6]*dense1.value[20]
    +dense0.value[7]*dense1.value[21]
    -dense0.value[8]*dense1.value[26]
    +dense0.value[9]*dense1.value[27]
    +dense0.value[10]*dense1.value[24]
    +dense0.value[11]*dense1.value[25]
    -dense0.value[12]*dense1.value[30]
    -dense0.value[13]*dense1.value[31]
    -dense0.value[14]*dense1.value[28]
    +dense0.value[15]*dense1.value[29]
    -dense0.value[16]*dense1.value[2]
    +dense0.value[17]*dense1.value[3]
    +dense0.value[18]*dense1.value[0]
    +dense0.value[19]*dense1.value[1]
    -dense0.value[20]*dense1.value[6]
    -dense0.value[21]*dense1.value[7]
    -dense0.value[22]*dense1.value[4]
    +dense0.value[23]*dense1.value[5]
    -dense0.value[24]*dense1.value[10]
    -dense0.value[25]*dense1.value[11]
    -dense0.value[26]*dense1.value[8]
    +dense0.value[27]*dense1.value[9]
    +dense0.value[28]*dense1.value[14]
    -dense0.value[29]*dense1.value[15]
    -dense0.value[30]*dense1.value[12]
    -dense0.value[31]*dense1.value[13]
;
    dense.value[19] =
    +dense0.value[0]*dense1.value[19]
    +dense0.value[1]*dense1.value[18]
    -dense0.value[2]*dense1.value[17]
    +dense0.value[3]*dense1.value[16]
    +dense0.value[4]*dense1.value[23]
    -dense0.value[5]*dense1.value[22]
    +dense0.value[6]*dense1.value[21]
    +dense0.value[7]*dense1.value[20]
    +dense0.value[8]*dense1.value[27]
    -dense0.value[9]*dense1.value[26]
    +dense0.value[10]*dense1.value[25]
    +dense0.value[11]*dense1.value[24]
    -dense0.value[12]*dense1.value[31]
    -dense0.value[13]*dense1.value[30]
    +dense0.value[14]*dense1.value[29]
    -dense0.value[15]*dense1.value[28]
    +dense0.value[16]*dense1.value[3]
    -dense0.value[17]*dense1.value[2]
    +dense0.value[18]*dense1.value[1]
    +dense0.value[19]*dense1.value[0]
    -dense0.value[20]*dense1.value[7]
    -dense0.value[21]*dense1.value[6]
    +dense0.value[22]*dense1.value[5]
    -dense0.value[23]*dense1.value[4]
    -dense0.value[24]*dense1.value[11]
    -dense0.value[25]*dense1.value[10]
    +dense0.value[26]*dense1.value[9]
    -dense0.value[27]*dense1.value[8]
    -dense0.value[28]*dense1.value[15]
    +dense0.value[29]*dense1.value[14]
    -dense0.value[30]*dense1.value[13]
    -dense0.value[31]*dense1.value[12]
;
    dense.value[20] =
    +dense0.value[0]*dense1.value[20]
    +dense0.value[1]*dense1.value[21]
    +dense0.value[2]*dense1.value[22]
    -dense0.value[3]*dense1.value[23]
    +dense0.value[4]*dense1.value[16]
    -dense0.value[5]*dense1.value[17]
    -dense0.value[6]*dense1.value[18]
    -dense0.value[7]*dense1.value[19]
    -dense0.value[8]*dense1.value[28]
    +dense0.value[9]*dense1.value[29]
    +dense0.value[10]*dense1.value[30]
    +dense0.value[11]*dense1.value[31]
    +dense0.value[12]*dense1.value[24]
    +dense0.value[13]*dense1.value[25]
    +dense0.value[14]*dense1.value[26]
    -dense0.value[15]*dense1.value[27]
    -dense0.value[16]*dense1.value[4]
    +dense0.value[17]*dense1.value[5]
    +dense0.value[18]*dense1.value[6]
    +dense0.value[19]*dense1.value[7]
    +dense0.value[20]*dense1.value[0]
    +dense0.value[21]*dense1.value[1]
    +dense0.value[22]*dense1.value[2]
    -dense0.value[23]*dense1.value[3]
    -dense0.value[24]*dense1.value[12]
    -dense0.value[25]*dense1.value[13]
    -dense0.value[26]*dense1.value[14]
    +dense0.value[27]*dense1.value[15]
    -dense0.value[28]*dense1.value[8]
    +dense0.value[29]*dense1.value[9]
    +dense0.value[30]*dense1.value[10]
    +dense0.value[31]*dense1.value[11]
;
    dense.value[21] =
    +dense0.value[0]*dense1.value[21]
    +dense0.value[1]*dense1.value[20]
    -dense0.value[2]*dense1.value[23]
    +dense0.value[3]*dense1.value[22]
    -dense0.value[4]*dense1.value[17]
    +dense0.value[5]*dense1.value[16]
    -dense0.value[6]*dense1.value[19]
    -dense0.value[7]*dense1.value[18]
    +dense0.value[8]*dense1.value[29]
    -dense0.value[9]*dense1.value[28]
    +dense0.value[10]*dense1.value[31]
    +dense0.value[11]*dense1.value[30]
    +dense0.value[12]*dense1.value[25]
    +dense0.value[13]*dense1.value[24]
    -dense0.value[14]*dense1.value[27]
    +dense0.value[15]*dense1.value[26]
    +dense0.value[16]*dense1.value[5]
    -dense0.value[17]*dense1.value[4]
    +dense0.value[18]*dense1.value[7]
    +dense0.value[19]*dense1.value[6]
    +dense0.value[20]*dense1.value[1]
    +dense0.value[21]*dense1.value[0]
    -dense0.value[22]*dense1.value[3]
    +dense0.value[23]*dense1.value[2]
    -dense0.value[24]*dense1.value[13]
    -dense0.value[25]*dense1.value[12]
    +dense0.value[26]*dense1.value[15]
    -dense0.value[27]*dense1.value[14]
    +dense0.value[28]*dense1.value[9]
    -dense0.value[29]*dense1.value[8]
    +dense0.value[30]*dense1.value[11]
    +dense0.value[31]*dense1.value[10]
;
    dense.value[22] =
    +dense0.value[0]*dense1.value[22]
    +dense0.value[1]*dense1.value[23]
    +dense0.value[2]*dense1.value[20]
    -dense0.value[3]*dense1.value[21]
    -dense0.value[4]*dense1.value[18]
    +dense0.value[5]*dense1.value[19]
    +dense0.value[6]*dense1.value[16]
    +dense0.value[7]*dense1.value[17]
    +dense0.value[8]*dense1.value[30]
    -dense0.value[9]*dense1.value[31]
    -dense0.value[10]*dense1.value[28]
    -dense0.value[11]*dense1.value[29]
    +dense0.value[12]*dense1.value[26]
    +dense0.value[13]*dense1.value[27]
    +dense0.value[14]*dense1.value[24]
    -dense0.value[15]*dense1.value[25]
    +dense0.value[16]*dense1.value[6]
    -dense0.value[17]*dense1.value[7]
    -dense0.value[18]*dense1.value[4]
    -dense0.value[19]*dense1.value[5]
    +dense0.value[20]*dense1.value[2]
    +dense0.value[21]*dense1.value[3]
    +dense0.value[22]*dense1.value[0]
    -dense0.value[23]*dense1.value[1]
    -dense0.value[24]*dense1.value[14]
    -dense0.value[25]*dense1.value[15]
    -dense0.value[26]*dense1.value[12]
    +dense0.value[27]*dense1.value[13]
    +dense0.value[28]*dense1.value[10]
    -dense0.value[29]*dense1.value[11]
    -dense0.value[30]*dense1.value[8]
    -dense0.value[31]*dense1.value[9]
;
    dense.value[23] =
    +dense0.value[0]*dense1.value[23]
    +dense0.value[1]*dense1.value[22]
    -dense0.value[2]*dense1.value[21]
    +dense0.value[3]*dense1.value[20]
    +dense0.value[4]*dense1.value[19]
    -dense0.value[5]*dense1.value[18]
    +dense0.value[6]*dense1.value[17]
    +dense0.value[7]*dense1.value[16]
    -dense0.value[8]*dense1.value[31]
    +dense0.value[9]*dense1.value[30]
    -dense0.value[10]*dense1.value[29]
    -dense0.value[11]*dense1.value[28]
    +dense0.value[12]*dense1.value[27]
    +dense0.value[13]*dense1.value[26]
    -dense0.value[14]*dense1.value[25]
    +dense0.value[15]*dense1.value[24]
    -dense0.value[16]*dense1.value[7]
    +dense0.value[17]*dense1.value[6]
    -dense0.value[18]*dense1.value[5]
    -dense0.value[19]*dense1.value[4]
    +dense0.value[20]*dense1.value[3]
    +dense0.value[21]*dense1.value[2]
    -dense0.value[22]*dense1.value[1]
    +dense0.value[23]*dense1.value[0]
    -dense0.value[24]*dense1.value[15]
    -dense0.value[25]*dense1.value[14]
    +dense0.value[26]*dense1.value[13]
    -dense0.value[27]*dense1.value[12]
    -dense0.value[28]*dense1.value[11]
    +dense0.value[29]*dense1.value[10]
    -dense0.value[30]*dense1.value[9]
    -dense0.value[31]*dense1.value[8]
;
    dense.value[24] =
    +dense0.value[0]*dense1.value[24]
    +dense0.value[1]*dense1.value[25]
    +dense0.value[2]*dense1.value[26]
    -dense0.value[3]*dense1.value[27]
    +dense0.value[4]*dense1.value[28]
    -dense0.value[5]*dense1.value[29]
    -dense0.value[6]*dense1.value[30]
    -dense0.value[7]*dense1.value[31]
    +dense0.value[8]*dense1.value[16]
    -dense0.value[9]*dense1.value[17]
    -dense0.value[10]*dense1.value[18]
    -dense0.value[11]*dense1.value[19]
    -dense0.value[12]*dense1.value[20]
    -dense0.value[13]*dense1.value[21]
    -dense0.value[14]*dense1.value[22]
    +dense0.value[15]*dense1.value[23]
    -dense0.value[16]*dense1.value[8]
    +dense0.value[17]*dense1.value[9]
    +dense0.value[18]*dense1.value[10]
    +dense0.value[19]*dense1.value[11]
    +dense0.value[20]*dense1.value[12]
    +dense0.value[21]*dense1.value[13]
    +dense0.value[22]*dense1.value[14]
    -dense0.value[23]*dense1.value[15]
    +dense0.value[24]*dense1.value[0]
    +dense0.value[25]*dense1.value[1]
    +dense0.value[26]*dense1.value[2]
    -dense0.value[27]*dense1.value[3]
    +dense0.value[28]*dense1.value[4]
    -dense0.value[29]*dense1.value[5]
    -dense0.value[30]*dense1.value[6]
    -dense0.value[31]*dense1.value[7]
;
    dense.value[25] =
    +dense0.value[0]*dense1.value[25]
    +dense0.value[1]*dense1.value[24]
    -dense0.value[2]*dense1.value[27]
    +dense0.value[3]*dense1.value[26]
    -dense0.value[4]*dense1.value[29]
    +dense0.value[5]*dense1.value[28]
    -dense0.value[6]*dense1.value[31]
    -dense0.value[7]*dense1.value[30]
    -dense0.value[8]*dense1.value[17]
    +dense0.value[9]*dense1.value[16]
    -dense0.value[10]*dense1.value[19]
    -dense0.value[11]*dense1.value[18]
    -dense0.value[12]*dense1.value[21]
    -dense0.value[13]*dense1.value[20]
    +dense0.value[14]*dense1.value[23]
    -dense0.value[15]*dense1.value[22]
    +dense0.value[16]*dense1.value[9]
    -dense0.value[17]*dense1.value[8]
    +dense0.value[18]*dense1.value[11]
    +dense0.value[19]*dense1.value[10]
    +dense0.value[20]*dense1.value[13]
    +dense0.value[21]*dense1.value[12]
    -dense0.value[22]*dense1.value[15]
    +dense0.value[23]*dense1.value[14]
    +dense0.value[24]*dense1.value[1]
    +dense0.value[25]*dense1.value[0]
    -dense0.value[26]*dense1.value[3]
    +dense0.value[27]*dense1.value[2]
    -dense0.value[28]*dense1.value[5]
    +dense0.value[29]*dense1.value[4]
    -dense0.value[30]*dense1.value[7]
    -dense0.value[31]*dense1.value[6]
;
    dense.value[26] =
    +dense0.value[0]*dense1.value[26]
    +dense0.value[1]*dense1.value[27]
    +dense0.value[2]*dense1.value[24]
    -dense0.value[3]*dense1.value[25]
    -dense0.value[4]*dense1.value[30]
    +dense0.value[5]*dense1.value[31]
    +dense0.value[6]*dense1.value[28]
    +dense0.value[7]*dense1.value[29]
    -dense0.value[8]*dense1.value[18]
    +dense0.value[9]*dense1.value[19]
    +dense0.value[10]*dense1.value[16]
    +dense0.value[11]*dense1.value[17]
    -dense0.value[12]*dense1.value[22]
    -dense0.value[13]*dense1.value[23]
    -dense0.value[14]*dense1.value[20]
    +dense0.value[15]*dense1.value[21]
    +dense0.value[16]*dense1.value[10]
    -dense0.value[17]*dense1.value[11]
    -dense0.value[18]*dense1.value[8]
    -dense0.value[19]*dense1.value[9]
    +dense0.value[20]*dense1.value[14]
    +dense0.value[21]*dense1.value[15]
    +dense0.value[22]*dense1.value[12]
    -dense0.value[23]*dense1.value[13]
    +dense0.value[24]*dense1.value[2]
    +dense0.value[25]*dense1.value[3]
    +dense0.value[26]*dense1.value[0]
    -dense0.value[27]*dense1.value[1]
    -dense0.value[28]*dense1.value[6]
    +dense0.value[29]*dense1.value[7]
    +dense0.value[30]*dense1.value[4]
    +dense0.value[31]*dense1.value[5]
;
    dense.value[27] =
    +dense0.value[0]*dense1.value[27]
    +dense0.value[1]*dense1.value[26]
    -dense0.value[2]*dense1.value[25]
    +dense0.value[3]*dense1.value[24]
    +dense0.value[4]*dense1.value[31]
    -dense0.value[5]*dense1.value[30]
    +dense0.value[6]*dense1.value[29]
    +dense0.value[7]*dense1.value[28]
    +dense0.value[8]*dense1.value[19]
    -dense0.value[9]*dense1.value[18]
    +dense0.value[10]*dense1.value[17]
    +dense0.value[11]*dense1.value[16]
    -dense0.value[12]*dense1.value[23]
    -dense0.value[13]*dense1.value[22]
    +dense0.value[14]*dense1.value[21]
    -dense0.value[15]*dense1.value[20]
    -dense0.value[16]*dense1.value[11]
    +dense0.value[17]*dense1.value[10]
    -dense0.value[18]*dense1.value[9]
    -dense0.value[19]*dense1.value[8]
    +dense0.value[20]*dense1.value[15]
    +dense0.value[21]*dense1.value[14]
    -dense0.value[22]*dense1.value[13]
    +dense0.value[23]*dense1.value[12]
    +dense0.value[24]*dense1.value[3]
    +dense0.value[25]*dense1.value[2]
    -dense0.value[26]*dense1.value[1]
    +dense0.value[27]*dense1.value[0]
    +dense0.value[28]*dense1.value[7]
    -dense0.value[29]*dense1.value[6]
    +dense0.value[30]*dense1.value[5]
    +dense0.value[31]*dense1.value[4]
;
    dense.value[28] =
    +dense0.value[0]*dense1.value[28]
    +dense0.value[1]*dense1.value[29]
    +dense0.value[2]*dense1.value[30]
    -dense0.value[3]*dense1.value[31]
    +dense0.value[4]*dense1.value[24]
    -dense0.value[5]*dense1.value[25]
    -dense0.value[6]*dense1.value[26]
    -dense0.value[7]*dense1.value[27]
    -dense0.value[8]*dense1.value[20]
    +dense0.value[9]*dense1.value[21]
    +dense0.value[10]*dense1.value[22]
    +dense0.value[11]*dense1.value[23]
    +dense0.value[12]*dense1.value[16]
    +dense0.value[13]*dense1.value[17]
    +dense0.value[14]*dense1.value[18]
    -dense0.value[15]*dense1.value[19]
    +dense0.value[16]*dense1.value[12]
    -dense0.value[17]*dense1.value[13]
    -dense0.value[18]*dense1.value[14]
    -dense0.value[19]*dense1.value[15]
    -dense0.value[20]*dense1.value[8]
    -dense0.value[21]*dense1.value[9]
    -dense0.value[22]*dense1.value[10]
    +dense0.value[23]*dense1.value[11]
    +dense0.value[24]*dense1.value[4]
    +dense0.value[25]*dense1.value[5]
    +dense0.value[26]*dense1.value[6]
    -dense0.value[27]*dense1.value[7]
    +dense0.value[28]*dense1.value[0]
    -dense0.value[29]*dense1.value[1]
    -dense0.value[30]*dense1.value[2]
    -dense0.value[31]*dense1.value[3]
;
    dense.value[29] =
    +dense0.value[0]*dense1.value[29]
    +dense0.value[1]*dense1.value[28]
    -dense0.value[2]*dense1.value[31]
    +dense0.value[3]*dense1.value[30]
    -dense0.value[4]*dense1.value[25]
    +dense0.value[5]*dense1.value[24]
    -dense0.value[6]*dense1.value[27]
    -dense0.value[7]*dense1.value[26]
    +dense0.value[8]*dense1.value[21]
    -dense0.value[9]*dense1.value[20]
    +dense0.value[10]*dense1.value[23]
    +dense0.value[11]*dense1.value[22]
    +dense0.value[12]*dense1.value[17]
    +dense0.value[13]*dense1.value[16]
    -dense0.value[14]*dense1.value[19]
    +dense0.value[15]*dense1.value[18]
    -dense0.value[16]*dense1.value[13]
    +dense0.value[17]*dense1.value[12]
    -dense0.value[18]*dense1.value[15]
    -dense0.value[19]*dense1.value[14]
    -dense0.value[20]*dense1.value[9]
    -dense0.value[21]*dense1.value[8]
    +dense0.value[22]*dense1.value[11]
    -dense0.value[23]*dense1.value[10]
    +dense0.value[24]*dense1.value[5]
    +dense0.value[25]*dense1.value[4]
    -dense0.value[26]*dense1.value[7]
    +dense0.value[27]*dense1.value[6]
    -dense0.value[28]*dense1.value[1]
    +dense0.value[29]*dense1.value[0]
    -dense0.value[30]*dense1.value[3]
    -dense0.value[31]*dense1.value[2]
;
    dense.value[30] =
    +dense0.value[0]*dense1.value[30]
    +dense0.value[1]*dense1.value[31]
    +dense0.value[2]*dense1.value[28]
    -dense0.value[3]*dense1.value[29]
    -dense0.value[4]*dense1.value[26]
    +dense0.value[5]*dense1.value[27]
    +dense0.value[6]*dense1.value[24]
    +dense0.value[7]*dense1.value[25]
    +dense0.value[8]*dense1.value[22]
    -dense0.value[9]*dense1.value[23]
    -dense0.value[10]*dense1.value[20]
    -dense0.value[11]*dense1.value[21]
    +dense0.value[12]*dense1.value[18]
    +dense0.value[13]*dense1.value[19]
    +dense0.value[14]*dense1.value[16]
    -dense0.value[15]*dense1.value[17]
    -dense0.value[16]*dense1.value[14]
    +dense0.value[17]*dense1.value[15]
    +dense0.value[18]*dense1.value[12]
    +dense0.value[19]*dense1.value[13]
    -dense0.value[20]*dense1.value[10]
    -dense0.value[21]*dense1.value[11]
    -dense0.value[22]*dense1.value[8]
    +dense0.value[23]*dense1.value[9]
    +dense0.value[24]*dense1.value[6]
    +dense0.value[25]*dense1.value[7]
    +dense0.value[26]*dense1.value[4]
    -dense0.value[27]*dense1.value[5]
    -dense0.value[28]*dense1.value[2]
    +dense0.value[29]*dense1.value[3]
    +dense0.value[30]*dense1.value[0]
    +dense0.value[31]*dense1.value[1]
;
    dense.value[31] =
    +dense0.value[0]*dense1.value[31]
    +dense0.value[1]*dense1.value[30]
    -dense0.value[2]*dense1.value[29]
    +dense0.value[3]*dense1.value[28]
    +dense0.value[4]*dense1.value[27]
    -dense0.value[5]*dense1.value[26]
    +dense0.value[6]*dense1.value[25]
    +dense0.value[7]*dense1.value[24]
    -dense0.value[8]*dense1.value[23]
    +dense0.value[9]*dense1.value[22]
    -dense0.value[10]*dense1.value[21]
    -dense0.value[11]*dense1.value[20]
    +dense0.value[12]*dense1.value[19]
    +dense0.value[13]*dense1.value[18]
    -dense0.value[14]*dense1.value[17]
    +dense0.value[15]*dense1.value[16]
    +dense0.value[16]*dense1.value[15]
    -dense0.value[17]*dense1.value[14]
    +dense0.value[18]*dense1.value[13]
    +dense0.value[19]*dense1.value[12]
    -dense0.value[20]*dense1.value[11]
    -dense0.value[21]*dense1.value[10]
    +dense0.value[22]*dense1.value[9]
    -dense0.value[23]*dense1.value[8]
    +dense0.value[24]*dense1.value[7]
    +dense0.value[25]*dense1.value[6]
    -dense0.value[26]*dense1.value[5]
    +dense0.value[27]*dense1.value[4]
    +dense0.value[28]*dense1.value[3]
    -dense0.value[29]*dense1.value[2]
    +dense0.value[30]*dense1.value[1]
    +dense0.value[31]*dense1.value[0]
;
    return dense;
}


static gen1_DenseMultivector gen1_dense_innerproduct(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1){
    gen1_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[1]*dense1.value[1]
    +dense0.value[2]*dense1.value[2]
    -dense0.value[3]*dense1.value[3]
    +dense0.value[4]*dense1.value[4]
    -dense0.value[5]*dense1.value[5]
    -dense0.value[6]*dense1.value[6]
    -dense0.value[7]*dense1.value[7]
    +dense0.value[8]*dense1.value[8]
    -dense0.value[9]*dense1.value[9]
    -dense0.value[10]*dense1.value[10]
    -dense0.value[11]*dense1.value[11]
    -dense0.value[12]*dense1.value[12]
    -dense0.value[13]*dense1.value[13]
    -dense0.value[14]*dense1.value[14]
    +dense0.value[15]*dense1.value[15]
    -dense0.value[16]*dense1.value[16]
    +dense0.value[17]*dense1.value[17]
    +dense0.value[18]*dense1.value[18]
    +dense0.value[19]*dense1.value[19]
    +dense0.value[20]*dense1.value[20]
    +dense0.value[21]*dense1.value[21]
    +dense0.value[22]*dense1.value[22]
    -dense0.value[23]*dense1.value[23]
    +dense0.value[24]*dense1.value[24]
    +dense0.value[25]*dense1.value[25]
    +dense0.value[26]*dense1.value[26]
    -dense0.value[27]*dense1.value[27]
    +dense0.value[28]*dense1.value[28]
    -dense0.value[29]*dense1.value[29]
    -dense0.value[30]*dense1.value[30]
    -dense0.value[31]*dense1.value[31]
;
    dense.value[1] =
    -dense0.value[2]*dense1.value[3]
    +dense0.value[3]*dense1.value[2]
    -dense0.value[4]*dense1.value[5]
    +dense0.value[5]*dense1.value[4]
    -dense0.value[6]*dense1.value[7]
    -dense0.value[7]*dense1.value[6]
    -dense0.value[8]*dense1.value[9]
    +dense0.value[9]*dense1.value[8]
    -dense0.value[10]*dense1.value[11]
    -dense0.value[11]*dense1.value[10]
    -dense0.value[12]*dense1.value[13]
    -dense0.value[13]*dense1.value[12]
    +dense0.value[14]*dense1.value[15]
    -dense0.value[15]*dense1.value[14]
    +dense0.value[16]*dense1.value[17]
    -dense0.value[17]*dense1.value[16]
    +dense0.value[18]*dense1.value[19]
    +dense0.value[19]*dense1.value[18]
    +dense0.value[20]*dense1.value[21]
    +dense0.value[21]*dense1.value[20]
    -dense0.value[22]*dense1.value[23]
    +dense0.value[23]*dense1.value[22]
    +dense0.value[24]*dense1.value[25]
    +dense0.value[25]*dense1.value[24]
    -dense0.value[26]*dense1.value[27]
    +dense0.value[27]*dense1.value[26]
    -dense0.value[28]*dense1.value[29]
    +dense0.value[29]*dense1.value[28]
    -dense0.value[30]*dense1.value[31]
    -dense0.value[31]*dense1.value[30]
;
    dense.value[2] =
    +dense0.value[1]*dense1.value[3]
    -dense0.value[3]*dense1.value[1]
    -dense0.value[4]*dense1.value[6]
    +dense0.value[5]*dense1.value[7]
    +dense0.value[6]*dense1.value[4]
    +dense0.value[7]*dense1.value[5]
    -dense0.value[8]*dense1.value[10]
    +dense0.value[9]*dense1.value[11]
    +dense0.value[10]*dense1.value[8]
    +dense0.value[11]*dense1.value[9]
    -dense0.value[12]*dense1.value[14]
    -dense0.value[13]*dense1.value[15]
    -dense0.value[14]*dense1.value[12]
    +dense0.value[15]*dense1.value[13]
    +dense0.value[16]*dense1.value[18]
    -dense0.value[17]*dense1.value[19]
    -dense0.value[18]*dense1.value[16]
    -dense0.value[19]*dense1.value[17]
    +dense0.value[20]*dense1.value[22]
    +dense0.value[21]*dense1.value[23]
    +dense0.value[22]*dense1.value[20]
    -dense0.value[23]*dense1.value[21]
    +dense0.value[24]*dense1.value[26]
    +dense0.value[25]*dense1.value[27]
    +dense0.value[26]*dense1.value[24]
    -dense0.value[27]*dense1.value[25]
    -dense0.value[28]*dense1.value[30]
    +dense0.value[29]*dense1.value[31]
    +dense0.value[30]*dense1.value[28]
    +dense0.value[31]*dense1.value[29]
;
    dense.value[3] =
    +dense0.value[4]*dense1.value[7]
    +dense0.value[7]*dense1.value[4]
    +dense0.value[8]*dense1.value[11]
    +dense0.value[11]*dense1.value[8]
    -dense0.value[12]*dense1.value[15]
    -dense0.value[15]*dense1.value[12]
    -dense0.value[16]*dense1.value[19]
    -dense0.value[19]*dense1.value[16]
    +dense0.value[20]*dense1.value[23]
    +dense0.value[23]*dense1.value[20]
    +dense0.value[24]*dense1.value[27]
    +dense0.value[27]*dense1.value[24]
    +dense0.value[28]*dense1.value[31]
    +dense0.value[31]*dense1.value[28]
;
    dense.value[4] =
    +dense0.value[1]*dense1.value[5]
    +dense0.value[2]*dense1.value[6]
    -dense0.value[3]*dense1.value[7]
    -dense0.value[5]*dense1.value[1]
    -dense0.value[6]*dense1.value[2]
    -dense0.value[7]*dense1.value[3]
    -dense0.value[8]*dense1.value[12]
    +dense0.value[9]*dense1.value[13]
    +dense0.value[10]*dense1.value[14]
    +dense0.value[11]*dense1.value[15]
    +dense0.value[12]*dense1.value[8]
    +dense0.value[13]*dense1.value[9]
    +dense0.value[14]*dense1.value[10]
    -dense0.value[15]*dense1.value[11]
    +dense0.value[16]*dense1.value[20]
    -dense0.value[17]*dense1.value[21]
    -dense0.value[18]*dense1.value[22]
    -dense0.value[19]*dense1.value[23]
    -dense0.value[20]*dense1.value[16]
    -dense0.value[21]*dense1.value[17]
    -dense0.value[22]*dense1.value[18]
    +dense0.value[23]*dense1.value[19]
    +dense0.value[24]*dense1.value[28]
    +dense0.value[25]*dense1.value[29]
    +dense0.value[26]*dense1.value[30]
    -dense0.value[27]*dense1.value[31]
    +dense0.value[28]*dense1.value[24]
    -dense0.value[29]*dense1.value[25]
    -dense0.value[30]*dense1.value[26]
    -dense0.value[31]*dense1.value[27]
;
    dense.value[5] =
    -dense0.value[2]*dense1.value[7]
    -dense0.value[7]*dense1.value[2]
    +dense0.value[8]*dense1.value[13]
    +dense0.value[10]*dense1.value[15]
    +dense0.value[13]*dense1.value[8]
    +dense0.value[15]*dense1.value[10]
    -dense0.value[16]*dense1.value[21]
    -dense0.value[18]*dense1.value[23]
    -dense0.value[21]*dense1.value[16]
    -dense0.value[23]*dense1.value[18]
    +dense0.value[24]*dense1.value[29]
    -dense0.value[26]*dense1.value[31]
    +dense0.value[29]*dense1.value[24]
    -dense0.value[31]*dense1.value[26]
;
    dense.value[6] =
    +dense0.value[1]*dense1.value[7]
    +dense0.value[7]*dense1.value[1]
    +dense0.value[8]*dense1.value[14]
    -dense0.value[9]*dense1.value[15]
    +dense0.value[14]*dense1.value[8]
    -dense0.value[15]*dense1.value[9]
    -dense0.value[16]*dense1.value[22]
    +dense0.value[17]*dense1.value[23]
    -dense0.value[22]*dense1.value[16]
    +dense0.value[23]*dense1.value[17]
    +dense0.value[24]*dense1.value[30]
    +dense0.value[25]*dense1.value[31]
    +dense0.value[30]*dense1.value[24]
    +dense0.value[31]*dense1.value[25]
;
    dense.value[7] =
    -dense0.value[8]*dense1.value[15]
    +dense0.value[15]*dense1.value[8]
    +dense0.value[16]*dense1.value[23]
    -dense0.value[23]*dense1.value[16]
    +dense0.value[24]*dense1.value[31]
    +dense0.value[31]*dense1.value[24]
;
    dense.value[8] =
    +dense0.value[1]*dense1.value[9]
    +dense0.value[2]*dense1.value[10]
    -dense0.value[3]*dense1.value[11]
    +dense0.value[4]*dense1.value[12]
    -dense0.value[5]*dense1.value[13]
    -dense0.value[6]*dense1.value[14]
    -dense0.value[7]*dense1.value[15]
    -dense0.value[9]*dense1.value[1]
    -dense0.value[10]*dense1.value[2]
    -dense0.value[11]*dense1.value[3]
    -dense0.value[12]*dense1.value[4]
    -dense0.value[13]*dense1.value[5]
    -dense0.value[14]*dense1.value[6]
    +dense0.value[15]*dense1.value[7]
    +dense0.value[16]*dense1.value[24]
    -dense0.value[17]*dense1.value[25]
    -dense0.value[18]*dense1.value[26]
    -dense0.value[19]*dense1.value[27]
    -dense0.value[20]*dense1.value[28]
    -dense0.value[21]*dense1.value[29]
    -dense0.value[22]*dense1.value[30]
    +dense0.value[23]*dense1.value[31]
    -dense0.value[24]*dense1.value[16]
    -dense0.value[25]*dense1.value[17]
    -dense0.value[26]*dense1.value[18]
    +dense0.value[27]*dense1.value[19]
    -dense0.value[28]*dense1.value[20]
    +dense0.value[29]*dense1.value[21]
    +dense0.value[30]*dense1.value[22]
    +dense0.value[31]*dense1.value[23]
;
    dense.value[9] =
    -dense0.value[2]*dense1.value[11]
    -dense0.value[4]*dense1.value[13]
    -dense0.value[6]*dense1.value[15]
    -dense0.value[11]*dense1.value[2]
    -dense0.value[13]*dense1.value[4]
    -dense0.value[15]*dense1.value[6]
    -dense0.value[16]*dense1.value[25]
    -dense0.value[18]*dense1.value[27]
    -dense0.value[20]*dense1.value[29]
    +dense0.value[22]*dense1.value[31]
    -dense0.value[25]*dense1.value[16]
    -dense0.value[27]*dense1.value[18]
    -dense0.value[29]*dense1.value[20]
    +dense0.value[31]*dense1.value[22]
;
    dense.value[10] =
    +dense0.value[1]*dense1.value[11]
    -dense0.value[4]*dense1.value[14]
    +dense0.value[5]*dense1.value[15]
    +dense0.value[11]*dense1.value[1]
    -dense0.value[14]*dense1.value[4]
    +dense0.value[15]*dense1.value[5]
    -dense0.value[16]*dense1.value[26]
    +dense0.value[17]*dense1.value[27]
    -dense0.value[20]*dense1.value[30]
    -dense0.value[21]*dense1.value[31]
    -dense0.value[26]*dense1.value[16]
    +dense0.value[27]*dense1.value[17]
    -dense0.value[30]*dense1.value[20]
    -dense0.value[31]*dense1.value[21]
;
    dense.value[11] =
    +dense0.value[4]*dense1.value[15]
    -dense0.value[15]*dense1.value[4]
    +dense0.value[16]*dense1.value[27]
    -dense0.value[20]*dense1.value[31]
    -dense0.value[27]*dense1.value[16]
    -dense0.value[31]*dense1.value[20]
;
    dense.value[12] =
    +dense0.value[1]*dense1.value[13]
    +dense0.value[2]*dense1.value[14]
    -dense0.value[3]*dense1.value[15]
    +dense0.value[13]*dense1.value[1]
    +dense0.value[14]*dense1.value[2]
    -dense0.value[15]*dense1.value[3]
    -dense0.value[16]*dense1.value[28]
    +dense0.value[17]*dense1.value[29]
    +dense0.value[18]*dense1.value[30]
    +dense0.value[19]*dense1.value[31]
    -dense0.value[28]*dense1.value[16]
    +dense0.value[29]*dense1.value[17]
    +dense0.value[30]*dense1.value[18]
    +dense0.value[31]*dense1.value[19]
;
    dense.value[13] =
    -dense0.value[2]*dense1.value[15]
    +dense0.value[15]*dense1.value[2]
    +dense0.value[16]*dense1.value[29]
    +dense0.value[18]*dense1.value[31]
    -dense0.value[29]*dense1.value[16]
    +dense0.value[31]*dense1.value[18]
;
    dense.value[14] =
    +dense0.value[1]*dense1.value[15]
    -dense0.value[15]*dense1.value[1]
    +dense0.value[16]*dense1.value[30]
    -dense0.value[17]*dense1.value[31]
    -dense0.value[30]*dense1.value[16]
    -dense0.value[31]*dense1.value[17]
;
    dense.value[15] =
    -dense0.value[16]*dense1.value[31]
    -dense0.value[31]*dense1.value[16]
;
    dense.value[16] =
    +dense0.value[1]*dense1.value[17]
    +dense0.value[2]*dense1.value[18]
    -dense0.value[3]*dense1.value[19]
    +dense0.value[4]*dense1.value[20]
    -dense0.value[5]*dense1.value[21]
    -dense0.value[6]*dense1.value[22]
    -dense0.value[7]*dense1.value[23]
    +dense0.value[8]*dense1.value[24]
    -dense0.value[9]*dense1.value[25]
    -dense0.value[10]*dense1.value[26]
    -dense0.value[11]*dense1.value[27]
    -dense0.value[12]*dense1.value[28]
    -dense0.value[13]*dense1.value[29]
    -dense0.value[14]*dense1.value[30]
    +dense0.value[15]*dense1.value[31]
    -dense0.value[17]*dense1.value[1]
    -dense0.value[18]*dense1.value[2]
    -dense0.value[19]*dense1.value[3]
    -dense0.value[20]*dense1.value[4]
    -dense0.value[21]*dense1.value[5]
    -dense0.value[22]*dense1.value[6]
    +dense0.value[23]*dense1.value[7]
    -dense0.value[24]*dense1.value[8]
    -dense0.value[25]*dense1.value[9]
    -dense0.value[26]*dense1.value[10]
    +dense0.value[27]*dense1.value[11]
    -dense0.value[28]*dense1.value[12]
    +dense0.value[29]*dense1.value[13]
    +dense0.value[30]*dense1.value[14]
    +dense0.value[31]*dense1.value[15]
;
    dense.value[17] =
    -dense0.value[2]*dense1.value[19]
    -dense0.value[4]*dense1.value[21]
    -dense0.value[6]*dense1.value[23]
    -dense0.value[8]*dense1.value[25]
    -dense0.value[10]*dense1.value[27]
    -dense0.value[12]*dense1.value[29]
    +dense0.value[14]*dense1.value[31]
    -dense0.value[19]*dense1.value[2]
    -dense0.value[21]*dense1.value[4]
    -dense0.value[23]*dense1.value[6]
    -dense0.value[25]*dense1.value[8]
    -dense0.value[27]*dense1.value[10]
    -dense0.value[29]*dense1.value[12]
    +dense0.value[31]*dense1.value[14]
;
    dense.value[18] =
    +dense0.value[1]*dense1.value[19]
    -dense0.value[4]*dense1.value[22]
    +dense0.value[5]*dense1.value[23]
    -dense0.value[8]*dense1.value[26]
    +dense0.value[9]*dense1.value[27]
    -dense0.value[12]*dense1.value[30]
    -dense0.value[13]*dense1.value[31]
    +dense0.value[19]*dense1.value[1]
    -dense0.value[22]*dense1.value[4]
    +dense0.value[23]*dense1.value[5]
    -dense0.value[26]*dense1.value[8]
    +dense0.value[27]*dense1.value[9]
    -dense0.value[30]*dense1.value[12]
    -dense0.value[31]*dense1.value[13]
;
    dense.value[19] =
    +dense0.value[4]*dense1.value[23]
    +dense0.value[8]*dense1.value[27]
    -dense0.value[12]*dense1.value[31]
    -dense0.value[23]*dense1.value[4]
    -dense0.value[27]*dense1.value[8]
    -dense0.value[31]*dense1.value[12]
;
    dense.value[20] =
    +dense0.value[1]*dense1.value[21]
    +dense0.value[2]*dense1.value[22]
    -dense0.value[3]*dense1.value[23]
    -dense0.value[8]*dense1.value[28]
    +dense0.value[9]*dense1.value[29]
    +dense0.value[10]*dense1.value[30]
    +dense0.value[11]*dense1.value[31]
    +dense0.value[21]*dense1.value[1]
    +dense0.value[22]*dense1.value[2]
    -dense0.value[23]*dense1.value[3]
    -dense0.value[28]*dense1.value[8]
    +dense0.value[29]*dense1.value[9]
    +dense0.value[30]*dense1.value[10]
    +dense0.value[31]*dense1.value[11]
;
    dense.value[21] =
    -dense0.value[2]*dense1.value[23]
    +dense0.value[8]*dense1.value[29]
    +dense0.value[10]*dense1.value[31]
    +dense0.value[23]*dense1.value[2]
    -dense0.value[29]*dense1.value[8]
    +dense0.value[31]*dense1.value[10]
;
    dense.value[22] =
    +dense0.value[1]*dense1.value[23]
    +dense0.value[8]*dense1.value[30]
    -dense0.value[9]*dense1.value[31]
    -dense0.value[23]*dense1.value[1]
    -dense0.value[30]*dense1.value[8]
    -dense0.value[31]*dense1.value[9]
;
    dense.value[23] =
    -dense0.value[8]*dense1.value[31]
    -dense0.value[31]*dense1.value[8]
;
    dense.value[24] =
    +dense0.value[1]*dense1.value[25]
    +dense0.value[2]*dense1.value[26]
    -dense0.value[3]*dense1.value[27]
    +dense0.value[4]*dense1.value[28]
    -dense0.value[5]*dense1.value[29]
    -dense0.value[6]*dense1.value[30]
    -dense0.value[7]*dense1.value[31]
    +dense0.value[25]*dense1.value[1]
    +dense0.value[26]*dense1.value[2]
    -dense0.value[27]*dense1.value[3]
    +dense0.value[28]*dense1.value[4]
    -dense0.value[29]*dense1.value[5]
    -dense0.value[30]*dense1.value[6]
    -dense0.value[31]*dense1.value[7]
;
    dense.value[25] =
    -dense0.value[2]*dense1.value[27]
    -dense0.value[4]*dense1.value[29]
    -dense0.value[6]*dense1.value[31]
    +dense0.value[27]*dense1.value[2]
    +dense0.value[29]*dense1.value[4]
    -dense0.value[31]*dense1.value[6]
;
    dense.value[26] =
    +dense0.value[1]*dense1.value[27]
    -dense0.value[4]*dense1.value[30]
    +dense0.value[5]*dense1.value[31]
    -dense0.value[27]*dense1.value[1]
    +dense0.value[30]*dense1.value[4]
    +dense0.value[31]*dense1.value[5]
;
    dense.value[27] =
    +dense0.value[4]*dense1.value[31]
    +dense0.value[31]*dense1.value[4]
;
    dense.value[28] =
    +dense0.value[1]*dense1.value[29]
    +dense0.value[2]*dense1.value[30]
    -dense0.value[3]*dense1.value[31]
    -dense0.value[29]*dense1.value[1]
    -dense0.value[30]*dense1.value[2]
    -dense0.value[31]*dense1.value[3]
;
    dense.value[29] =
    -dense0.value[2]*dense1.value[31]
    -dense0.value[31]*dense1.value[2]
;
    dense.value[30] =
    +dense0.value[1]*dense1.value[31]
    +dense0.value[31]*dense1.value[1]
;
    return dense;
}


static gen1_DenseMultivector gen1_dense_outerproduct(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1){
    gen1_DenseMultivector dense = {{0}};
    dense.value[0] =
    +dense0.value[0]*dense1.value[0]
;
    dense.value[1] =
    +dense0.value[0]*dense1.value[1]
    +dense0.value[1]*dense1.value[0]
;
    dense.value[2] =
    +dense0.value[0]*dense1.value[2]
    +dense0.value[2]*dense1.value[0]
;
    dense.value[3] =
    +dense0.value[0]*dense1.value[3]
    +dense0.value[1]*dense1.value[2]
    -dense0.value[2]*dense1.value[1]
    +dense0.value[3]*dense1.value[0]
;
    dense.value[4] =
    +dense0.value[0]*dense1.value[4]
    +dense0.value[4]*dense1.value[0]
;
    dense.value[5] =
    +dense0.value[0]*dense1.value[5]
    +dense0.value[1]*dense1.value[4]
    -dense0.value[4]*dense1.value[1]
    +dense0.value[5]*dense1.value[0]
;
    dense.value[6] =
    +dense0.value[0]*dense1.value[6]
    +dense0.value[2]*dense1.value[4]
    -dense0.value[4]*dense1.value[2]
    +dense0.value[6]*dense1.value[0]
;
    dense.value[7] =
    +dense0.value[0]*dense1.value[7]
    +dense0.value[1]*dense1.value[6]
    -dense0.value[2]*dense1.value[5]
    +dense0.value[3]*dense1.value[4]
    +dense0.value[4]*dense1.value[3]
    -dense0.value[5]*dense1.value[2]
    +dense0.value[6]*dense1.value[1]
    +dense0.value[7]*dense1.value[0]
;
    dense.value[8] =
    +dense0.value[0]*dense1.value[8]
    +dense0.value[8]*dense1.value[0]
;
    dense.value[9] =
    +dense0.value[0]*dense1.value[9]
    +dense0.value[1]*dense1.value[8]
    -dense0.value[8]*dense1.value[1]
    +dense0.value[9]*dense1.value[0]
;
    dense.value[10] =
    +dense0.value[0]*dense1.value[10]
    +dense0.value[2]*dense1.value[8]
    -dense0.value[8]*dense1.value[2]
    +dense0.value[10]*dense1.value[0]
;
    dense.value[11] =
    +dense0.value[0]*dense1.value[11]
    +dense0.value[1]*dense1.value[10]
    -dense0.value[2]*dense1.value[9]
    +dense0.value[3]*dense1.value[8]
    +dense0.value[8]*dense1.value[3]
    -dense0.value[9]*dense1.value[2]
    +dense0.value[10]*dense1.value[1]
    +dense0.value[11]*dense1.value[0]
;
    dense.value[12] =
    +dense0.value[0]*dense1.value[12]
    +dense0.value[4]*dense1.value[8]
    -dense0.value[8]*dense1.value[4]
    +dense0.value[12]*dense1.value[0]
;
    dense.value[13] =
    +dense0.value[0]*dense1.value[13]
    +dense0.value[1]*dense1.value[12]
    -dense0.value[4]*dense1.value[9]
    +dense0.value[5]*dense1.value[8]
    +dense0.value[8]*dense1.value[5]
    -dense0.value[9]*dense1.value[4]
    +dense0.value[12]*dense1.value[1]
    +dense0.value[13]*dense1.value[0]
;
    dense.value[14] =
    +dense0.value[0]*dense1.value[14]
    +dense0.value[2]*dense1.value[12]
    -dense0.value[4]*dense1.value[10]
    +dense0.value[6]*dense1.value[8]
    +dense0.value[8]*dense1.value[6]
    -dense0.value[10]*dense1.value[4]
    +dense0.value[12]*dense1.value[2]
    +dense0.value[14]*dense1.value[0]
;
    dense.value[15] =
    +dense0.value[0]*dense1.value[15]
    +dense0.value[1]*dense1.value[14]
    -dense0.value[2]*dense1.value[13]
    +dense0.value[3]*dense1.value[12]
    +dense0.value[4]*dense1.value[11]
    -dense0.value[5]*dense1.value[10]
    +dense0.value[6]*dense1.value[9]
    +dense0.value[7]*dense1.value[8]
    -dense0.value[8]*dense1.value[7]
    +dense0.value[9]*dense1.value[6]
    -dense0.value[10]*dense1.value[5]
    -dense0.value[11]*dense1.value[4]
    +dense0.value[12]*dense1.value[3]
    +dense0.value[13]*dense1.value[2]
    -dense0.value[14]*dense1.value[1]
    +dense0.value[15]*dense1.value[0]
;
    dense.value[16] =
    +dense0.value[0]*dense1.value[16]
    +dense0.value[16]*dense1.value[0]
;
    dense.value[17] =
    +dense0.value[0]*dense1.value[17]
    +dense0.value[1]*dense1.value[16]
    -dense0.value[16]*dense1.value[1]
    +dense0.value[17]*dense1.value[0]
;
    dense.value[18] =
    +dense0.value[0]*dense1.value[18]
    +dense0.value[2]*dense1.value[16]
    -dense0.value[16]*dense1.value[2]
    +dense0.value[18]*dense1.value[0]
;
    dense.value[19] =
    +dense0.value[0]*dense1.value[19]
    +dense0.value[1]*dense1.value[18]
    -dense0.value[2]*dense1.value[17]
    +dense0.value[3]*dense1.value[16]
    +dense0.value[16]*dense1.value[3]
    -dense0.value[17]*dense1.value[2]
    +dense0.value[18]*dense1.value[1]
    +dense0.value[19]*dense1.value[0]
;
    dense.value[20] =
    +dense0.value[0]*dense1.value[20]
    +dense0.value[4]*dense1.value[16]
    -dense0.value[16]*dense1.value[4]
    +dense0.value[20]*dense1.value[0]
;
    dense.value[21] =
    +dense0.value[0]*dense1.value[21]
    +dense0.value[1]*dense1.value[20]
    -dense0.value[4]*dense1.value[17]
    +dense0.value[5]*dense1.value[16]
    +dense0.value[16]*dense1.value[5]
    -dense0.value[17]*dense1.value[4]
    +dense0.value[20]*dense1.value[1]
    +dense0.value[21]*dense1.value[0]
;
    dense.value[22] =
    +dense0.value[0]*dense1.value[22]
    +dense0.value[2]*dense1.value[20]
    -dense0.value[4]*dense1.value[18]
    +dense0.value[6]*dense1.value[16]
    +dense0.value[16]*dense1.value[6]
    -dense0.value[18]*dense1.value[4]
    +dense0.value[20]*dense1.value[2]
    +dense0.value[22]*dense1.value[0]
;
    dense.value[23] =
    +dense0.value[0]*dense1.value[23]
    +dense0.value[1]*dense1.value[22]
    -dense0.value[2]*dense1.value[21]
    +dense0.value[3]*dense1.value[20]
    +dense0.value[4]*dense1.value[19]
    -dense0.value[5]*dense1.value[18]
    +dense0.value[6]*dense1.value[17]
    +dense0.value[7]*dense1.value[16]
    -dense0.value[16]*dense1.value[7]
    +dense0.value[17]*dense1.value[6]
    -dense0.value[18]*dense1.value[5]
    -dense0.value[19]*dense1.value[4]
    +dense0.value[20]*dense1.value[3]
    +dense0.value[21]*dense1.value[2]
    -dense0.value[22]*dense1.value[1]
    +dense0.value[23]*dense1.value[0]
;
    dense.value[24] =
    +dense0.value[0]*dense1.value[24]
    +dense0.value[8]*dense1.value[16]
    -dense0.value[16]*dense1.value[8]
    +dense0.value[24]*dense1.value[0]
;
    dense.value[25] =
    +dense0.value[0]*dense1.value[25]
    +dense0.value[1]*dense1.value[24]
    -dense0.value[8]*dense1.value[17]
    +dense0.value[9]*dense1.value[16]
    +dense0.value[16]*dense1.value[9]
    -dense0.value[17]*dense1.value[8]
    +dense0.value[24]*dense1.value[1]
    +dense0.value[25]*dense1.value[0]
;
    dense.value[26] =
    +dense0.value[0]*dense1.value[26]
    +dense0.value[2]*dense1.value[24]
    -dense0.value[8]*dense1.value[18]
    +dense0.value[10]*dense1.value[16]
    +dense0.value[16]*dense1.value[10]
    -dense0.value[18]*dense1.value[8]
    +dense0.value[24]*dense1.value[2]
    +dense0.value[26]*dense1.value[0]
;
    dense.value[27] =
    +dense0.value[0]*dense1.value[27]
    +dense0.value[1]*dense1.value[26]
    -dense0.value[2]*dense1.value[25]
    +dense0.value[3]*dense1.value[24]
    +dense0.value[8]*dense1.value[19]
    -dense0.value[9]*dense1.value[18]
    +dense0.value[10]*dense1.value[17]
    +dense0.value[11]*dense1.value[16]
    -dense0.value[16]*dense1.value[11]
    +dense0.value[17]*dense1.value[10]
    -dense0.value[18]*dense1.value[9]
    -dense0.value[19]*dense1.value[8]
    +dense0.value[24]*dense1.value[3]
    +dense0.value[25]*dense1.value[2]
    -dense0.value[26]*dense1.value[1]
    +dense0.value[27]*dense1.value[0]
;
    dense.value[28] =
    +dense0.value[0]*dense1.value[28]
    +dense0.value[4]*dense1.value[24]
    -dense0.value[8]*dense1.value[20]
    +dense0.value[12]*dense1.value[16]
    +dense0.value[16]*dense1.value[12]
    -dense0.value[20]*dense1.value[8]
    +dense0.value[24]*dense1.value[4]
    +dense0.value[28]*dense1.value[0]
;
    dense.value[29] =
    +dense0.value[0]*dense1.value[29]
    +dense0.value[1]*dense1.value[28]
    -dense0.value[4]*dense1.value[25]
    +dense0.value[5]*dense1.value[24]
    +dense0.value[8]*dense1.value[21]
    -dense0.value[9]*dense1.value[20]
    +dense0.value[12]*dense1.value[17]
    +dense0.value[13]*dense1.value[16]
    -dense0.value[16]*dense1.value[13]
    +dense0.value[17]*dense1.value[12]
    -dense0.value[20]*dense1.value[9]
    -dense0.value[21]*dense1.value[8]
    +dense0.value[24]*dense1.value[5]
    +dense0.value[25]*dense1.value[4]
    -dense0.value[28]*dense1.value[1]
    +dense0.value[29]*dense1.value[0]
;
    dense.value[30] =
    +dense0.value[0]*dense1.value[30]
    +dense0.value[2]*dense1.value[28]
    -dense0.value[4]*dense1.value[26]
    +dense0.value[6]*dense1.value[24]
    +dense0.value[8]*dense1.value[22]
    -dense0.value[10]*dense1.value[20]
    +dense0.value[12]*dense1.value[18]
    +dense0.value[14]*dense1.value[16]
    -dense0.value[16]*dense1.value[14]
    +dense0.value[18]*dense1.value[12]
    -dense0.value[20]*dense1.value[10]
    -dense0.value[22]*dense1.value[8]
    +dense0.value[24]*dense1.value[6]
    +dense0.value[26]*dense1.value[4]
    -dense0.value[28]*dense1.value[2]
    +dense0.value[30]*dense1.value[0]
;
    dense.value[31] =
    +dense0.value[0]*dense1.value[31]
    +dense0.value[1]*dense1.value[30]
    -dense0.value[2]*dense1.value[29]
    +dense0.value[3]*dense1.value[28]
    +dense0.value[4]*dense1.value[27]
    -dense0.value[5]*dense1.value[26]
    +dense0.value[6]*dense1.value[25]
    +dense0.value[7]*dense1.value[24]
    -dense0.value[8]*dense1.value[23]
    +dense0.value[9]*dense1.value[22]
    -dense0.value[10]*dense1.value[21]
    -dense0.value[11]*dense1.value[20]
    +dense0.value[12]*dense1.value[19]
    +dense0.value[13]*dense1.value[18]
    -dense0.value[14]*dense1.value[17]
    +dense0.value[15]*dense1.value[16]
    +dense0.value[16]*dense1.value[15]
    -dense0.value[17]*dense1.value[14]
    +dense0.value[18]*dense1.value[13]
    +dense0.value[19]*dense1.value[12]
    -dense0.value[20]*dense1.value[11]
    -dense0.value[21]*dense1.value[10]
    +dense0.value[22]*dense1.value[9]
    -dense0.value[23]*dense1.value[8]
    +dense0.value[24]*dense1.value[7]
    +dense0.value[25]*dense1.value[6]
    -dense0.value[26]*dense1.value[5]
    +dense0.value[27]*dense1.value[4]
    +dense0.value[28]*dense1.value[3]
    -dense0.value[29]*dense1.value[2]
    +dense0.value[30]*dense1.value[1]
    +dense0.value[31]*dense1.value[0]
;
    return dense;
}


static gen1_DenseMultivector gen1_dense_atomicadd(gen1_DenseMultivector *dense_array, Py_ssize_t size){
    gen1_DenseMultivector dense = {{0}};

    for(Py_ssize_t i = 0; i < size; i++){
        dense.value[0] += dense_array[i].value[0];
        dense.value[1] += dense_array[i].value[1];
        dense.value[2] += dense_array[i].value[2];
        dense.value[3] += dense_array[i].value[3];
        dense.value[4] += dense_array[i].value[4];
        dense.value[5] += dense_array[i].value[5];
        dense.value[6] += dense_array[i].value[6];
        dense.value[7] += dense_array[i].value[7];
        dense.value[8] += dense_array[i].value[8];
        dense.value[9] += dense_array[i].value[9];
        dense.value[10] += dense_array[i].value[10];
        dense.value[11] += dense_array[i].value[11];
        dense.value[12] += dense_array[i].value[12];
        dense.value[13] += dense_array[i].value[13];
        dense.value[14] += dense_array[i].value[14];
        dense.value[15] += dense_array[i].value[15];
        dense.value[16] += dense_array[i].value[16];
        dense.value[17] += dense_array[i].value[17];
        dense.value[18] += dense_array[i].value[18];
        dense.value[19] += dense_array[i].value[19];
        dense.value[20] += dense_array[i].value[20];
        dense.value[21] += dense_array[i].value[21];
        dense.value[22] += dense_array[i].value[22];
        dense.value[23] += dense_array[i].value[23];
        dense.value[24] += dense_array[i].value[24];
        dense.value[25] += dense_array[i].value[25];
        dense.value[26] += dense_array[i].value[26];
        dense.value[27] += dense_array[i].value[27];
        dense.value[28] += dense_array[i].value[28];
        dense.value[29] += dense_array[i].value[29];
        dense.value[30] += dense_array[i].value[30];
        dense.value[31] += dense_array[i].value[31];
    }

    return dense;
}




static gen1_DenseMultivector gen1_dense_add(gen1_DenseMultivector dense0, gen1_DenseMultivector dense1, int sign){
    gen1_DenseMultivector dense = {{0}};
    if(sign == -1){
        dense.value[0] = dense0.value[0] - dense1.value[0];
        dense.value[1] = dense0.value[1] - dense1.value[1];
        dense.value[2] = dense0.value[2] - dense1.value[2];
        dense.value[3] = dense0.value[3] - dense1.value[3];
        dense.value[4] = dense0.value[4] - dense1.value[4];
        dense.value[5] = dense0.value[5] - dense1.value[5];
        dense.value[6] = dense0.value[6] - dense1.value[6];
        dense.value[7] = dense0.value[7] - dense1.value[7];
        dense.value[8] = dense0.value[8] - dense1.value[8];
        dense.value[9] = dense0.value[9] - dense1.value[9];
        dense.value[10] = dense0.value[10] - dense1.value[10];
        dense.value[11] = dense0.value[11] - dense1.value[11];
        dense.value[12] = dense0.value[12] - dense1.value[12];
        dense.value[13] = dense0.value[13] - dense1.value[13];
        dense.value[14] = dense0.value[14] - dense1.value[14];
        dense.value[15] = dense0.value[15] - dense1.value[15];
        dense.value[16] = dense0.value[16] - dense1.value[16];
        dense.value[17] = dense0.value[17] - dense1.value[17];
        dense.value[18] = dense0.value[18] - dense1.value[18];
        dense.value[19] = dense0.value[19] - dense1.value[19];
        dense.value[20] = dense0.value[20] - dense1.value[20];
        dense.value[21] = dense0.value[21] - dense1.value[21];
        dense.value[22] = dense0.value[22] - dense1.value[22];
        dense.value[23] = dense0.value[23] - dense1.value[23];
        dense.value[24] = dense0.value[24] - dense1.value[24];
        dense.value[25] = dense0.value[25] - dense1.value[25];
        dense.value[26] = dense0.value[26] - dense1.value[26];
        dense.value[27] = dense0.value[27] - dense1.value[27];
        dense.value[28] = dense0.value[28] - dense1.value[28];
        dense.value[29] = dense0.value[29] - dense1.value[29];
        dense.value[30] = dense0.value[30] - dense1.value[30];
        dense.value[31] = dense0.value[31] - dense1.value[31];
    }else if(sign == 1){
        dense.value[0] = dense0.value[0] + dense1.value[0];
        dense.value[1] = dense0.value[1] + dense1.value[1];
        dense.value[2] = dense0.value[2] + dense1.value[2];
        dense.value[3] = dense0.value[3] + dense1.value[3];
        dense.value[4] = dense0.value[4] + dense1.value[4];
        dense.value[5] = dense0.value[5] + dense1.value[5];
        dense.value[6] = dense0.value[6] + dense1.value[6];
        dense.value[7] = dense0.value[7] + dense1.value[7];
        dense.value[8] = dense0.value[8] + dense1.value[8];
        dense.value[9] = dense0.value[9] + dense1.value[9];
        dense.value[10] = dense0.value[10] + dense1.value[10];
        dense.value[11] = dense0.value[11] + dense1.value[11];
        dense.value[12] = dense0.value[12] + dense1.value[12];
        dense.value[13] = dense0.value[13] + dense1.value[13];
        dense.value[14] = dense0.value[14] + dense1.value[14];
        dense.value[15] = dense0.value[15] + dense1.value[15];
        dense.value[16] = dense0.value[16] + dense1.value[16];
        dense.value[17] = dense0.value[17] + dense1.value[17];
        dense.value[18] = dense0.value[18] + dense1.value[18];
        dense.value[19] = dense0.value[19] + dense1.value[19];
        dense.value[20] = dense0.value[20] + dense1.value[20];
        dense.value[21] = dense0.value[21] + dense1.value[21];
        dense.value[22] = dense0.value[22] + dense1.value[22];
        dense.value[23] = dense0.value[23] + dense1.value[23];
        dense.value[24] = dense0.value[24] + dense1.value[24];
        dense.value[25] = dense0.value[25] + dense1.value[25];
        dense.value[26] = dense0.value[26] + dense1.value[26];
        dense.value[27] = dense0.value[27] + dense1.value[27];
        dense.value[28] = dense0.value[28] + dense1.value[28];
        dense.value[29] = dense0.value[29] + dense1.value[29];
        dense.value[30] = dense0.value[30] + dense1.value[30];
        dense.value[31] = dense0.value[31] + dense1.value[31];
    } else{
        dense.value[0] = dense0.value[0] + sign*dense1.value[0];
        dense.value[1] = dense0.value[1] + sign*dense1.value[1];
        dense.value[2] = dense0.value[2] + sign*dense1.value[2];
        dense.value[3] = dense0.value[3] + sign*dense1.value[3];
        dense.value[4] = dense0.value[4] + sign*dense1.value[4];
        dense.value[5] = dense0.value[5] + sign*dense1.value[5];
        dense.value[6] = dense0.value[6] + sign*dense1.value[6];
        dense.value[7] = dense0.value[7] + sign*dense1.value[7];
        dense.value[8] = dense0.value[8] + sign*dense1.value[8];
        dense.value[9] = dense0.value[9] + sign*dense1.value[9];
        dense.value[10] = dense0.value[10] + sign*dense1.value[10];
        dense.value[11] = dense0.value[11] + sign*dense1.value[11];
        dense.value[12] = dense0.value[12] + sign*dense1.value[12];
        dense.value[13] = dense0.value[13] + sign*dense1.value[13];
        dense.value[14] = dense0.value[14] + sign*dense1.value[14];
        dense.value[15] = dense0.value[15] + sign*dense1.value[15];
        dense.value[16] = dense0.value[16] + sign*dense1.value[16];
        dense.value[17] = dense0.value[17] + sign*dense1.value[17];
        dense.value[18] = dense0.value[18] + sign*dense1.value[18];
        dense.value[19] = dense0.value[19] + sign*dense1.value[19];
        dense.value[20] = dense0.value[20] + sign*dense1.value[20];
        dense.value[21] = dense0.value[21] + sign*dense1.value[21];
        dense.value[22] = dense0.value[22] + sign*dense1.value[22];
        dense.value[23] = dense0.value[23] + sign*dense1.value[23];
        dense.value[24] = dense0.value[24] + sign*dense1.value[24];
        dense.value[25] = dense0.value[25] + sign*dense1.value[25];
        dense.value[26] = dense0.value[26] + sign*dense1.value[26];
        dense.value[27] = dense0.value[27] + sign*dense1.value[27];
        dense.value[28] = dense0.value[28] + sign*dense1.value[28];
        dense.value[29] = dense0.value[29] + sign*dense1.value[29];
        dense.value[30] = dense0.value[30] + sign*dense1.value[30];
        dense.value[31] = dense0.value[31] + sign*dense1.value[31];
    }
    return dense;
}


static gen1_DenseMultivector gen1_dense_scalaradd(gen1_DenseMultivector dense0, ga_float value, int sign){
    gen1_DenseMultivector dense = {{0}};
    if(sign == -1){
        dense.value[0] = -dense0.value[0];
        dense.value[1] = -dense0.value[1];
        dense.value[2] = -dense0.value[2];
        dense.value[3] = -dense0.value[3];
        dense.value[4] = -dense0.value[4];
        dense.value[5] = -dense0.value[5];
        dense.value[6] = -dense0.value[6];
        dense.value[7] = -dense0.value[7];
        dense.value[8] = -dense0.value[8];
        dense.value[9] = -dense0.value[9];
        dense.value[10] = -dense0.value[10];
        dense.value[11] = -dense0.value[11];
        dense.value[12] = -dense0.value[12];
        dense.value[13] = -dense0.value[13];
        dense.value[14] = -dense0.value[14];
        dense.value[15] = -dense0.value[15];
        dense.value[16] = -dense0.value[16];
        dense.value[17] = -dense0.value[17];
        dense.value[18] = -dense0.value[18];
        dense.value[19] = -dense0.value[19];
        dense.value[20] = -dense0.value[20];
        dense.value[21] = -dense0.value[21];
        dense.value[22] = -dense0.value[22];
        dense.value[23] = -dense0.value[23];
        dense.value[24] = -dense0.value[24];
        dense.value[25] = -dense0.value[25];
        dense.value[26] = -dense0.value[26];
        dense.value[27] = -dense0.value[27];
        dense.value[28] = -dense0.value[28];
        dense.value[29] = -dense0.value[29];
        dense.value[30] = -dense0.value[30];
        dense.value[31] = -dense0.value[31];
    }else if(sign == 1){
        dense.value[0] = dense0.value[0];
        dense.value[1] = dense0.value[1];
        dense.value[2] = dense0.value[2];
        dense.value[3] = dense0.value[3];
        dense.value[4] = dense0.value[4];
        dense.value[5] = dense0.value[5];
        dense.value[6] = dense0.value[6];
        dense.value[7] = dense0.value[7];
        dense.value[8] = dense0.value[8];
        dense.value[9] = dense0.value[9];
        dense.value[10] = dense0.value[10];
        dense.value[11] = dense0.value[11];
        dense.value[12] = dense0.value[12];
        dense.value[13] = dense0.value[13];
        dense.value[14] = dense0.value[14];
        dense.value[15] = dense0.value[15];
        dense.value[16] = dense0.value[16];
        dense.value[17] = dense0.value[17];
        dense.value[18] = dense0.value[18];
        dense.value[19] = dense0.value[19];
        dense.value[20] = dense0.value[20];
        dense.value[21] = dense0.value[21];
        dense.value[22] = dense0.value[22];
        dense.value[23] = dense0.value[23];
        dense.value[24] = dense0.value[24];
        dense.value[25] = dense0.value[25];
        dense.value[26] = dense0.value[26];
        dense.value[27] = dense0.value[27];
        dense.value[28] = dense0.value[28];
        dense.value[29] = dense0.value[29];
        dense.value[30] = dense0.value[30];
        dense.value[31] = dense0.value[31];
    } else{
        dense.value[0] = sign*dense0.value[0];
        dense.value[1] = sign*dense0.value[1];
        dense.value[2] = sign*dense0.value[2];
        dense.value[3] = sign*dense0.value[3];
        dense.value[4] = sign*dense0.value[4];
        dense.value[5] = sign*dense0.value[5];
        dense.value[6] = sign*dense0.value[6];
        dense.value[7] = sign*dense0.value[7];
        dense.value[8] = sign*dense0.value[8];
        dense.value[9] = sign*dense0.value[9];
        dense.value[10] = sign*dense0.value[10];
        dense.value[11] = sign*dense0.value[11];
        dense.value[12] = sign*dense0.value[12];
        dense.value[13] = sign*dense0.value[13];
        dense.value[14] = sign*dense0.value[14];
        dense.value[15] = sign*dense0.value[15];
        dense.value[16] = sign*dense0.value[16];
        dense.value[17] = sign*dense0.value[17];
        dense.value[18] = sign*dense0.value[18];
        dense.value[19] = sign*dense0.value[19];
        dense.value[20] = sign*dense0.value[20];
        dense.value[21] = sign*dense0.value[21];
        dense.value[22] = sign*dense0.value[22];
        dense.value[23] = sign*dense0.value[23];
        dense.value[24] = sign*dense0.value[24];
        dense.value[25] = sign*dense0.value[25];
        dense.value[26] = sign*dense0.value[26];
        dense.value[27] = sign*dense0.value[27];
        dense.value[28] = sign*dense0.value[28];
        dense.value[29] = sign*dense0.value[29];
        dense.value[30] = sign*dense0.value[30];
        dense.value[31] = sign*dense0.value[31];
    }
    dense.value[0] += value;
    return dense;
}

static gen1_DenseMultivector gen1_dense_scalarproduct(gen1_DenseMultivector dense0, ga_float value){
    gen1_DenseMultivector dense = {{0}};

    dense.value[0] = value*dense0.value[0];
    dense.value[1] = value*dense0.value[1];
    dense.value[2] = value*dense0.value[2];
    dense.value[3] = value*dense0.value[3];
    dense.value[4] = value*dense0.value[4];
    dense.value[5] = value*dense0.value[5];
    dense.value[6] = value*dense0.value[6];
    dense.value[7] = value*dense0.value[7];
    dense.value[8] = value*dense0.value[8];
    dense.value[9] = value*dense0.value[9];
    dense.value[10] = value*dense0.value[10];
    dense.value[11] = value*dense0.value[11];
    dense.value[12] = value*dense0.value[12];
    dense.value[13] = value*dense0.value[13];
    dense.value[14] = value*dense0.value[14];
    dense.value[15] = value*dense0.value[15];
    dense.value[16] = value*dense0.value[16];
    dense.value[17] = value*dense0.value[17];
    dense.value[18] = value*dense0.value[18];
    dense.value[19] = value*dense0.value[19];
    dense.value[20] = value*dense0.value[20];
    dense.value[21] = value*dense0.value[21];
    dense.value[22] = value*dense0.value[22];
    dense.value[23] = value*dense0.value[23];
    dense.value[24] = value*dense0.value[24];
    dense.value[25] = value*dense0.value[25];
    dense.value[26] = value*dense0.value[26];
    dense.value[27] = value*dense0.value[27];
    dense.value[28] = value*dense0.value[28];
    dense.value[29] = value*dense0.value[29];
    dense.value[30] = value*dense0.value[30];
    dense.value[31] = value*dense0.value[31];

    return dense;
}

static gen1_DenseMultivector gen1_dense_reverse(gen1_DenseMultivector dense0){
    gen1_DenseMultivector dense = {{0}};

    dense.value[0] = dense0.value[0];
    dense.value[1] = dense0.value[1];
    dense.value[2] = dense0.value[2];
    dense.value[3] = -dense0.value[3];
    dense.value[4] = dense0.value[4];
    dense.value[5] = -dense0.value[5];
    dense.value[6] = -dense0.value[6];
    dense.value[7] = -dense0.value[7];
    dense.value[8] = dense0.value[8];
    dense.value[9] = -dense0.value[9];
    dense.value[10] = -dense0.value[10];
    dense.value[11] = -dense0.value[11];
    dense.value[12] = -dense0.value[12];
    dense.value[13] = -dense0.value[13];
    dense.value[14] = -dense0.value[14];
    dense.value[15] = dense0.value[15];
    dense.value[16] = dense0.value[16];
    dense.value[17] = -dense0.value[17];
    dense.value[18] = -dense0.value[18];
    dense.value[19] = -dense0.value[19];
    dense.value[20] = -dense0.value[20];
    dense.value[21] = -dense0.value[21];
    dense.value[22] = -dense0.value[22];
    dense.value[23] = dense0.value[23];
    dense.value[24] = -dense0.value[24];
    dense.value[25] = -dense0.value[25];
    dense.value[26] = -dense0.value[26];
    dense.value[27] = dense0.value[27];
    dense.value[28] = -dense0.value[28];
    dense.value[29] = dense0.value[29];
    dense.value[30] = dense0.value[30];
    dense.value[31] = dense0.value[31];

    return dense;
}



static gen1_BladesMultivector gen1_blades_geometricproduct(gen1_BladesMultivector blades0, gen1_BladesMultivector blades1){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    blades.value0[0] =
    +blades0.value0[0]*blades1.value0[0]
    +blades0.value1[0]*blades1.value1[0]
    +blades0.value1[1]*blades1.value1[1]
    -blades0.value2[0]*blades1.value2[0]
    +blades0.value1[2]*blades1.value1[2]
    -blades0.value2[1]*blades1.value2[1]
    -blades0.value2[2]*blades1.value2[2]
    -blades0.value3[0]*blades1.value3[0]
    +blades0.value1[3]*blades1.value1[3]
    -blades0.value2[3]*blades1.value2[3]
    -blades0.value2[4]*blades1.value2[4]
    -blades0.value3[1]*blades1.value3[1]
    -blades0.value2[5]*blades1.value2[5]
    -blades0.value3[2]*blades1.value3[2]
    -blades0.value3[3]*blades1.value3[3]
    +blades0.value4[0]*blades1.value4[0]
    -blades0.value1[4]*blades1.value1[4]
    +blades0.value2[6]*blades1.value2[6]
    +blades0.value2[7]*blades1.value2[7]
    +blades0.value3[4]*blades1.value3[4]
    +blades0.value2[8]*blades1.value2[8]
    +blades0.value3[5]*blades1.value3[5]
    +blades0.value3[6]*blades1.value3[6]
    -blades0.value4[1]*blades1.value4[1]
    +blades0.value2[9]*blades1.value2[9]
    +blades0.value3[7]*blades1.value3[7]
    +blades0.value3[8]*blades1.value3[8]
    -blades0.value4[2]*blades1.value4[2]
    +blades0.value3[9]*blades1.value3[9]
    -blades0.value4[3]*blades1.value4[3]
    -blades0.value4[4]*blades1.value4[4]
    -blades0.value5[0]*blades1.value5[0]
;
    blades.value1[0] =
    +blades0.value0[0]*blades1.value1[0]
    +blades0.value1[0]*blades1.value0[0]
    -blades0.value1[1]*blades1.value2[0]
    +blades0.value2[0]*blades1.value1[1]
    -blades0.value1[2]*blades1.value2[1]
    +blades0.value2[1]*blades1.value1[2]
    -blades0.value2[2]*blades1.value3[0]
    -blades0.value3[0]*blades1.value2[2]
    -blades0.value1[3]*blades1.value2[3]
    +blades0.value2[3]*blades1.value1[3]
    -blades0.value2[4]*blades1.value3[1]
    -blades0.value3[1]*blades1.value2[4]
    -blades0.value2[5]*blades1.value3[2]
    -blades0.value3[2]*blades1.value2[5]
    +blades0.value3[3]*blades1.value4[0]
    -blades0.value4[0]*blades1.value3[3]
    +blades0.value1[4]*blades1.value2[6]
    -blades0.value2[6]*blades1.value1[4]
    +blades0.value2[7]*blades1.value3[4]
    +blades0.value3[4]*blades1.value2[7]
    +blades0.value2[8]*blades1.value3[5]
    +blades0.value3[5]*blades1.value2[8]
    -blades0.value3[6]*blades1.value4[1]
    +blades0.value4[1]*blades1.value3[6]
    +blades0.value2[9]*blades1.value3[7]
    +blades0.value3[7]*blades1.value2[9]
    -blades0.value3[8]*blades1.value4[2]
    +blades0.value4[2]*blades1.value3[8]
    -blades0.value3[9]*blades1.value4[3]
    +blades0.value4[3]*blades1.value3[9]
    -blades0.value4[4]*blades1.value5[0]
    -blades0.value5[0]*blades1.value4[4]
;
    blades.value1[1] =
    +blades0.value0[0]*blades1.value1[1]
    +blades0.value1[0]*blades1.value2[0]
    +blades0.value1[1]*blades1.value0[0]
    -blades0.value2[0]*blades1.value1[0]
    -blades0.value1[2]*blades1.value2[2]
    +blades0.value2[1]*blades1.value3[0]
    +blades0.value2[2]*blades1.value1[2]
    +blades0.value3[0]*blades1.value2[1]
    -blades0.value1[3]*blades1.value2[4]
    +blades0.value2[3]*blades1.value3[1]
    +blades0.value2[4]*blades1.value1[3]
    +blades0.value3[1]*blades1.value2[3]
    -blades0.value2[5]*blades1.value3[3]
    -blades0.value3[2]*blades1.value4[0]
    -blades0.value3[3]*blades1.value2[5]
    +blades0.value4[0]*blades1.value3[2]
    +blades0.value1[4]*blades1.value2[7]
    -blades0.value2[6]*blades1.value3[4]
    -blades0.value2[7]*blades1.value1[4]
    -blades0.value3[4]*blades1.value2[6]
    +blades0.value2[8]*blades1.value3[6]
    +blades0.value3[5]*blades1.value4[1]
    +blades0.value3[6]*blades1.value2[8]
    -blades0.value4[1]*blades1.value3[5]
    +blades0.value2[9]*blades1.value3[8]
    +blades0.value3[7]*blades1.value4[2]
    +blades0.value3[8]*blades1.value2[9]
    -blades0.value4[2]*blades1.value3[7]
    -blades0.value3[9]*blades1.value4[4]
    +blades0.value4[3]*blades1.value5[0]
    +blades0.value4[4]*blades1.value3[9]
    +blades0.value5[0]*blades1.value4[3]
;
    blades.value2[0] =
    +blades0.value0[0]*blades1.value2[0]
    +blades0.value1[0]*blades1.value1[1]
    -blades0.value1[1]*blades1.value1[0]
    +blades0.value2[0]*blades1.value0[0]
    +blades0.value1[2]*blades1.value3[0]
    -blades0.value2[1]*blades1.value2[2]
    +blades0.value2[2]*blades1.value2[1]
    +blades0.value3[0]*blades1.value1[2]
    +blades0.value1[3]*blades1.value3[1]
    -blades0.value2[3]*blades1.value2[4]
    +blades0.value2[4]*blades1.value2[3]
    +blades0.value3[1]*blades1.value1[3]
    -blades0.value2[5]*blades1.value4[0]
    -blades0.value3[2]*blades1.value3[3]
    +blades0.value3[3]*blades1.value3[2]
    -blades0.value4[0]*blades1.value2[5]
    -blades0.value1[4]*blades1.value3[4]
    +blades0.value2[6]*blades1.value2[7]
    -blades0.value2[7]*blades1.value2[6]
    -blades0.value3[4]*blades1.value1[4]
    +blades0.value2[8]*blades1.value4[1]
    +blades0.value3[5]*blades1.value3[6]
    -blades0.value3[6]*blades1.value3[5]
    +blades0.value4[1]*blades1.value2[8]
    +blades0.value2[9]*blades1.value4[2]
    +blades0.value3[7]*blades1.value3[8]
    -blades0.value3[8]*blades1.value3[7]
    +blades0.value4[2]*blades1.value2[9]
    +blades0.value3[9]*blades1.value5[0]
    -blades0.value4[3]*blades1.value4[4]
    +blades0.value4[4]*blades1.value4[3]
    +blades0.value5[0]*blades1.value3[9]
;
    blades.value1[2] =
    +blades0.value0[0]*blades1.value1[2]
    +blades0.value1[0]*blades1.value2[1]
    +blades0.value1[1]*blades1.value2[2]
    -blades0.value2[0]*blades1.value3[0]
    +blades0.value1[2]*blades1.value0[0]
    -blades0.value2[1]*blades1.value1[0]
    -blades0.value2[2]*blades1.value1[1]
    -blades0.value3[0]*blades1.value2[0]
    -blades0.value1[3]*blades1.value2[5]
    +blades0.value2[3]*blades1.value3[2]
    +blades0.value2[4]*blades1.value3[3]
    +blades0.value3[1]*blades1.value4[0]
    +blades0.value2[5]*blades1.value1[3]
    +blades0.value3[2]*blades1.value2[3]
    +blades0.value3[3]*blades1.value2[4]
    -blades0.value4[0]*blades1.value3[1]
    +blades0.value1[4]*blades1.value2[8]
    -blades0.value2[6]*blades1.value3[5]
    -blades0.value2[7]*blades1.value3[6]
    -blades0.value3[4]*blades1.value4[1]
    -blades0.value2[8]*blades1.value1[4]
    -blades0.value3[5]*blades1.value2[6]
    -blades0.value3[6]*blades1.value2[7]
    +blades0.value4[1]*blades1.value3[4]
    +blades0.value2[9]*blades1.value3[9]
    +blades0.value3[7]*blades1.value4[3]
    +blades0.value3[8]*blades1.value4[4]
    -blades0.value4[2]*blades1.value5[0]
    +blades0.value3[9]*blades1.value2[9]
    -blades0.value4[3]*blades1.value3[7]
    -blades0.value4[4]*blades1.value3[8]
    -blades0.value5[0]*blades1.value4[2]
;
    blades.value2[1] =
    +blades0.value0[0]*blades1.value2[1]
    +blades0.value1[0]*blades1.value1[2]
    -blades0.value1[1]*blades1.value3[0]
    +blades0.value2[0]*blades1.value2[2]
    -blades0.value1[2]*blades1.value1[0]
    +blades0.value2[1]*blades1.value0[0]
    -blades0.value2[2]*blades1.value2[0]
    -blades0.value3[0]*blades1.value1[1]
    +blades0.value1[3]*blades1.value3[2]
    -blades0.value2[3]*blades1.value2[5]
    +blades0.value2[4]*blades1.value4[0]
    +blades0.value3[1]*blades1.value3[3]
    +blades0.value2[5]*blades1.value2[3]
    +blades0.value3[2]*blades1.value1[3]
    -blades0.value3[3]*blades1.value3[1]
    +blades0.value4[0]*blades1.value2[4]
    -blades0.value1[4]*blades1.value3[5]
    +blades0.value2[6]*blades1.value2[8]
    -blades0.value2[7]*blades1.value4[1]
    -blades0.value3[4]*blades1.value3[6]
    -blades0.value2[8]*blades1.value2[6]
    -blades0.value3[5]*blades1.value1[4]
    +blades0.value3[6]*blades1.value3[4]
    -blades0.value4[1]*blades1.value2[7]
    +blades0.value2[9]*blades1.value4[3]
    +blades0.value3[7]*blades1.value3[9]
    -blades0.value3[8]*blades1.value5[0]
    +blades0.value4[2]*blades1.value4[4]
    -blades0.value3[9]*blades1.value3[7]
    +blades0.value4[3]*blades1.value2[9]
    -blades0.value4[4]*blades1.value4[2]
    -blades0.value5[0]*blades1.value3[8]
;
    blades.value2[2] =
    +blades0.value0[0]*blades1.value2[2]
    +blades0.value1[0]*blades1.value3[0]
    +blades0.value1[1]*blades1.value1[2]
    -blades0.value2[0]*blades1.value2[1]
    -blades0.value1[2]*blades1.value1[1]
    +blades0.value2[1]*blades1.value2[0]
    +blades0.value2[2]*blades1.value0[0]
    +blades0.value3[0]*blades1.value1[0]
    +blades0.value1[3]*blades1.value3[3]
    -blades0.value2[3]*blades1.value4[0]
    -blades0.value2[4]*blades1.value2[5]
    -blades0.value3[1]*blades1.value3[2]
    +blades0.value2[5]*blades1.value2[4]
    +blades0.value3[2]*blades1.value3[1]
    +blades0.value3[3]*blades1.value1[3]
    -blades0.value4[0]*blades1.value2[3]
    -blades0.value1[4]*blades1.value3[6]
    +blades0.value2[6]*blades1.value4[1]
    +blades0.value2[7]*blades1.value2[8]
    +blades0.value3[4]*blades1.value3[5]
    -blades0.value2[8]*blades1.value2[7]
    -blades0.value3[5]*blades1.value3[4]
    -blades0.value3[6]*blades1.value1[4]
    +blades0.value4[1]*blades1.value2[6]
    +blades0.value2[9]*blades1.value4[4]
    +blades0.value3[7]*blades1.value5[0]
    +blades0.value3[8]*blades1.value3[9]
    -blades0.value4[2]*blades1.value4[3]
    -blades0.value3[9]*blades1.value3[8]
    +blades0.value4[3]*blades1.value4[2]
    +blades0.value4[4]*blades1.value2[9]
    +blades0.value5[0]*blades1.value3[7]
;
    blades.value3[0] =
    +blades0.value0[0]*blades1.value3[0]
    +blades0.value1[0]*blades1.value2[2]
    -blades0.value1[1]*blades1.value2[1]
    +blades0.value2[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value2[0]
    -blades0.value2[1]*blades1.value1[1]
    +blades0.value2[2]*blades1.value1[0]
    +blades0.value3[0]*blades1.value0[0]
    -blades0.value1[3]*blades1.value4[0]
    +blades0.value2[3]*blades1.value3[3]
    -blades0.value2[4]*blades1.value3[2]
    -blades0.value3[1]*blades1.value2[5]
    +blades0.value2[5]*blades1.value3[1]
    +blades0.value3[2]*blades1.value2[4]
    -blades0.value3[3]*blades1.value2[3]
    +blades0.value4[0]*blades1.value1[3]
    +blades0.value1[4]*blades1.value4[1]
    -blades0.value2[6]*blades1.value3[6]
    +blades0.value2[7]*blades1.value3[5]
    +blades0.value3[4]*blades1.value2[8]
    -blades0.value2[8]*blades1.value3[4]
    -blades0.value3[5]*blades1.value2[7]
    +blades0.value3[6]*blades1.value2[6]
    -blades0.value4[1]*blades1.value1[4]
    +blades0.value2[9]*blades1.value5[0]
    +blades0.value3[7]*blades1.value4[4]
    -blades0.value3[8]*blades1.value4[3]
    +blades0.value4[2]*blades1.value3[9]
    +blades0.value3[9]*blades1.value4[2]
    -blades0.value4[3]*blades1.value3[8]
    +blades0.value4[4]*blades1.value3[7]
    +blades0.value5[0]*blades1.value2[9]
;
    blades.value1[3] =
    +blades0.value0[0]*blades1.value1[3]
    +blades0.value1[0]*blades1.value2[3]
    +blades0.value1[1]*blades1.value2[4]
    -blades0.value2[0]*blades1.value3[1]
    +blades0.value1[2]*blades1.value2[5]
    -blades0.value2[1]*blades1.value3[2]
    -blades0.value2[2]*blades1.value3[3]
    -blades0.value3[0]*blades1.value4[0]
    +blades0.value1[3]*blades1.value0[0]
    -blades0.value2[3]*blades1.value1[0]
    -blades0.value2[4]*blades1.value1[1]
    -blades0.value3[1]*blades1.value2[0]
    -blades0.value2[5]*blades1.value1[2]
    -blades0.value3[2]*blades1.value2[1]
    -blades0.value3[3]*blades1.value2[2]
    +blades0.value4[0]*blades1.value3[0]
    +blades0.value1[4]*blades1.value2[9]
    -blades0.value2[6]*blades1.value3[7]
    -blades0.value2[7]*blades1.value3[8]
    -blades0.value3[4]*blades1.value4[2]
    -blades0.value2[8]*blades1.value3[9]
    -blades0.value3[5]*blades1.value4[3]
    -blades0.value3[6]*blades1.value4[4]
    +blades0.value4[1]*blades1.value5[0]
    -blades0.value2[9]*blades1.value1[4]
    -blades0.value3[7]*blades1.value2[6]
    -blades0.value3[8]*blades1.value2[7]
    +blades0.value4[2]*blades1.value3[4]
    -blades0.value3[9]*blades1.value2[8]
    +blades0.value4[3]*blades1.value3[5]
    +blades0.value4[4]*blades1.value3[6]
    +blades0.value5[0]*blades1.value4[1]
;
    blades.value2[3] =
    +blades0.value0[0]*blades1.value2[3]
    +blades0.value1[0]*blades1.value1[3]
    -blades0.value1[1]*blades1.value3[1]
    +blades0.value2[0]*blades1.value2[4]
    -blades0.value1[2]*blades1.value3[2]
    +blades0.value2[1]*blades1.value2[5]
    -blades0.value2[2]*blades1.value4[0]
    -blades0.value3[0]*blades1.value3[3]
    -blades0.value1[3]*blades1.value1[0]
    +blades0.value2[3]*blades1.value0[0]
    -blades0.value2[4]*blades1.value2[0]
    -blades0.value3[1]*blades1.value1[1]
    -blades0.value2[5]*blades1.value2[1]
    -blades0.value3[2]*blades1.value1[2]
    +blades0.value3[3]*blades1.value3[0]
    -blades0.value4[0]*blades1.value2[2]
    -blades0.value1[4]*blades1.value3[7]
    +blades0.value2[6]*blades1.value2[9]
    -blades0.value2[7]*blades1.value4[2]
    -blades0.value3[4]*blades1.value3[8]
    -blades0.value2[8]*blades1.value4[3]
    -blades0.value3[5]*blades1.value3[9]
    +blades0.value3[6]*blades1.value5[0]
    -blades0.value4[1]*blades1.value4[4]
    -blades0.value2[9]*blades1.value2[6]
    -blades0.value3[7]*blades1.value1[4]
    +blades0.value3[8]*blades1.value3[4]
    -blades0.value4[2]*blades1.value2[7]
    +blades0.value3[9]*blades1.value3[5]
    -blades0.value4[3]*blades1.value2[8]
    +blades0.value4[4]*blades1.value4[1]
    +blades0.value5[0]*blades1.value3[6]
;
    blades.value2[4] =
    +blades0.value0[0]*blades1.value2[4]
    +blades0.value1[0]*blades1.value3[1]
    +blades0.value1[1]*blades1.value1[3]
    -blades0.value2[0]*blades1.value2[3]
    -blades0.value1[2]*blades1.value3[3]
    +blades0.value2[1]*blades1.value4[0]
    +blades0.value2[2]*blades1.value2[5]
    +blades0.value3[0]*blades1.value3[2]
    -blades0.value1[3]*blades1.value1[1]
    +blades0.value2[3]*blades1.value2[0]
    +blades0.value2[4]*blades1.value0[0]
    +blades0.value3[1]*blades1.value1[0]
    -blades0.value2[5]*blades1.value2[2]
    -blades0.value3[2]*blades1.value3[0]
    -blades0.value3[3]*blades1.value1[2]
    +blades0.value4[0]*blades1.value2[1]
    -blades0.value1[4]*blades1.value3[8]
    +blades0.value2[6]*blades1.value4[2]
    +blades0.value2[7]*blades1.value2[9]
    +blades0.value3[4]*blades1.value3[7]
    -blades0.value2[8]*blades1.value4[4]
    -blades0.value3[5]*blades1.value5[0]
    -blades0.value3[6]*blades1.value3[9]
    +blades0.value4[1]*blades1.value4[3]
    -blades0.value2[9]*blades1.value2[7]
    -blades0.value3[7]*blades1.value3[4]
    -blades0.value3[8]*blades1.value1[4]
    +blades0.value4[2]*blades1.value2[6]
    +blades0.value3[9]*blades1.value3[6]
    -blades0.value4[3]*blades1.value4[1]
    -blades0.value4[4]*blades1.value2[8]
    -blades0.value5[0]*blades1.value3[5]
;
    blades.value3[1] =
    +blades0.value0[0]*blades1.value3[1]
    +blades0.value1[0]*blades1.value2[4]
    -blades0.value1[1]*blades1.value2[3]
    +blades0.value2[0]*blades1.value1[3]
    +blades0.value1[2]*blades1.value4[0]
    -blades0.value2[1]*blades1.value3[3]
    +blades0.value2[2]*blades1.value3[2]
    +blades0.value3[0]*blades1.value2[5]
    +blades0.value1[3]*blades1.value2[0]
    -blades0.value2[3]*blades1.value1[1]
    +blades0.value2[4]*blades1.value1[0]
    +blades0.value3[1]*blades1.value0[0]
    -blades0.value2[5]*blades1.value3[0]
    -blades0.value3[2]*blades1.value2[2]
    +blades0.value3[3]*blades1.value2[1]
    -blades0.value4[0]*blades1.value1[2]
    +blades0.value1[4]*blades1.value4[2]
    -blades0.value2[6]*blades1.value3[8]
    +blades0.value2[7]*blades1.value3[7]
    +blades0.value3[4]*blades1.value2[9]
    -blades0.value2[8]*blades1.value5[0]
    -blades0.value3[5]*blades1.value4[4]
    +blades0.value3[6]*blades1.value4[3]
    -blades0.value4[1]*blades1.value3[9]
    -blades0.value2[9]*blades1.value3[4]
    -blades0.value3[7]*blades1.value2[7]
    +blades0.value3[8]*blades1.value2[6]
    -blades0.value4[2]*blades1.value1[4]
    -blades0.value3[9]*blades1.value4[1]
    +blades0.value4[3]*blades1.value3[6]
    -blades0.value4[4]*blades1.value3[5]
    -blades0.value5[0]*blades1.value2[8]
;
    blades.value2[5] =
    +blades0.value0[0]*blades1.value2[5]
    +blades0.value1[0]*blades1.value3[2]
    +blades0.value1[1]*blades1.value3[3]
    -blades0.value2[0]*blades1.value4[0]
    +blades0.value1[2]*blades1.value1[3]
    -blades0.value2[1]*blades1.value2[3]
    -blades0.value2[2]*blades1.value2[4]
    -blades0.value3[0]*blades1.value3[1]
    -blades0.value1[3]*blades1.value1[2]
    +blades0.value2[3]*blades1.value2[1]
    +blades0.value2[4]*blades1.value2[2]
    +blades0.value3[1]*blades1.value3[0]
    +blades0.value2[5]*blades1.value0[0]
    +blades0.value3[2]*blades1.value1[0]
    +blades0.value3[3]*blades1.value1[1]
    -blades0.value4[0]*blades1.value2[0]
    -blades0.value1[4]*blades1.value3[9]
    +blades0.value2[6]*blades1.value4[3]
    +blades0.value2[7]*blades1.value4[4]
    +blades0.value3[4]*blades1.value5[0]
    +blades0.value2[8]*blades1.value2[9]
    +blades0.value3[5]*blades1.value3[7]
    +blades0.value3[6]*blades1.value3[8]
    -blades0.value4[1]*blades1.value4[2]
    -blades0.value2[9]*blades1.value2[8]
    -blades0.value3[7]*blades1.value3[5]
    -blades0.value3[8]*blades1.value3[6]
    +blades0.value4[2]*blades1.value4[1]
    -blades0.value3[9]*blades1.value1[4]
    +blades0.value4[3]*blades1.value2[6]
    +blades0.value4[4]*blades1.value2[7]
    +blades0.value5[0]*blades1.value3[4]
;
    blades.value3[2] =
    +blades0.value0[0]*blades1.value3[2]
    +blades0.value1[0]*blades1.value2[5]
    -blades0.value1[1]*blades1.value4[0]
    +blades0.value2[0]*blades1.value3[3]
    -blades0.value1[2]*blades1.value2[3]
    +blades0.value2[1]*blades1.value1[3]
    -blades0.value2[2]*blades1.value3[1]
    -blades0.value3[0]*blades1.value2[4]
    +blades0.value1[3]*blades1.value2[1]
    -blades0.value2[3]*blades1.value1[2]
    +blades0.value2[4]*blades1.value3[0]
    +blades0.value3[1]*blades1.value2[2]
    +blades0.value2[5]*blades1.value1[0]
    +blades0.value3[2]*blades1.value0[0]
    -blades0.value3[3]*blades1.value2[0]
    +blades0.value4[0]*blades1.value1[1]
    +blades0.value1[4]*blades1.value4[3]
    -blades0.value2[6]*blades1.value3[9]
    +blades0.value2[7]*blades1.value5[0]
    +blades0.value3[4]*blades1.value4[4]
    +blades0.value2[8]*blades1.value3[7]
    +blades0.value3[5]*blades1.value2[9]
    -blades0.value3[6]*blades1.value4[2]
    +blades0.value4[1]*blades1.value3[8]
    -blades0.value2[9]*blades1.value3[5]
    -blades0.value3[7]*blades1.value2[8]
    +blades0.value3[8]*blades1.value4[1]
    -blades0.value4[2]*blades1.value3[6]
    +blades0.value3[9]*blades1.value2[6]
    -blades0.value4[3]*blades1.value1[4]
    +blades0.value4[4]*blades1.value3[4]
    +blades0.value5[0]*blades1.value2[7]
;
    blades.value3[3] =
    +blades0.value0[0]*blades1.value3[3]
    +blades0.value1[0]*blades1.value4[0]
    +blades0.value1[1]*blades1.value2[5]
    -blades0.value2[0]*blades1.value3[2]
    -blades0.value1[2]*blades1.value2[4]
    +blades0.value2[1]*blades1.value3[1]
    +blades0.value2[2]*blades1.value1[3]
    +blades0.value3[0]*blades1.value2[3]
    +blades0.value1[3]*blades1.value2[2]
    -blades0.value2[3]*blades1.value3[0]
    -blades0.value2[4]*blades1.value1[2]
    -blades0.value3[1]*blades1.value2[1]
    +blades0.value2[5]*blades1.value1[1]
    +blades0.value3[2]*blades1.value2[0]
    +blades0.value3[3]*blades1.value0[0]
    -blades0.value4[0]*blades1.value1[0]
    +blades0.value1[4]*blades1.value4[4]
    -blades0.value2[6]*blades1.value5[0]
    -blades0.value2[7]*blades1.value3[9]
    -blades0.value3[4]*blades1.value4[3]
    +blades0.value2[8]*blades1.value3[8]
    +blades0.value3[5]*blades1.value4[2]
    +blades0.value3[6]*blades1.value2[9]
    -blades0.value4[1]*blades1.value3[7]
    -blades0.value2[9]*blades1.value3[6]
    -blades0.value3[7]*blades1.value4[1]
    -blades0.value3[8]*blades1.value2[8]
    +blades0.value4[2]*blades1.value3[5]
    +blades0.value3[9]*blades1.value2[7]
    -blades0.value4[3]*blades1.value3[4]
    -blades0.value4[4]*blades1.value1[4]
    -blades0.value5[0]*blades1.value2[6]
;
    blades.value4[0] =
    +blades0.value0[0]*blades1.value4[0]
    +blades0.value1[0]*blades1.value3[3]
    -blades0.value1[1]*blades1.value3[2]
    +blades0.value2[0]*blades1.value2[5]
    +blades0.value1[2]*blades1.value3[1]
    -blades0.value2[1]*blades1.value2[4]
    +blades0.value2[2]*blades1.value2[3]
    +blades0.value3[0]*blades1.value1[3]
    -blades0.value1[3]*blades1.value3[0]
    +blades0.value2[3]*blades1.value2[2]
    -blades0.value2[4]*blades1.value2[1]
    -blades0.value3[1]*blades1.value1[2]
    +blades0.value2[5]*blades1.value2[0]
    +blades0.value3[2]*blades1.value1[1]
    -blades0.value3[3]*blades1.value1[0]
    +blades0.value4[0]*blades1.value0[0]
    -blades0.value1[4]*blades1.value5[0]
    +blades0.value2[6]*blades1.value4[4]
    -blades0.value2[7]*blades1.value4[3]
    -blades0.value3[4]*blades1.value3[9]
    +blades0.value2[8]*blades1.value4[2]
    +blades0.value3[5]*blades1.value3[8]
    -blades0.value3[6]*blades1.value3[7]
    +blades0.value4[1]*blades1.value2[9]
    -blades0.value2[9]*blades1.value4[1]
    -blades0.value3[7]*blades1.value3[6]
    +blades0.value3[8]*blades1.value3[5]
    -blades0.value4[2]*blades1.value2[8]
    -blades0.value3[9]*blades1.value3[4]
    +blades0.value4[3]*blades1.value2[7]
    -blades0.value4[4]*blades1.value2[6]
    -blades0.value5[0]*blades1.value1[4]
;
    blades.value1[4] =
    +blades0.value0[0]*blades1.value1[4]
    +blades0.value1[0]*blades1.value2[6]
    +blades0.value1[1]*blades1.value2[7]
    -blades0.value2[0]*blades1.value3[4]
    +blades0.value1[2]*blades1.value2[8]
    -blades0.value2[1]*blades1.value3[5]
    -blades0.value2[2]*blades1.value3[6]
    -blades0.value3[0]*blades1.value4[1]
    +blades0.value1[3]*blades1.value2[9]
    -blades0.value2[3]*blades1.value3[7]
    -blades0.value2[4]*blades1.value3[8]
    -blades0.value3[1]*blades1.value4[2]
    -blades0.value2[5]*blades1.value3[9]
    -blades0.value3[2]*blades1.value4[3]
    -blades0.value3[3]*blades1.value4[4]
    +blades0.value4[0]*blades1.value5[0]
    +blades0.value1[4]*blades1.value0[0]
    -blades0.value2[6]*blades1.value1[0]
    -blades0.value2[7]*blades1.value1[1]
    -blades0.value3[4]*blades1.value2[0]
    -blades0.value2[8]*blades1.value1[2]
    -blades0.value3[5]*blades1.value2[1]
    -blades0.value3[6]*blades1.value2[2]
    +blades0.value4[1]*blades1.value3[0]
    -blades0.value2[9]*blades1.value1[3]
    -blades0.value3[7]*blades1.value2[3]
    -blades0.value3[8]*blades1.value2[4]
    +blades0.value4[2]*blades1.value3[1]
    -blades0.value3[9]*blades1.value2[5]
    +blades0.value4[3]*blades1.value3[2]
    +blades0.value4[4]*blades1.value3[3]
    +blades0.value5[0]*blades1.value4[0]
;
    blades.value2[6] =
    +blades0.value0[0]*blades1.value2[6]
    +blades0.value1[0]*blades1.value1[4]
    -blades0.value1[1]*blades1.value3[4]
    +blades0.value2[0]*blades1.value2[7]
    -blades0.value1[2]*blades1.value3[5]
    +blades0.value2[1]*blades1.value2[8]
    -blades0.value2[2]*blades1.value4[1]
    -blades0.value3[0]*blades1.value3[6]
    -blades0.value1[3]*blades1.value3[7]
    +blades0.value2[3]*blades1.value2[9]
    -blades0.value2[4]*blades1.value4[2]
    -blades0.value3[1]*blades1.value3[8]
    -blades0.value2[5]*blades1.value4[3]
    -blades0.value3[2]*blades1.value3[9]
    +blades0.value3[3]*blades1.value5[0]
    -blades0.value4[0]*blades1.value4[4]
    -blades0.value1[4]*blades1.value1[0]
    +blades0.value2[6]*blades1.value0[0]
    -blades0.value2[7]*blades1.value2[0]
    -blades0.value3[4]*blades1.value1[1]
    -blades0.value2[8]*blades1.value2[1]
    -blades0.value3[5]*blades1.value1[2]
    +blades0.value3[6]*blades1.value3[0]
    -blades0.value4[1]*blades1.value2[2]
    -blades0.value2[9]*blades1.value2[3]
    -blades0.value3[7]*blades1.value1[3]
    +blades0.value3[8]*blades1.value3[1]
    -blades0.value4[2]*blades1.value2[4]
    +blades0.value3[9]*blades1.value3[2]
    -blades0.value4[3]*blades1.value2[5]
    +blades0.value4[4]*blades1.value4[0]
    +blades0.value5[0]*blades1.value3[3]
;
    blades.value2[7] =
    +blades0.value0[0]*blades1.value2[7]
    +blades0.value1[0]*blades1.value3[4]
    +blades0.value1[1]*blades1.value1[4]
    -blades0.value2[0]*blades1.value2[6]
    -blades0.value1[2]*blades1.value3[6]
    +blades0.value2[1]*blades1.value4[1]
    +blades0.value2[2]*blades1.value2[8]
    +blades0.value3[0]*blades1.value3[5]
    -blades0.value1[3]*blades1.value3[8]
    +blades0.value2[3]*blades1.value4[2]
    +blades0.value2[4]*blades1.value2[9]
    +blades0.value3[1]*blades1.value3[7]
    -blades0.value2[5]*blades1.value4[4]
    -blades0.value3[2]*blades1.value5[0]
    -blades0.value3[3]*blades1.value3[9]
    +blades0.value4[0]*blades1.value4[3]
    -blades0.value1[4]*blades1.value1[1]
    +blades0.value2[6]*blades1.value2[0]
    +blades0.value2[7]*blades1.value0[0]
    +blades0.value3[4]*blades1.value1[0]
    -blades0.value2[8]*blades1.value2[2]
    -blades0.value3[5]*blades1.value3[0]
    -blades0.value3[6]*blades1.value1[2]
    +blades0.value4[1]*blades1.value2[1]
    -blades0.value2[9]*blades1.value2[4]
    -blades0.value3[7]*blades1.value3[1]
    -blades0.value3[8]*blades1.value1[3]
    +blades0.value4[2]*blades1.value2[3]
    +blades0.value3[9]*blades1.value3[3]
    -blades0.value4[3]*blades1.value4[0]
    -blades0.value4[4]*blades1.value2[5]
    -blades0.value5[0]*blades1.value3[2]
;
    blades.value3[4] =
    +blades0.value0[0]*blades1.value3[4]
    +blades0.value1[0]*blades1.value2[7]
    -blades0.value1[1]*blades1.value2[6]
    +blades0.value2[0]*blades1.value1[4]
    +blades0.value1[2]*blades1.value4[1]
    -blades0.value2[1]*blades1.value3[6]
    +blades0.value2[2]*blades1.value3[5]
    +blades0.value3[0]*blades1.value2[8]
    +blades0.value1[3]*blades1.value4[2]
    -blades0.value2[3]*blades1.value3[8]
    +blades0.value2[4]*blades1.value3[7]
    +blades0.value3[1]*blades1.value2[9]
    -blades0.value2[5]*blades1.value5[0]
    -blades0.value3[2]*blades1.value4[4]
    +blades0.value3[3]*blades1.value4[3]
    -blades0.value4[0]*blades1.value3[9]
    +blades0.value1[4]*blades1.value2[0]
    -blades0.value2[6]*blades1.value1[1]
    +blades0.value2[7]*blades1.value1[0]
    +blades0.value3[4]*blades1.value0[0]
    -blades0.value2[8]*blades1.value3[0]
    -blades0.value3[5]*blades1.value2[2]
    +blades0.value3[6]*blades1.value2[1]
    -blades0.value4[1]*blades1.value1[2]
    -blades0.value2[9]*blades1.value3[1]
    -blades0.value3[7]*blades1.value2[4]
    +blades0.value3[8]*blades1.value2[3]
    -blades0.value4[2]*blades1.value1[3]
    -blades0.value3[9]*blades1.value4[0]
    +blades0.value4[3]*blades1.value3[3]
    -blades0.value4[4]*blades1.value3[2]
    -blades0.value5[0]*blades1.value2[5]
;
    blades.value2[8] =
    +blades0.value0[0]*blades1.value2[8]
    +blades0.value1[0]*blades1.value3[5]
    +blades0.value1[1]*blades1.value3[6]
    -blades0.value2[0]*blades1.value4[1]
    +blades0.value1[2]*blades1.value1[4]
    -blades0.value2[1]*blades1.value2[6]
    -blades0.value2[2]*blades1.value2[7]
    -blades0.value3[0]*blades1.value3[4]
    -blades0.value1[3]*blades1.value3[9]
    +blades0.value2[3]*blades1.value4[3]
    +blades0.value2[4]*blades1.value4[4]
    +blades0.value3[1]*blades1.value5[0]
    +blades0.value2[5]*blades1.value2[9]
    +blades0.value3[2]*blades1.value3[7]
    +blades0.value3[3]*blades1.value3[8]
    -blades0.value4[0]*blades1.value4[2]
    -blades0.value1[4]*blades1.value1[2]
    +blades0.value2[6]*blades1.value2[1]
    +blades0.value2[7]*blades1.value2[2]
    +blades0.value3[4]*blades1.value3[0]
    +blades0.value2[8]*blades1.value0[0]
    +blades0.value3[5]*blades1.value1[0]
    +blades0.value3[6]*blades1.value1[1]
    -blades0.value4[1]*blades1.value2[0]
    -blades0.value2[9]*blades1.value2[5]
    -blades0.value3[7]*blades1.value3[2]
    -blades0.value3[8]*blades1.value3[3]
    +blades0.value4[2]*blades1.value4[0]
    -blades0.value3[9]*blades1.value1[3]
    +blades0.value4[3]*blades1.value2[3]
    +blades0.value4[4]*blades1.value2[4]
    +blades0.value5[0]*blades1.value3[1]
;
    blades.value3[5] =
    +blades0.value0[0]*blades1.value3[5]
    +blades0.value1[0]*blades1.value2[8]
    -blades0.value1[1]*blades1.value4[1]
    +blades0.value2[0]*blades1.value3[6]
    -blades0.value1[2]*blades1.value2[6]
    +blades0.value2[1]*blades1.value1[4]
    -blades0.value2[2]*blades1.value3[4]
    -blades0.value3[0]*blades1.value2[7]
    +blades0.value1[3]*blades1.value4[3]
    -blades0.value2[3]*blades1.value3[9]
    +blades0.value2[4]*blades1.value5[0]
    +blades0.value3[1]*blades1.value4[4]
    +blades0.value2[5]*blades1.value3[7]
    +blades0.value3[2]*blades1.value2[9]
    -blades0.value3[3]*blades1.value4[2]
    +blades0.value4[0]*blades1.value3[8]
    +blades0.value1[4]*blades1.value2[1]
    -blades0.value2[6]*blades1.value1[2]
    +blades0.value2[7]*blades1.value3[0]
    +blades0.value3[4]*blades1.value2[2]
    +blades0.value2[8]*blades1.value1[0]
    +blades0.value3[5]*blades1.value0[0]
    -blades0.value3[6]*blades1.value2[0]
    +blades0.value4[1]*blades1.value1[1]
    -blades0.value2[9]*blades1.value3[2]
    -blades0.value3[7]*blades1.value2[5]
    +blades0.value3[8]*blades1.value4[0]
    -blades0.value4[2]*blades1.value3[3]
    +blades0.value3[9]*blades1.value2[3]
    -blades0.value4[3]*blades1.value1[3]
    +blades0.value4[4]*blades1.value3[1]
    +blades0.value5[0]*blades1.value2[4]
;
    blades.value3[6] =
    +blades0.value0[0]*blades1.value3[6]
    +blades0.value1[0]*blades1.value4[1]
    +blades0.value1[1]*blades1.value2[8]
    -blades0.value2[0]*blades1.value3[5]
    -blades0.value1[2]*blades1.value2[7]
    +blades0.value2[1]*blades1.value3[4]
    +blades0.value2[2]*blades1.value1[4]
    +blades0.value3[0]*blades1.value2[6]
    +blades0.value1[3]*blades1.value4[4]
    -blades0.value2[3]*blades1.value5[0]
    -blades0.value2[4]*blades1.value3[9]
    -blades0.value3[1]*blades1.value4[3]
    +blades0.value2[5]*blades1.value3[8]
    +blades0.value3[2]*blades1.value4[2]
    +blades0.value3[3]*blades1.value2[9]
    -blades0.value4[0]*blades1.value3[7]
    +blades0.value1[4]*blades1.value2[2]
    -blades0.value2[6]*blades1.value3[0]
    -blades0.value2[7]*blades1.value1[2]
    -blades0.value3[4]*blades1.value2[1]
    +blades0.value2[8]*blades1.value1[1]
    +blades0.value3[5]*blades1.value2[0]
    +blades0.value3[6]*blades1.value0[0]
    -blades0.value4[1]*blades1.value1[0]
    -blades0.value2[9]*blades1.value3[3]
    -blades0.value3[7]*blades1.value4[0]
    -blades0.value3[8]*blades1.value2[5]
    +blades0.value4[2]*blades1.value3[2]
    +blades0.value3[9]*blades1.value2[4]
    -blades0.value4[3]*blades1.value3[1]
    -blades0.value4[4]*blades1.value1[3]
    -blades0.value5[0]*blades1.value2[3]
;
    blades.value4[1] =
    +blades0.value0[0]*blades1.value4[1]
    +blades0.value1[0]*blades1.value3[6]
    -blades0.value1[1]*blades1.value3[5]
    +blades0.value2[0]*blades1.value2[8]
    +blades0.value1[2]*blades1.value3[4]
    -blades0.value2[1]*blades1.value2[7]
    +blades0.value2[2]*blades1.value2[6]
    +blades0.value3[0]*blades1.value1[4]
    -blades0.value1[3]*blades1.value5[0]
    +blades0.value2[3]*blades1.value4[4]
    -blades0.value2[4]*blades1.value4[3]
    -blades0.value3[1]*blades1.value3[9]
    +blades0.value2[5]*blades1.value4[2]
    +blades0.value3[2]*blades1.value3[8]
    -blades0.value3[3]*blades1.value3[7]
    +blades0.value4[0]*blades1.value2[9]
    -blades0.value1[4]*blades1.value3[0]
    +blades0.value2[6]*blades1.value2[2]
    -blades0.value2[7]*blades1.value2[1]
    -blades0.value3[4]*blades1.value1[2]
    +blades0.value2[8]*blades1.value2[0]
    +blades0.value3[5]*blades1.value1[1]
    -blades0.value3[6]*blades1.value1[0]
    +blades0.value4[1]*blades1.value0[0]
    -blades0.value2[9]*blades1.value4[0]
    -blades0.value3[7]*blades1.value3[3]
    +blades0.value3[8]*blades1.value3[2]
    -blades0.value4[2]*blades1.value2[5]
    -blades0.value3[9]*blades1.value3[1]
    +blades0.value4[3]*blades1.value2[4]
    -blades0.value4[4]*blades1.value2[3]
    -blades0.value5[0]*blades1.value1[3]
;
    blades.value2[9] =
    +blades0.value0[0]*blades1.value2[9]
    +blades0.value1[0]*blades1.value3[7]
    +blades0.value1[1]*blades1.value3[8]
    -blades0.value2[0]*blades1.value4[2]
    +blades0.value1[2]*blades1.value3[9]
    -blades0.value2[1]*blades1.value4[3]
    -blades0.value2[2]*blades1.value4[4]
    -blades0.value3[0]*blades1.value5[0]
    +blades0.value1[3]*blades1.value1[4]
    -blades0.value2[3]*blades1.value2[6]
    -blades0.value2[4]*blades1.value2[7]
    -blades0.value3[1]*blades1.value3[4]
    -blades0.value2[5]*blades1.value2[8]
    -blades0.value3[2]*blades1.value3[5]
    -blades0.value3[3]*blades1.value3[6]
    +blades0.value4[0]*blades1.value4[1]
    -blades0.value1[4]*blades1.value1[3]
    +blades0.value2[6]*blades1.value2[3]
    +blades0.value2[7]*blades1.value2[4]
    +blades0.value3[4]*blades1.value3[1]
    +blades0.value2[8]*blades1.value2[5]
    +blades0.value3[5]*blades1.value3[2]
    +blades0.value3[6]*blades1.value3[3]
    -blades0.value4[1]*blades1.value4[0]
    +blades0.value2[9]*blades1.value0[0]
    +blades0.value3[7]*blades1.value1[0]
    +blades0.value3[8]*blades1.value1[1]
    -blades0.value4[2]*blades1.value2[0]
    +blades0.value3[9]*blades1.value1[2]
    -blades0.value4[3]*blades1.value2[1]
    -blades0.value4[4]*blades1.value2[2]
    -blades0.value5[0]*blades1.value3[0]
;
    blades.value3[7] =
    +blades0.value0[0]*blades1.value3[7]
    +blades0.value1[0]*blades1.value2[9]
    -blades0.value1[1]*blades1.value4[2]
    +blades0.value2[0]*blades1.value3[8]
    -blades0.value1[2]*blades1.value4[3]
    +blades0.value2[1]*blades1.value3[9]
    -blades0.value2[2]*blades1.value5[0]
    -blades0.value3[0]*blades1.value4[4]
    -blades0.value1[3]*blades1.value2[6]
    +blades0.value2[3]*blades1.value1[4]
    -blades0.value2[4]*blades1.value3[4]
    -blades0.value3[1]*blades1.value2[7]
    -blades0.value2[5]*blades1.value3[5]
    -blades0.value3[2]*blades1.value2[8]
    +blades0.value3[3]*blades1.value4[1]
    -blades0.value4[0]*blades1.value3[6]
    +blades0.value1[4]*blades1.value2[3]
    -blades0.value2[6]*blades1.value1[3]
    +blades0.value2[7]*blades1.value3[1]
    +blades0.value3[4]*blades1.value2[4]
    +blades0.value2[8]*blades1.value3[2]
    +blades0.value3[5]*blades1.value2[5]
    -blades0.value3[6]*blades1.value4[0]
    +blades0.value4[1]*blades1.value3[3]
    +blades0.value2[9]*blades1.value1[0]
    +blades0.value3[7]*blades1.value0[0]
    -blades0.value3[8]*blades1.value2[0]
    +blades0.value4[2]*blades1.value1[1]
    -blades0.value3[9]*blades1.value2[1]
    +blades0.value4[3]*blades1.value1[2]
    -blades0.value4[4]*blades1.value3[0]
    -blades0.value5[0]*blades1.value2[2]
;
    blades.value3[8] =
    +blades0.value0[0]*blades1.value3[8]
    +blades0.value1[0]*blades1.value4[2]
    +blades0.value1[1]*blades1.value2[9]
    -blades0.value2[0]*blades1.value3[7]
    -blades0.value1[2]*blades1.value4[4]
    +blades0.value2[1]*blades1.value5[0]
    +blades0.value2[2]*blades1.value3[9]
    +blades0.value3[0]*blades1.value4[3]
    -blades0.value1[3]*blades1.value2[7]
    +blades0.value2[3]*blades1.value3[4]
    +blades0.value2[4]*blades1.value1[4]
    +blades0.value3[1]*blades1.value2[6]
    -blades0.value2[5]*blades1.value3[6]
    -blades0.value3[2]*blades1.value4[1]
    -blades0.value3[3]*blades1.value2[8]
    +blades0.value4[0]*blades1.value3[5]
    +blades0.value1[4]*blades1.value2[4]
    -blades0.value2[6]*blades1.value3[1]
    -blades0.value2[7]*blades1.value1[3]
    -blades0.value3[4]*blades1.value2[3]
    +blades0.value2[8]*blades1.value3[3]
    +blades0.value3[5]*blades1.value4[0]
    +blades0.value3[6]*blades1.value2[5]
    -blades0.value4[1]*blades1.value3[2]
    +blades0.value2[9]*blades1.value1[1]
    +blades0.value3[7]*blades1.value2[0]
    +blades0.value3[8]*blades1.value0[0]
    -blades0.value4[2]*blades1.value1[0]
    -blades0.value3[9]*blades1.value2[2]
    +blades0.value4[3]*blades1.value3[0]
    +blades0.value4[4]*blades1.value1[2]
    +blades0.value5[0]*blades1.value2[1]
;
    blades.value4[2] =
    +blades0.value0[0]*blades1.value4[2]
    +blades0.value1[0]*blades1.value3[8]
    -blades0.value1[1]*blades1.value3[7]
    +blades0.value2[0]*blades1.value2[9]
    +blades0.value1[2]*blades1.value5[0]
    -blades0.value2[1]*blades1.value4[4]
    +blades0.value2[2]*blades1.value4[3]
    +blades0.value3[0]*blades1.value3[9]
    +blades0.value1[3]*blades1.value3[4]
    -blades0.value2[3]*blades1.value2[7]
    +blades0.value2[4]*blades1.value2[6]
    +blades0.value3[1]*blades1.value1[4]
    -blades0.value2[5]*blades1.value4[1]
    -blades0.value3[2]*blades1.value3[6]
    +blades0.value3[3]*blades1.value3[5]
    -blades0.value4[0]*blades1.value2[8]
    -blades0.value1[4]*blades1.value3[1]
    +blades0.value2[6]*blades1.value2[4]
    -blades0.value2[7]*blades1.value2[3]
    -blades0.value3[4]*blades1.value1[3]
    +blades0.value2[8]*blades1.value4[0]
    +blades0.value3[5]*blades1.value3[3]
    -blades0.value3[6]*blades1.value3[2]
    +blades0.value4[1]*blades1.value2[5]
    +blades0.value2[9]*blades1.value2[0]
    +blades0.value3[7]*blades1.value1[1]
    -blades0.value3[8]*blades1.value1[0]
    +blades0.value4[2]*blades1.value0[0]
    +blades0.value3[9]*blades1.value3[0]
    -blades0.value4[3]*blades1.value2[2]
    +blades0.value4[4]*blades1.value2[1]
    +blades0.value5[0]*blades1.value1[2]
;
    blades.value3[9] =
    +blades0.value0[0]*blades1.value3[9]
    +blades0.value1[0]*blades1.value4[3]
    +blades0.value1[1]*blades1.value4[4]
    -blades0.value2[0]*blades1.value5[0]
    +blades0.value1[2]*blades1.value2[9]
    -blades0.value2[1]*blades1.value3[7]
    -blades0.value2[2]*blades1.value3[8]
    -blades0.value3[0]*blades1.value4[2]
    -blades0.value1[3]*blades1.value2[8]
    +blades0.value2[3]*blades1.value3[5]
    +blades0.value2[4]*blades1.value3[6]
    +blades0.value3[1]*blades1.value4[1]
    +blades0.value2[5]*blades1.value1[4]
    +blades0.value3[2]*blades1.value2[6]
    +blades0.value3[3]*blades1.value2[7]
    -blades0.value4[0]*blades1.value3[4]
    +blades0.value1[4]*blades1.value2[5]
    -blades0.value2[6]*blades1.value3[2]
    -blades0.value2[7]*blades1.value3[3]
    -blades0.value3[4]*blades1.value4[0]
    -blades0.value2[8]*blades1.value1[3]
    -blades0.value3[5]*blades1.value2[3]
    -blades0.value3[6]*blades1.value2[4]
    +blades0.value4[1]*blades1.value3[1]
    +blades0.value2[9]*blades1.value1[2]
    +blades0.value3[7]*blades1.value2[1]
    +blades0.value3[8]*blades1.value2[2]
    -blades0.value4[2]*blades1.value3[0]
    +blades0.value3[9]*blades1.value0[0]
    -blades0.value4[3]*blades1.value1[0]
    -blades0.value4[4]*blades1.value1[1]
    -blades0.value5[0]*blades1.value2[0]
;
    blades.value4[3] =
    +blades0.value0[0]*blades1.value4[3]
    +blades0.value1[0]*blades1.value3[9]
    -blades0.value1[1]*blades1.value5[0]
    +blades0.value2[0]*blades1.value4[4]
    -blades0.value1[2]*blades1.value3[7]
    +blades0.value2[1]*blades1.value2[9]
    -blades0.value2[2]*blades1.value4[2]
    -blades0.value3[0]*blades1.value3[8]
    +blades0.value1[3]*blades1.value3[5]
    -blades0.value2[3]*blades1.value2[8]
    +blades0.value2[4]*blades1.value4[1]
    +blades0.value3[1]*blades1.value3[6]
    +blades0.value2[5]*blades1.value2[6]
    +blades0.value3[2]*blades1.value1[4]
    -blades0.value3[3]*blades1.value3[4]
    +blades0.value4[0]*blades1.value2[7]
    -blades0.value1[4]*blades1.value3[2]
    +blades0.value2[6]*blades1.value2[5]
    -blades0.value2[7]*blades1.value4[0]
    -blades0.value3[4]*blades1.value3[3]
    -blades0.value2[8]*blades1.value2[3]
    -blades0.value3[5]*blades1.value1[3]
    +blades0.value3[6]*blades1.value3[1]
    -blades0.value4[1]*blades1.value2[4]
    +blades0.value2[9]*blades1.value2[1]
    +blades0.value3[7]*blades1.value1[2]
    -blades0.value3[8]*blades1.value3[0]
    +blades0.value4[2]*blades1.value2[2]
    -blades0.value3[9]*blades1.value1[0]
    +blades0.value4[3]*blades1.value0[0]
    -blades0.value4[4]*blades1.value2[0]
    -blades0.value5[0]*blades1.value1[1]
;
    blades.value4[4] =
    +blades0.value0[0]*blades1.value4[4]
    +blades0.value1[0]*blades1.value5[0]
    +blades0.value1[1]*blades1.value3[9]
    -blades0.value2[0]*blades1.value4[3]
    -blades0.value1[2]*blades1.value3[8]
    +blades0.value2[1]*blades1.value4[2]
    +blades0.value2[2]*blades1.value2[9]
    +blades0.value3[0]*blades1.value3[7]
    +blades0.value1[3]*blades1.value3[6]
    -blades0.value2[3]*blades1.value4[1]
    -blades0.value2[4]*blades1.value2[8]
    -blades0.value3[1]*blades1.value3[5]
    +blades0.value2[5]*blades1.value2[7]
    +blades0.value3[2]*blades1.value3[4]
    +blades0.value3[3]*blades1.value1[4]
    -blades0.value4[0]*blades1.value2[6]
    -blades0.value1[4]*blades1.value3[3]
    +blades0.value2[6]*blades1.value4[0]
    +blades0.value2[7]*blades1.value2[5]
    +blades0.value3[4]*blades1.value3[2]
    -blades0.value2[8]*blades1.value2[4]
    -blades0.value3[5]*blades1.value3[1]
    -blades0.value3[6]*blades1.value1[3]
    +blades0.value4[1]*blades1.value2[3]
    +blades0.value2[9]*blades1.value2[2]
    +blades0.value3[7]*blades1.value3[0]
    +blades0.value3[8]*blades1.value1[2]
    -blades0.value4[2]*blades1.value2[1]
    -blades0.value3[9]*blades1.value1[1]
    +blades0.value4[3]*blades1.value2[0]
    +blades0.value4[4]*blades1.value0[0]
    +blades0.value5[0]*blades1.value1[0]
;
    blades.value5[0] =
    +blades0.value0[0]*blades1.value5[0]
    +blades0.value1[0]*blades1.value4[4]
    -blades0.value1[1]*blades1.value4[3]
    +blades0.value2[0]*blades1.value3[9]
    +blades0.value1[2]*blades1.value4[2]
    -blades0.value2[1]*blades1.value3[8]
    +blades0.value2[2]*blades1.value3[7]
    +blades0.value3[0]*blades1.value2[9]
    -blades0.value1[3]*blades1.value4[1]
    +blades0.value2[3]*blades1.value3[6]
    -blades0.value2[4]*blades1.value3[5]
    -blades0.value3[1]*blades1.value2[8]
    +blades0.value2[5]*blades1.value3[4]
    +blades0.value3[2]*blades1.value2[7]
    -blades0.value3[3]*blades1.value2[6]
    +blades0.value4[0]*blades1.value1[4]
    +blades0.value1[4]*blades1.value4[0]
    -blades0.value2[6]*blades1.value3[3]
    +blades0.value2[7]*blades1.value3[2]
    +blades0.value3[4]*blades1.value2[5]
    -blades0.value2[8]*blades1.value3[1]
    -blades0.value3[5]*blades1.value2[4]
    +blades0.value3[6]*blades1.value2[3]
    -blades0.value4[1]*blades1.value1[3]
    +blades0.value2[9]*blades1.value3[0]
    +blades0.value3[7]*blades1.value2[2]
    -blades0.value3[8]*blades1.value2[1]
    +blades0.value4[2]*blades1.value1[2]
    +blades0.value3[9]*blades1.value2[0]
    -blades0.value4[3]*blades1.value1[1]
    +blades0.value4[4]*blades1.value1[0]
    +blades0.value5[0]*blades1.value0[0]
;
    return blades;
}


static gen1_BladesMultivector gen1_blades_innerproduct(gen1_BladesMultivector blades0, gen1_BladesMultivector blades1){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    blades.value0[0] =
    +blades0.value1[0]*blades1.value1[0]
    +blades0.value1[1]*blades1.value1[1]
    -blades0.value2[0]*blades1.value2[0]
    +blades0.value1[2]*blades1.value1[2]
    -blades0.value2[1]*blades1.value2[1]
    -blades0.value2[2]*blades1.value2[2]
    -blades0.value3[0]*blades1.value3[0]
    +blades0.value1[3]*blades1.value1[3]
    -blades0.value2[3]*blades1.value2[3]
    -blades0.value2[4]*blades1.value2[4]
    -blades0.value3[1]*blades1.value3[1]
    -blades0.value2[5]*blades1.value2[5]
    -blades0.value3[2]*blades1.value3[2]
    -blades0.value3[3]*blades1.value3[3]
    +blades0.value4[0]*blades1.value4[0]
    -blades0.value1[4]*blades1.value1[4]
    +blades0.value2[6]*blades1.value2[6]
    +blades0.value2[7]*blades1.value2[7]
    +blades0.value3[4]*blades1.value3[4]
    +blades0.value2[8]*blades1.value2[8]
    +blades0.value3[5]*blades1.value3[5]
    +blades0.value3[6]*blades1.value3[6]
    -blades0.value4[1]*blades1.value4[1]
    +blades0.value2[9]*blades1.value2[9]
    +blades0.value3[7]*blades1.value3[7]
    +blades0.value3[8]*blades1.value3[8]
    -blades0.value4[2]*blades1.value4[2]
    +blades0.value3[9]*blades1.value3[9]
    -blades0.value4[3]*blades1.value4[3]
    -blades0.value4[4]*blades1.value4[4]
    -blades0.value5[0]*blades1.value5[0]
;
    blades.value1[0] =
    -blades0.value1[1]*blades1.value2[0]
    +blades0.value2[0]*blades1.value1[1]
    -blades0.value1[2]*blades1.value2[1]
    +blades0.value2[1]*blades1.value1[2]
    -blades0.value2[2]*blades1.value3[0]
    -blades0.value3[0]*blades1.value2[2]
    -blades0.value1[3]*blades1.value2[3]
    +blades0.value2[3]*blades1.value1[3]
    -blades0.value2[4]*blades1.value3[1]
    -blades0.value3[1]*blades1.value2[4]
    -blades0.value2[5]*blades1.value3[2]
    -blades0.value3[2]*blades1.value2[5]
    +blades0.value3[3]*blades1.value4[0]
    -blades0.value4[0]*blades1.value3[3]
    +blades0.value1[4]*blades1.value2[6]
    -blades0.value2[6]*blades1.value1[4]
    +blades0.value2[7]*blades1.value3[4]
    +blades0.value3[4]*blades1.value2[7]
    +blades0.value2[8]*blades1.value3[5]
    +blades0.value3[5]*blades1.value2[8]
    -blades0.value3[6]*blades1.value4[1]
    +blades0.value4[1]*blades1.value3[6]
    +blades0.value2[9]*blades1.value3[7]
    +blades0.value3[7]*blades1.value2[9]
    -blades0.value3[8]*blades1.value4[2]
    +blades0.value4[2]*blades1.value3[8]
    -blades0.value3[9]*blades1.value4[3]
    +blades0.value4[3]*blades1.value3[9]
    -blades0.value4[4]*blades1.value5[0]
    -blades0.value5[0]*blades1.value4[4]
;
    blades.value1[1] =
    +blades0.value1[0]*blades1.value2[0]
    -blades0.value2[0]*blades1.value1[0]
    -blades0.value1[2]*blades1.value2[2]
    +blades0.value2[1]*blades1.value3[0]
    +blades0.value2[2]*blades1.value1[2]
    +blades0.value3[0]*blades1.value2[1]
    -blades0.value1[3]*blades1.value2[4]
    +blades0.value2[3]*blades1.value3[1]
    +blades0.value2[4]*blades1.value1[3]
    +blades0.value3[1]*blades1.value2[3]
    -blades0.value2[5]*blades1.value3[3]
    -blades0.value3[2]*blades1.value4[0]
    -blades0.value3[3]*blades1.value2[5]
    +blades0.value4[0]*blades1.value3[2]
    +blades0.value1[4]*blades1.value2[7]
    -blades0.value2[6]*blades1.value3[4]
    -blades0.value2[7]*blades1.value1[4]
    -blades0.value3[4]*blades1.value2[6]
    +blades0.value2[8]*blades1.value3[6]
    +blades0.value3[5]*blades1.value4[1]
    +blades0.value3[6]*blades1.value2[8]
    -blades0.value4[1]*blades1.value3[5]
    +blades0.value2[9]*blades1.value3[8]
    +blades0.value3[7]*blades1.value4[2]
    +blades0.value3[8]*blades1.value2[9]
    -blades0.value4[2]*blades1.value3[7]
    -blades0.value3[9]*blades1.value4[4]
    +blades0.value4[3]*blades1.value5[0]
    +blades0.value4[4]*blades1.value3[9]
    +blades0.value5[0]*blades1.value4[3]
;
    blades.value2[0] =
    +blades0.value1[2]*blades1.value3[0]
    +blades0.value3[0]*blades1.value1[2]
    +blades0.value1[3]*blades1.value3[1]
    +blades0.value3[1]*blades1.value1[3]
    -blades0.value2[5]*blades1.value4[0]
    -blades0.value4[0]*blades1.value2[5]
    -blades0.value1[4]*blades1.value3[4]
    -blades0.value3[4]*blades1.value1[4]
    +blades0.value2[8]*blades1.value4[1]
    +blades0.value4[1]*blades1.value2[8]
    +blades0.value2[9]*blades1.value4[2]
    +blades0.value4[2]*blades1.value2[9]
    +blades0.value3[9]*blades1.value5[0]
    +blades0.value5[0]*blades1.value3[9]
;
    blades.value1[2] =
    +blades0.value1[0]*blades1.value2[1]
    +blades0.value1[1]*blades1.value2[2]
    -blades0.value2[0]*blades1.value3[0]
    -blades0.value2[1]*blades1.value1[0]
    -blades0.value2[2]*blades1.value1[1]
    -blades0.value3[0]*blades1.value2[0]
    -blades0.value1[3]*blades1.value2[5]
    +blades0.value2[3]*blades1.value3[2]
    +blades0.value2[4]*blades1.value3[3]
    +blades0.value3[1]*blades1.value4[0]
    +blades0.value2[5]*blades1.value1[3]
    +blades0.value3[2]*blades1.value2[3]
    +blades0.value3[3]*blades1.value2[4]
    -blades0.value4[0]*blades1.value3[1]
    +blades0.value1[4]*blades1.value2[8]
    -blades0.value2[6]*blades1.value3[5]
    -blades0.value2[7]*blades1.value3[6]
    -blades0.value3[4]*blades1.value4[1]
    -blades0.value2[8]*blades1.value1[4]
    -blades0.value3[5]*blades1.value2[6]
    -blades0.value3[6]*blades1.value2[7]
    +blades0.value4[1]*blades1.value3[4]
    +blades0.value2[9]*blades1.value3[9]
    +blades0.value3[7]*blades1.value4[3]
    +blades0.value3[8]*blades1.value4[4]
    -blades0.value4[2]*blades1.value5[0]
    +blades0.value3[9]*blades1.value2[9]
    -blades0.value4[3]*blades1.value3[7]
    -blades0.value4[4]*blades1.value3[8]
    -blades0.value5[0]*blades1.value4[2]
;
    blades.value2[1] =
    -blades0.value1[1]*blades1.value3[0]
    -blades0.value3[0]*blades1.value1[1]
    +blades0.value1[3]*blades1.value3[2]
    +blades0.value2[4]*blades1.value4[0]
    +blades0.value3[2]*blades1.value1[3]
    +blades0.value4[0]*blades1.value2[4]
    -blades0.value1[4]*blades1.value3[5]
    -blades0.value2[7]*blades1.value4[1]
    -blades0.value3[5]*blades1.value1[4]
    -blades0.value4[1]*blades1.value2[7]
    +blades0.value2[9]*blades1.value4[3]
    -blades0.value3[8]*blades1.value5[0]
    +blades0.value4[3]*blades1.value2[9]
    -blades0.value5[0]*blades1.value3[8]
;
    blades.value2[2] =
    +blades0.value1[0]*blades1.value3[0]
    +blades0.value3[0]*blades1.value1[0]
    +blades0.value1[3]*blades1.value3[3]
    -blades0.value2[3]*blades1.value4[0]
    +blades0.value3[3]*blades1.value1[3]
    -blades0.value4[0]*blades1.value2[3]
    -blades0.value1[4]*blades1.value3[6]
    +blades0.value2[6]*blades1.value4[1]
    -blades0.value3[6]*blades1.value1[4]
    +blades0.value4[1]*blades1.value2[6]
    +blades0.value2[9]*blades1.value4[4]
    +blades0.value3[7]*blades1.value5[0]
    +blades0.value4[4]*blades1.value2[9]
    +blades0.value5[0]*blades1.value3[7]
;
    blades.value3[0] =
    -blades0.value1[3]*blades1.value4[0]
    +blades0.value4[0]*blades1.value1[3]
    +blades0.value1[4]*blades1.value4[1]
    -blades0.value4[1]*blades1.value1[4]
    +blades0.value2[9]*blades1.value5[0]
    +blades0.value5[0]*blades1.value2[9]
;
    blades.value1[3] =
    +blades0.value1[0]*blades1.value2[3]
    +blades0.value1[1]*blades1.value2[4]
    -blades0.value2[0]*blades1.value3[1]
    +blades0.value1[2]*blades1.value2[5]
    -blades0.value2[1]*blades1.value3[2]
    -blades0.value2[2]*blades1.value3[3]
    -blades0.value3[0]*blades1.value4[0]
    -blades0.value2[3]*blades1.value1[0]
    -blades0.value2[4]*blades1.value1[1]
    -blades0.value3[1]*blades1.value2[0]
    -blades0.value2[5]*blades1.value1[2]
    -blades0.value3[2]*blades1.value2[1]
    -blades0.value3[3]*blades1.value2[2]
    +blades0.value4[0]*blades1.value3[0]
    +blades0.value1[4]*blades1.value2[9]
    -blades0.value2[6]*blades1.value3[7]
    -blades0.value2[7]*blades1.value3[8]
    -blades0.value3[4]*blades1.value4[2]
    -blades0.value2[8]*blades1.value3[9]
    -blades0.value3[5]*blades1.value4[3]
    -blades0.value3[6]*blades1.value4[4]
    +blades0.value4[1]*blades1.value5[0]
    -blades0.value2[9]*blades1.value1[4]
    -blades0.value3[7]*blades1.value2[6]
    -blades0.value3[8]*blades1.value2[7]
    +blades0.value4[2]*blades1.value3[4]
    -blades0.value3[9]*blades1.value2[8]
    +blades0.value4[3]*blades1.value3[5]
    +blades0.value4[4]*blades1.value3[6]
    +blades0.value5[0]*blades1.value4[1]
;
    blades.value2[3] =
    -blades0.value1[1]*blades1.value3[1]
    -blades0.value1[2]*blades1.value3[2]
    -blades0.value2[2]*blades1.value4[0]
    -blades0.value3[1]*blades1.value1[1]
    -blades0.value3[2]*blades1.value1[2]
    -blades0.value4[0]*blades1.value2[2]
    -blades0.value1[4]*blades1.value3[7]
    -blades0.value2[7]*blades1.value4[2]
    -blades0.value2[8]*blades1.value4[3]
    +blades0.value3[6]*blades1.value5[0]
    -blades0.value3[7]*blades1.value1[4]
    -blades0.value4[2]*blades1.value2[7]
    -blades0.value4[3]*blades1.value2[8]
    +blades0.value5[0]*blades1.value3[6]
;
    blades.value2[4] =
    +blades0.value1[0]*blades1.value3[1]
    -blades0.value1[2]*blades1.value3[3]
    +blades0.value2[1]*blades1.value4[0]
    +blades0.value3[1]*blades1.value1[0]
    -blades0.value3[3]*blades1.value1[2]
    +blades0.value4[0]*blades1.value2[1]
    -blades0.value1[4]*blades1.value3[8]
    +blades0.value2[6]*blades1.value4[2]
    -blades0.value2[8]*blades1.value4[4]
    -blades0.value3[5]*blades1.value5[0]
    -blades0.value3[8]*blades1.value1[4]
    +blades0.value4[2]*blades1.value2[6]
    -blades0.value4[4]*blades1.value2[8]
    -blades0.value5[0]*blades1.value3[5]
;
    blades.value3[1] =
    +blades0.value1[2]*blades1.value4[0]
    -blades0.value4[0]*blades1.value1[2]
    +blades0.value1[4]*blades1.value4[2]
    -blades0.value2[8]*blades1.value5[0]
    -blades0.value4[2]*blades1.value1[4]
    -blades0.value5[0]*blades1.value2[8]
;
    blades.value2[5] =
    +blades0.value1[0]*blades1.value3[2]
    +blades0.value1[1]*blades1.value3[3]
    -blades0.value2[0]*blades1.value4[0]
    +blades0.value3[2]*blades1.value1[0]
    +blades0.value3[3]*blades1.value1[1]
    -blades0.value4[0]*blades1.value2[0]
    -blades0.value1[4]*blades1.value3[9]
    +blades0.value2[6]*blades1.value4[3]
    +blades0.value2[7]*blades1.value4[4]
    +blades0.value3[4]*blades1.value5[0]
    -blades0.value3[9]*blades1.value1[4]
    +blades0.value4[3]*blades1.value2[6]
    +blades0.value4[4]*blades1.value2[7]
    +blades0.value5[0]*blades1.value3[4]
;
    blades.value3[2] =
    -blades0.value1[1]*blades1.value4[0]
    +blades0.value4[0]*blades1.value1[1]
    +blades0.value1[4]*blades1.value4[3]
    +blades0.value2[7]*blades1.value5[0]
    -blades0.value4[3]*blades1.value1[4]
    +blades0.value5[0]*blades1.value2[7]
;
    blades.value3[3] =
    +blades0.value1[0]*blades1.value4[0]
    -blades0.value4[0]*blades1.value1[0]
    +blades0.value1[4]*blades1.value4[4]
    -blades0.value2[6]*blades1.value5[0]
    -blades0.value4[4]*blades1.value1[4]
    -blades0.value5[0]*blades1.value2[6]
;
    blades.value4[0] =
    -blades0.value1[4]*blades1.value5[0]
    -blades0.value5[0]*blades1.value1[4]
;
    blades.value1[4] =
    +blades0.value1[0]*blades1.value2[6]
    +blades0.value1[1]*blades1.value2[7]
    -blades0.value2[0]*blades1.value3[4]
    +blades0.value1[2]*blades1.value2[8]
    -blades0.value2[1]*blades1.value3[5]
    -blades0.value2[2]*blades1.value3[6]
    -blades0.value3[0]*blades1.value4[1]
    +blades0.value1[3]*blades1.value2[9]
    -blades0.value2[3]*blades1.value3[7]
    -blades0.value2[4]*blades1.value3[8]
    -blades0.value3[1]*blades1.value4[2]
    -blades0.value2[5]*blades1.value3[9]
    -blades0.value3[2]*blades1.value4[3]
    -blades0.value3[3]*blades1.value4[4]
    +blades0.value4[0]*blades1.value5[0]
    -blades0.value2[6]*blades1.value1[0]
    -blades0.value2[7]*blades1.value1[1]
    -blades0.value3[4]*blades1.value2[0]
    -blades0.value2[8]*blades1.value1[2]
    -blades0.value3[5]*blades1.value2[1]
    -blades0.value3[6]*blades1.value2[2]
    +blades0.value4[1]*blades1.value3[0]
    -blades0.value2[9]*blades1.value1[3]
    -blades0.value3[7]*blades1.value2[3]
    -blades0.value3[8]*blades1.value2[4]
    +blades0.value4[2]*blades1.value3[1]
    -blades0.value3[9]*blades1.value2[5]
    +blades0.value4[3]*blades1.value3[2]
    +blades0.value4[4]*blades1.value3[3]
    +blades0.value5[0]*blades1.value4[0]
;
    blades.value2[6] =
    -blades0.value1[1]*blades1.value3[4]
    -blades0.value1[2]*blades1.value3[5]
    -blades0.value2[2]*blades1.value4[1]
    -blades0.value1[3]*blades1.value3[7]
    -blades0.value2[4]*blades1.value4[2]
    -blades0.value2[5]*blades1.value4[3]
    +blades0.value3[3]*blades1.value5[0]
    -blades0.value3[4]*blades1.value1[1]
    -blades0.value3[5]*blades1.value1[2]
    -blades0.value4[1]*blades1.value2[2]
    -blades0.value3[7]*blades1.value1[3]
    -blades0.value4[2]*blades1.value2[4]
    -blades0.value4[3]*blades1.value2[5]
    +blades0.value5[0]*blades1.value3[3]
;
    blades.value2[7] =
    +blades0.value1[0]*blades1.value3[4]
    -blades0.value1[2]*blades1.value3[6]
    +blades0.value2[1]*blades1.value4[1]
    -blades0.value1[3]*blades1.value3[8]
    +blades0.value2[3]*blades1.value4[2]
    -blades0.value2[5]*blades1.value4[4]
    -blades0.value3[2]*blades1.value5[0]
    +blades0.value3[4]*blades1.value1[0]
    -blades0.value3[6]*blades1.value1[2]
    +blades0.value4[1]*blades1.value2[1]
    -blades0.value3[8]*blades1.value1[3]
    +blades0.value4[2]*blades1.value2[3]
    -blades0.value4[4]*blades1.value2[5]
    -blades0.value5[0]*blades1.value3[2]
;
    blades.value3[4] =
    +blades0.value1[2]*blades1.value4[1]
    +blades0.value1[3]*blades1.value4[2]
    -blades0.value2[5]*blades1.value5[0]
    -blades0.value4[1]*blades1.value1[2]
    -blades0.value4[2]*blades1.value1[3]
    -blades0.value5[0]*blades1.value2[5]
;
    blades.value2[8] =
    +blades0.value1[0]*blades1.value3[5]
    +blades0.value1[1]*blades1.value3[6]
    -blades0.value2[0]*blades1.value4[1]
    -blades0.value1[3]*blades1.value3[9]
    +blades0.value2[3]*blades1.value4[3]
    +blades0.value2[4]*blades1.value4[4]
    +blades0.value3[1]*blades1.value5[0]
    +blades0.value3[5]*blades1.value1[0]
    +blades0.value3[6]*blades1.value1[1]
    -blades0.value4[1]*blades1.value2[0]
    -blades0.value3[9]*blades1.value1[3]
    +blades0.value4[3]*blades1.value2[3]
    +blades0.value4[4]*blades1.value2[4]
    +blades0.value5[0]*blades1.value3[1]
;
    blades.value3[5] =
    -blades0.value1[1]*blades1.value4[1]
    +blades0.value1[3]*blades1.value4[3]
    +blades0.value2[4]*blades1.value5[0]
    +blades0.value4[1]*blades1.value1[1]
    -blades0.value4[3]*blades1.value1[3]
    +blades0.value5[0]*blades1.value2[4]
;
    blades.value3[6] =
    +blades0.value1[0]*blades1.value4[1]
    +blades0.value1[3]*blades1.value4[4]
    -blades0.value2[3]*blades1.value5[0]
    -blades0.value4[1]*blades1.value1[0]
    -blades0.value4[4]*blades1.value1[3]
    -blades0.value5[0]*blades1.value2[3]
;
    blades.value4[1] =
    -blades0.value1[3]*blades1.value5[0]
    -blades0.value5[0]*blades1.value1[3]
;
    blades.value2[9] =
    +blades0.value1[0]*blades1.value3[7]
    +blades0.value1[1]*blades1.value3[8]
    -blades0.value2[0]*blades1.value4[2]
    +blades0.value1[2]*blades1.value3[9]
    -blades0.value2[1]*blades1.value4[3]
    -blades0.value2[2]*blades1.value4[4]
    -blades0.value3[0]*blades1.value5[0]
    +blades0.value3[7]*blades1.value1[0]
    +blades0.value3[8]*blades1.value1[1]
    -blades0.value4[2]*blades1.value2[0]
    +blades0.value3[9]*blades1.value1[2]
    -blades0.value4[3]*blades1.value2[1]
    -blades0.value4[4]*blades1.value2[2]
    -blades0.value5[0]*blades1.value3[0]
;
    blades.value3[7] =
    -blades0.value1[1]*blades1.value4[2]
    -blades0.value1[2]*blades1.value4[3]
    -blades0.value2[2]*blades1.value5[0]
    +blades0.value4[2]*blades1.value1[1]
    +blades0.value4[3]*blades1.value1[2]
    -blades0.value5[0]*blades1.value2[2]
;
    blades.value3[8] =
    +blades0.value1[0]*blades1.value4[2]
    -blades0.value1[2]*blades1.value4[4]
    +blades0.value2[1]*blades1.value5[0]
    -blades0.value4[2]*blades1.value1[0]
    +blades0.value4[4]*blades1.value1[2]
    +blades0.value5[0]*blades1.value2[1]
;
    blades.value4[2] =
    +blades0.value1[2]*blades1.value5[0]
    +blades0.value5[0]*blades1.value1[2]
;
    blades.value3[9] =
    +blades0.value1[0]*blades1.value4[3]
    +blades0.value1[1]*blades1.value4[4]
    -blades0.value2[0]*blades1.value5[0]
    -blades0.value4[3]*blades1.value1[0]
    -blades0.value4[4]*blades1.value1[1]
    -blades0.value5[0]*blades1.value2[0]
;
    blades.value4[3] =
    -blades0.value1[1]*blades1.value5[0]
    -blades0.value5[0]*blades1.value1[1]
;
    blades.value4[4] =
    +blades0.value1[0]*blades1.value5[0]
    +blades0.value5[0]*blades1.value1[0]
;
    return blades;
}


static gen1_BladesMultivector gen1_blades_outerproduct(gen1_BladesMultivector blades0, gen1_BladesMultivector blades1){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    blades.value0[0] =
    +blades0.value0[0]*blades1.value0[0]
;
    blades.value1[0] =
    +blades0.value0[0]*blades1.value1[0]
    +blades0.value1[0]*blades1.value0[0]
;
    blades.value1[1] =
    +blades0.value0[0]*blades1.value1[1]
    +blades0.value1[1]*blades1.value0[0]
;
    blades.value2[0] =
    +blades0.value0[0]*blades1.value2[0]
    +blades0.value1[0]*blades1.value1[1]
    -blades0.value1[1]*blades1.value1[0]
    +blades0.value2[0]*blades1.value0[0]
;
    blades.value1[2] =
    +blades0.value0[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value0[0]
;
    blades.value2[1] =
    +blades0.value0[0]*blades1.value2[1]
    +blades0.value1[0]*blades1.value1[2]
    -blades0.value1[2]*blades1.value1[0]
    +blades0.value2[1]*blades1.value0[0]
;
    blades.value2[2] =
    +blades0.value0[0]*blades1.value2[2]
    +blades0.value1[1]*blades1.value1[2]
    -blades0.value1[2]*blades1.value1[1]
    +blades0.value2[2]*blades1.value0[0]
;
    blades.value3[0] =
    +blades0.value0[0]*blades1.value3[0]
    +blades0.value1[0]*blades1.value2[2]
    -blades0.value1[1]*blades1.value2[1]
    +blades0.value2[0]*blades1.value1[2]
    +blades0.value1[2]*blades1.value2[0]
    -blades0.value2[1]*blades1.value1[1]
    +blades0.value2[2]*blades1.value1[0]
    +blades0.value3[0]*blades1.value0[0]
;
    blades.value1[3] =
    +blades0.value0[0]*blades1.value1[3]
    +blades0.value1[3]*blades1.value0[0]
;
    blades.value2[3] =
    +blades0.value0[0]*blades1.value2[3]
    +blades0.value1[0]*blades1.value1[3]
    -blades0.value1[3]*blades1.value1[0]
    +blades0.value2[3]*blades1.value0[0]
;
    blades.value2[4] =
    +blades0.value0[0]*blades1.value2[4]
    +blades0.value1[1]*blades1.value1[3]
    -blades0.value1[3]*blades1.value1[1]
    +blades0.value2[4]*blades1.value0[0]
;
    blades.value3[1] =
    +blades0.value0[0]*blades1.value3[1]
    +blades0.value1[0]*blades1.value2[4]
    -blades0.value1[1]*blades1.value2[3]
    +blades0.value2[0]*blades1.value1[3]
    +blades0.value1[3]*blades1.value2[0]
    -blades0.value2[3]*blades1.value1[1]
    +blades0.value2[4]*blades1.value1[0]
    +blades0.value3[1]*blades1.value0[0]
;
    blades.value2[5] =
    +blades0.value0[0]*blades1.value2[5]
    +blades0.value1[2]*blades1.value1[3]
    -blades0.value1[3]*blades1.value1[2]
    +blades0.value2[5]*blades1.value0[0]
;
    blades.value3[2] =
    +blades0.value0[0]*blades1.value3[2]
    +blades0.value1[0]*blades1.value2[5]
    -blades0.value1[2]*blades1.value2[3]
    +blades0.value2[1]*blades1.value1[3]
    +blades0.value1[3]*blades1.value2[1]
    -blades0.value2[3]*blades1.value1[2]
    +blades0.value2[5]*blades1.value1[0]
    +blades0.value3[2]*blades1.value0[0]
;
    blades.value3[3] =
    +blades0.value0[0]*blades1.value3[3]
    +blades0.value1[1]*blades1.value2[5]
    -blades0.value1[2]*blades1.value2[4]
    +blades0.value2[2]*blades1.value1[3]
    +blades0.value1[3]*blades1.value2[2]
    -blades0.value2[4]*blades1.value1[2]
    +blades0.value2[5]*blades1.value1[1]
    +blades0.value3[3]*blades1.value0[0]
;
    blades.value4[0] =
    +blades0.value0[0]*blades1.value4[0]
    +blades0.value1[0]*blades1.value3[3]
    -blades0.value1[1]*blades1.value3[2]
    +blades0.value2[0]*blades1.value2[5]
    +blades0.value1[2]*blades1.value3[1]
    -blades0.value2[1]*blades1.value2[4]
    +blades0.value2[2]*blades1.value2[3]
    +blades0.value3[0]*blades1.value1[3]
    -blades0.value1[3]*blades1.value3[0]
    +blades0.value2[3]*blades1.value2[2]
    -blades0.value2[4]*blades1.value2[1]
    -blades0.value3[1]*blades1.value1[2]
    +blades0.value2[5]*blades1.value2[0]
    +blades0.value3[2]*blades1.value1[1]
    -blades0.value3[3]*blades1.value1[0]
    +blades0.value4[0]*blades1.value0[0]
;
    blades.value1[4] =
    +blades0.value0[0]*blades1.value1[4]
    +blades0.value1[4]*blades1.value0[0]
;
    blades.value2[6] =
    +blades0.value0[0]*blades1.value2[6]
    +blades0.value1[0]*blades1.value1[4]
    -blades0.value1[4]*blades1.value1[0]
    +blades0.value2[6]*blades1.value0[0]
;
    blades.value2[7] =
    +blades0.value0[0]*blades1.value2[7]
    +blades0.value1[1]*blades1.value1[4]
    -blades0.value1[4]*blades1.value1[1]
    +blades0.value2[7]*blades1.value0[0]
;
    blades.value3[4] =
    +blades0.value0[0]*blades1.value3[4]
    +blades0.value1[0]*blades1.value2[7]
    -blades0.value1[1]*blades1.value2[6]
    +blades0.value2[0]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[0]
    -blades0.value2[6]*blades1.value1[1]
    +blades0.value2[7]*blades1.value1[0]
    +blades0.value3[4]*blades1.value0[0]
;
    blades.value2[8] =
    +blades0.value0[0]*blades1.value2[8]
    +blades0.value1[2]*blades1.value1[4]
    -blades0.value1[4]*blades1.value1[2]
    +blades0.value2[8]*blades1.value0[0]
;
    blades.value3[5] =
    +blades0.value0[0]*blades1.value3[5]
    +blades0.value1[0]*blades1.value2[8]
    -blades0.value1[2]*blades1.value2[6]
    +blades0.value2[1]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[1]
    -blades0.value2[6]*blades1.value1[2]
    +blades0.value2[8]*blades1.value1[0]
    +blades0.value3[5]*blades1.value0[0]
;
    blades.value3[6] =
    +blades0.value0[0]*blades1.value3[6]
    +blades0.value1[1]*blades1.value2[8]
    -blades0.value1[2]*blades1.value2[7]
    +blades0.value2[2]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[2]
    -blades0.value2[7]*blades1.value1[2]
    +blades0.value2[8]*blades1.value1[1]
    +blades0.value3[6]*blades1.value0[0]
;
    blades.value4[1] =
    +blades0.value0[0]*blades1.value4[1]
    +blades0.value1[0]*blades1.value3[6]
    -blades0.value1[1]*blades1.value3[5]
    +blades0.value2[0]*blades1.value2[8]
    +blades0.value1[2]*blades1.value3[4]
    -blades0.value2[1]*blades1.value2[7]
    +blades0.value2[2]*blades1.value2[6]
    +blades0.value3[0]*blades1.value1[4]
    -blades0.value1[4]*blades1.value3[0]
    +blades0.value2[6]*blades1.value2[2]
    -blades0.value2[7]*blades1.value2[1]
    -blades0.value3[4]*blades1.value1[2]
    +blades0.value2[8]*blades1.value2[0]
    +blades0.value3[5]*blades1.value1[1]
    -blades0.value3[6]*blades1.value1[0]
    +blades0.value4[1]*blades1.value0[0]
;
    blades.value2[9] =
    +blades0.value0[0]*blades1.value2[9]
    +blades0.value1[3]*blades1.value1[4]
    -blades0.value1[4]*blades1.value1[3]
    +blades0.value2[9]*blades1.value0[0]
;
    blades.value3[7] =
    +blades0.value0[0]*blades1.value3[7]
    +blades0.value1[0]*blades1.value2[9]
    -blades0.value1[3]*blades1.value2[6]
    +blades0.value2[3]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[3]
    -blades0.value2[6]*blades1.value1[3]
    +blades0.value2[9]*blades1.value1[0]
    +blades0.value3[7]*blades1.value0[0]
;
    blades.value3[8] =
    +blades0.value0[0]*blades1.value3[8]
    +blades0.value1[1]*blades1.value2[9]
    -blades0.value1[3]*blades1.value2[7]
    +blades0.value2[4]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[4]
    -blades0.value2[7]*blades1.value1[3]
    +blades0.value2[9]*blades1.value1[1]
    +blades0.value3[8]*blades1.value0[0]
;
    blades.value4[2] =
    +blades0.value0[0]*blades1.value4[2]
    +blades0.value1[0]*blades1.value3[8]
    -blades0.value1[1]*blades1.value3[7]
    +blades0.value2[0]*blades1.value2[9]
    +blades0.value1[3]*blades1.value3[4]
    -blades0.value2[3]*blades1.value2[7]
    +blades0.value2[4]*blades1.value2[6]
    +blades0.value3[1]*blades1.value1[4]
    -blades0.value1[4]*blades1.value3[1]
    +blades0.value2[6]*blades1.value2[4]
    -blades0.value2[7]*blades1.value2[3]
    -blades0.value3[4]*blades1.value1[3]
    +blades0.value2[9]*blades1.value2[0]
    +blades0.value3[7]*blades1.value1[1]
    -blades0.value3[8]*blades1.value1[0]
    +blades0.value4[2]*blades1.value0[0]
;
    blades.value3[9] =
    +blades0.value0[0]*blades1.value3[9]
    +blades0.value1[2]*blades1.value2[9]
    -blades0.value1[3]*blades1.value2[8]
    +blades0.value2[5]*blades1.value1[4]
    +blades0.value1[4]*blades1.value2[5]
    -blades0.value2[8]*blades1.value1[3]
    +blades0.value2[9]*blades1.value1[2]
    +blades0.value3[9]*blades1.value0[0]
;
    blades.value4[3] =
    +blades0.value0[0]*blades1.value4[3]
    +blades0.value1[0]*blades1.value3[9]
    -blades0.value1[2]*blades1.value3[7]
    +blades0.value2[1]*blades1.value2[9]
    +blades0.value1[3]*blades1.value3[5]
    -blades0.value2[3]*blades1.value2[8]
    +blades0.value2[5]*blades1.value2[6]
    +blades0.value3[2]*blades1.value1[4]
    -blades0.value1[4]*blades1.value3[2]
    +blades0.value2[6]*blades1.value2[5]
    -blades0.value2[8]*blades1.value2[3]
    -blades0.value3[5]*blades1.value1[3]
    +blades0.value2[9]*blades1.value2[1]
    +blades0.value3[7]*blades1.value1[2]
    -blades0.value3[9]*blades1.value1[0]
    +blades0.value4[3]*blades1.value0[0]
;
    blades.value4[4] =
    +blades0.value0[0]*blades1.value4[4]
    +blades0.value1[1]*blades1.value3[9]
    -blades0.value1[2]*blades1.value3[8]
    +blades0.value2[2]*blades1.value2[9]
    +blades0.value1[3]*blades1.value3[6]
    -blades0.value2[4]*blades1.value2[8]
    +blades0.value2[5]*blades1.value2[7]
    +blades0.value3[3]*blades1.value1[4]
    -blades0.value1[4]*blades1.value3[3]
    +blades0.value2[7]*blades1.value2[5]
    -blades0.value2[8]*blades1.value2[4]
    -blades0.value3[6]*blades1.value1[3]
    +blades0.value2[9]*blades1.value2[2]
    +blades0.value3[8]*blades1.value1[2]
    -blades0.value3[9]*blades1.value1[1]
    +blades0.value4[4]*blades1.value0[0]
;
    blades.value5[0] =
    +blades0.value0[0]*blades1.value5[0]
    +blades0.value1[0]*blades1.value4[4]
    -blades0.value1[1]*blades1.value4[3]
    +blades0.value2[0]*blades1.value3[9]
    +blades0.value1[2]*blades1.value4[2]
    -blades0.value2[1]*blades1.value3[8]
    +blades0.value2[2]*blades1.value3[7]
    +blades0.value3[0]*blades1.value2[9]
    -blades0.value1[3]*blades1.value4[1]
    +blades0.value2[3]*blades1.value3[6]
    -blades0.value2[4]*blades1.value3[5]
    -blades0.value3[1]*blades1.value2[8]
    +blades0.value2[5]*blades1.value3[4]
    +blades0.value3[2]*blades1.value2[7]
    -blades0.value3[3]*blades1.value2[6]
    +blades0.value4[0]*blades1.value1[4]
    +blades0.value1[4]*blades1.value4[0]
    -blades0.value2[6]*blades1.value3[3]
    +blades0.value2[7]*blades1.value3[2]
    +blades0.value3[4]*blades1.value2[5]
    -blades0.value2[8]*blades1.value3[1]
    -blades0.value3[5]*blades1.value2[4]
    +blades0.value3[6]*blades1.value2[3]
    -blades0.value4[1]*blades1.value1[3]
    +blades0.value2[9]*blades1.value3[0]
    +blades0.value3[7]*blades1.value2[2]
    -blades0.value3[8]*blades1.value2[1]
    +blades0.value4[2]*blades1.value1[2]
    +blades0.value3[9]*blades1.value2[0]
    -blades0.value4[3]*blades1.value1[1]
    +blades0.value4[4]*blades1.value1[0]
    +blades0.value5[0]*blades1.value0[0]
;
    return blades;
}


static gen1_BladesMultivector gen1_blades_atomicadd(gen1_BladesMultivector *blades_array, Py_ssize_t size){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    for(Py_ssize_t i = 0; i < size; i++){
       blades.value0[0] += blades_array[i].value0[0];
       blades.value1[0] += blades_array[i].value1[0];
       blades.value1[1] += blades_array[i].value1[1];
       blades.value1[2] += blades_array[i].value1[2];
       blades.value1[3] += blades_array[i].value1[3];
       blades.value1[4] += blades_array[i].value1[4];
       blades.value2[0] += blades_array[i].value2[0];
       blades.value2[1] += blades_array[i].value2[1];
       blades.value2[2] += blades_array[i].value2[2];
       blades.value2[3] += blades_array[i].value2[3];
       blades.value2[4] += blades_array[i].value2[4];
       blades.value2[5] += blades_array[i].value2[5];
       blades.value2[6] += blades_array[i].value2[6];
       blades.value2[7] += blades_array[i].value2[7];
       blades.value2[8] += blades_array[i].value2[8];
       blades.value2[9] += blades_array[i].value2[9];
       blades.value3[0] += blades_array[i].value3[0];
       blades.value3[1] += blades_array[i].value3[1];
       blades.value3[2] += blades_array[i].value3[2];
       blades.value3[3] += blades_array[i].value3[3];
       blades.value3[4] += blades_array[i].value3[4];
       blades.value3[5] += blades_array[i].value3[5];
       blades.value3[6] += blades_array[i].value3[6];
       blades.value3[7] += blades_array[i].value3[7];
       blades.value3[8] += blades_array[i].value3[8];
       blades.value3[9] += blades_array[i].value3[9];
       blades.value4[0] += blades_array[i].value4[0];
       blades.value4[1] += blades_array[i].value4[1];
       blades.value4[2] += blades_array[i].value4[2];
       blades.value4[3] += blades_array[i].value4[3];
       blades.value4[4] += blades_array[i].value4[4];
       blades.value5[0] += blades_array[i].value5[0];
    }
    return blades;
}

static gen1_BladesMultivector gen1_blades_add(gen1_BladesMultivector blades0, gen1_BladesMultivector blades1, int sign){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    if(sign == -1){
        blades.value0[0] = blades0.value0[0] - blades1.value0[0];
        blades.value1[0] = blades0.value1[0] - blades1.value1[0];
        blades.value1[1] = blades0.value1[1] - blades1.value1[1];
        blades.value1[2] = blades0.value1[2] - blades1.value1[2];
        blades.value1[3] = blades0.value1[3] - blades1.value1[3];
        blades.value1[4] = blades0.value1[4] - blades1.value1[4];
        blades.value2[0] = blades0.value2[0] - blades1.value2[0];
        blades.value2[1] = blades0.value2[1] - blades1.value2[1];
        blades.value2[2] = blades0.value2[2] - blades1.value2[2];
        blades.value2[3] = blades0.value2[3] - blades1.value2[3];
        blades.value2[4] = blades0.value2[4] - blades1.value2[4];
        blades.value2[5] = blades0.value2[5] - blades1.value2[5];
        blades.value2[6] = blades0.value2[6] - blades1.value2[6];
        blades.value2[7] = blades0.value2[7] - blades1.value2[7];
        blades.value2[8] = blades0.value2[8] - blades1.value2[8];
        blades.value2[9] = blades0.value2[9] - blades1.value2[9];
        blades.value3[0] = blades0.value3[0] - blades1.value3[0];
        blades.value3[1] = blades0.value3[1] - blades1.value3[1];
        blades.value3[2] = blades0.value3[2] - blades1.value3[2];
        blades.value3[3] = blades0.value3[3] - blades1.value3[3];
        blades.value3[4] = blades0.value3[4] - blades1.value3[4];
        blades.value3[5] = blades0.value3[5] - blades1.value3[5];
        blades.value3[6] = blades0.value3[6] - blades1.value3[6];
        blades.value3[7] = blades0.value3[7] - blades1.value3[7];
        blades.value3[8] = blades0.value3[8] - blades1.value3[8];
        blades.value3[9] = blades0.value3[9] - blades1.value3[9];
        blades.value4[0] = blades0.value4[0] - blades1.value4[0];
        blades.value4[1] = blades0.value4[1] - blades1.value4[1];
        blades.value4[2] = blades0.value4[2] - blades1.value4[2];
        blades.value4[3] = blades0.value4[3] - blades1.value4[3];
        blades.value4[4] = blades0.value4[4] - blades1.value4[4];
        blades.value5[0] = blades0.value5[0] - blades1.value5[0];
    }else if(sign == 1){
        blades.value0[0] = blades0.value0[0] + blades1.value0[0];
        blades.value1[0] = blades0.value1[0] + blades1.value1[0];
        blades.value1[1] = blades0.value1[1] + blades1.value1[1];
        blades.value1[2] = blades0.value1[2] + blades1.value1[2];
        blades.value1[3] = blades0.value1[3] + blades1.value1[3];
        blades.value1[4] = blades0.value1[4] + blades1.value1[4];
        blades.value2[0] = blades0.value2[0] + blades1.value2[0];
        blades.value2[1] = blades0.value2[1] + blades1.value2[1];
        blades.value2[2] = blades0.value2[2] + blades1.value2[2];
        blades.value2[3] = blades0.value2[3] + blades1.value2[3];
        blades.value2[4] = blades0.value2[4] + blades1.value2[4];
        blades.value2[5] = blades0.value2[5] + blades1.value2[5];
        blades.value2[6] = blades0.value2[6] + blades1.value2[6];
        blades.value2[7] = blades0.value2[7] + blades1.value2[7];
        blades.value2[8] = blades0.value2[8] + blades1.value2[8];
        blades.value2[9] = blades0.value2[9] + blades1.value2[9];
        blades.value3[0] = blades0.value3[0] + blades1.value3[0];
        blades.value3[1] = blades0.value3[1] + blades1.value3[1];
        blades.value3[2] = blades0.value3[2] + blades1.value3[2];
        blades.value3[3] = blades0.value3[3] + blades1.value3[3];
        blades.value3[4] = blades0.value3[4] + blades1.value3[4];
        blades.value3[5] = blades0.value3[5] + blades1.value3[5];
        blades.value3[6] = blades0.value3[6] + blades1.value3[6];
        blades.value3[7] = blades0.value3[7] + blades1.value3[7];
        blades.value3[8] = blades0.value3[8] + blades1.value3[8];
        blades.value3[9] = blades0.value3[9] + blades1.value3[9];
        blades.value4[0] = blades0.value4[0] + blades1.value4[0];
        blades.value4[1] = blades0.value4[1] + blades1.value4[1];
        blades.value4[2] = blades0.value4[2] + blades1.value4[2];
        blades.value4[3] = blades0.value4[3] + blades1.value4[3];
        blades.value4[4] = blades0.value4[4] + blades1.value4[4];
        blades.value5[0] = blades0.value5[0] + blades1.value5[0];
    }else{
        blades.value0[0] = blades0.value0[0] + sign*blades1.value0[0];
        blades.value1[0] = blades0.value1[0] + sign*blades1.value1[0];
        blades.value1[1] = blades0.value1[1] + sign*blades1.value1[1];
        blades.value1[2] = blades0.value1[2] + sign*blades1.value1[2];
        blades.value1[3] = blades0.value1[3] + sign*blades1.value1[3];
        blades.value1[4] = blades0.value1[4] + sign*blades1.value1[4];
        blades.value2[0] = blades0.value2[0] + sign*blades1.value2[0];
        blades.value2[1] = blades0.value2[1] + sign*blades1.value2[1];
        blades.value2[2] = blades0.value2[2] + sign*blades1.value2[2];
        blades.value2[3] = blades0.value2[3] + sign*blades1.value2[3];
        blades.value2[4] = blades0.value2[4] + sign*blades1.value2[4];
        blades.value2[5] = blades0.value2[5] + sign*blades1.value2[5];
        blades.value2[6] = blades0.value2[6] + sign*blades1.value2[6];
        blades.value2[7] = blades0.value2[7] + sign*blades1.value2[7];
        blades.value2[8] = blades0.value2[8] + sign*blades1.value2[8];
        blades.value2[9] = blades0.value2[9] + sign*blades1.value2[9];
        blades.value3[0] = blades0.value3[0] + sign*blades1.value3[0];
        blades.value3[1] = blades0.value3[1] + sign*blades1.value3[1];
        blades.value3[2] = blades0.value3[2] + sign*blades1.value3[2];
        blades.value3[3] = blades0.value3[3] + sign*blades1.value3[3];
        blades.value3[4] = blades0.value3[4] + sign*blades1.value3[4];
        blades.value3[5] = blades0.value3[5] + sign*blades1.value3[5];
        blades.value3[6] = blades0.value3[6] + sign*blades1.value3[6];
        blades.value3[7] = blades0.value3[7] + sign*blades1.value3[7];
        blades.value3[8] = blades0.value3[8] + sign*blades1.value3[8];
        blades.value3[9] = blades0.value3[9] + sign*blades1.value3[9];
        blades.value4[0] = blades0.value4[0] + sign*blades1.value4[0];
        blades.value4[1] = blades0.value4[1] + sign*blades1.value4[1];
        blades.value4[2] = blades0.value4[2] + sign*blades1.value4[2];
        blades.value4[3] = blades0.value4[3] + sign*blades1.value4[3];
        blades.value4[4] = blades0.value4[4] + sign*blades1.value4[4];
        blades.value5[0] = blades0.value5[0] + sign*blades1.value5[0];
    }
    return blades;
}


static gen1_BladesMultivector gen1_blades_scalaradd(gen1_BladesMultivector blades0, ga_float value, int sign){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    if(sign == -1){
        blades.value0[0] = -blades0.value0[0];
        blades.value1[0] = -blades0.value1[0];
        blades.value1[1] = -blades0.value1[1];
        blades.value1[2] = -blades0.value1[2];
        blades.value1[3] = -blades0.value1[3];
        blades.value1[4] = -blades0.value1[4];
        blades.value2[0] = -blades0.value2[0];
        blades.value2[1] = -blades0.value2[1];
        blades.value2[2] = -blades0.value2[2];
        blades.value2[3] = -blades0.value2[3];
        blades.value2[4] = -blades0.value2[4];
        blades.value2[5] = -blades0.value2[5];
        blades.value2[6] = -blades0.value2[6];
        blades.value2[7] = -blades0.value2[7];
        blades.value2[8] = -blades0.value2[8];
        blades.value2[9] = -blades0.value2[9];
        blades.value3[0] = -blades0.value3[0];
        blades.value3[1] = -blades0.value3[1];
        blades.value3[2] = -blades0.value3[2];
        blades.value3[3] = -blades0.value3[3];
        blades.value3[4] = -blades0.value3[4];
        blades.value3[5] = -blades0.value3[5];
        blades.value3[6] = -blades0.value3[6];
        blades.value3[7] = -blades0.value3[7];
        blades.value3[8] = -blades0.value3[8];
        blades.value3[9] = -blades0.value3[9];
        blades.value4[0] = -blades0.value4[0];
        blades.value4[1] = -blades0.value4[1];
        blades.value4[2] = -blades0.value4[2];
        blades.value4[3] = -blades0.value4[3];
        blades.value4[4] = -blades0.value4[4];
        blades.value5[0] = -blades0.value5[0];
    }else if(sign == 1){
        blades.value0[0] = blades0.value0[0];
        blades.value1[0] = blades0.value1[0];
        blades.value1[1] = blades0.value1[1];
        blades.value1[2] = blades0.value1[2];
        blades.value1[3] = blades0.value1[3];
        blades.value1[4] = blades0.value1[4];
        blades.value2[0] = blades0.value2[0];
        blades.value2[1] = blades0.value2[1];
        blades.value2[2] = blades0.value2[2];
        blades.value2[3] = blades0.value2[3];
        blades.value2[4] = blades0.value2[4];
        blades.value2[5] = blades0.value2[5];
        blades.value2[6] = blades0.value2[6];
        blades.value2[7] = blades0.value2[7];
        blades.value2[8] = blades0.value2[8];
        blades.value2[9] = blades0.value2[9];
        blades.value3[0] = blades0.value3[0];
        blades.value3[1] = blades0.value3[1];
        blades.value3[2] = blades0.value3[2];
        blades.value3[3] = blades0.value3[3];
        blades.value3[4] = blades0.value3[4];
        blades.value3[5] = blades0.value3[5];
        blades.value3[6] = blades0.value3[6];
        blades.value3[7] = blades0.value3[7];
        blades.value3[8] = blades0.value3[8];
        blades.value3[9] = blades0.value3[9];
        blades.value4[0] = blades0.value4[0];
        blades.value4[1] = blades0.value4[1];
        blades.value4[2] = blades0.value4[2];
        blades.value4[3] = blades0.value4[3];
        blades.value4[4] = blades0.value4[4];
        blades.value5[0] = blades0.value5[0];
    }else{
        blades.value0[0] = sign*blades0.value0[0];
        blades.value1[0] = sign*blades0.value1[0];
        blades.value1[1] = sign*blades0.value1[1];
        blades.value1[2] = sign*blades0.value1[2];
        blades.value1[3] = sign*blades0.value1[3];
        blades.value1[4] = sign*blades0.value1[4];
        blades.value2[0] = sign*blades0.value2[0];
        blades.value2[1] = sign*blades0.value2[1];
        blades.value2[2] = sign*blades0.value2[2];
        blades.value2[3] = sign*blades0.value2[3];
        blades.value2[4] = sign*blades0.value2[4];
        blades.value2[5] = sign*blades0.value2[5];
        blades.value2[6] = sign*blades0.value2[6];
        blades.value2[7] = sign*blades0.value2[7];
        blades.value2[8] = sign*blades0.value2[8];
        blades.value2[9] = sign*blades0.value2[9];
        blades.value3[0] = sign*blades0.value3[0];
        blades.value3[1] = sign*blades0.value3[1];
        blades.value3[2] = sign*blades0.value3[2];
        blades.value3[3] = sign*blades0.value3[3];
        blades.value3[4] = sign*blades0.value3[4];
        blades.value3[5] = sign*blades0.value3[5];
        blades.value3[6] = sign*blades0.value3[6];
        blades.value3[7] = sign*blades0.value3[7];
        blades.value3[8] = sign*blades0.value3[8];
        blades.value3[9] = sign*blades0.value3[9];
        blades.value4[0] = sign*blades0.value4[0];
        blades.value4[1] = sign*blades0.value4[1];
        blades.value4[2] = sign*blades0.value4[2];
        blades.value4[3] = sign*blades0.value4[3];
        blades.value4[4] = sign*blades0.value4[4];
        blades.value5[0] = sign*blades0.value5[0];
    }
    blades.value0[0] += value;
    return blades;
}


static gen1_BladesMultivector gen1_blades_scalarproduct(gen1_BladesMultivector blades0, ga_float value){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    blades.value0[0] = value*blades0.value0[0];
    blades.value1[0] = value*blades0.value1[0];
    blades.value1[1] = value*blades0.value1[1];
    blades.value1[2] = value*blades0.value1[2];
    blades.value1[3] = value*blades0.value1[3];
    blades.value1[4] = value*blades0.value1[4];
    blades.value2[0] = value*blades0.value2[0];
    blades.value2[1] = value*blades0.value2[1];
    blades.value2[2] = value*blades0.value2[2];
    blades.value2[3] = value*blades0.value2[3];
    blades.value2[4] = value*blades0.value2[4];
    blades.value2[5] = value*blades0.value2[5];
    blades.value2[6] = value*blades0.value2[6];
    blades.value2[7] = value*blades0.value2[7];
    blades.value2[8] = value*blades0.value2[8];
    blades.value2[9] = value*blades0.value2[9];
    blades.value3[0] = value*blades0.value3[0];
    blades.value3[1] = value*blades0.value3[1];
    blades.value3[2] = value*blades0.value3[2];
    blades.value3[3] = value*blades0.value3[3];
    blades.value3[4] = value*blades0.value3[4];
    blades.value3[5] = value*blades0.value3[5];
    blades.value3[6] = value*blades0.value3[6];
    blades.value3[7] = value*blades0.value3[7];
    blades.value3[8] = value*blades0.value3[8];
    blades.value3[9] = value*blades0.value3[9];
    blades.value4[0] = value*blades0.value4[0];
    blades.value4[1] = value*blades0.value4[1];
    blades.value4[2] = value*blades0.value4[2];
    blades.value4[3] = value*blades0.value4[3];
    blades.value4[4] = value*blades0.value4[4];
    blades.value5[0] = value*blades0.value5[0];
    return blades;
}

static gen1_BladesMultivector gen1_blades_reverse(gen1_BladesMultivector blades0){
    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};

    blades.value0[0] = blades0.value0[0];
    blades.value1[0] = blades0.value1[0];
    blades.value1[1] = blades0.value1[1];
    blades.value1[2] = blades0.value1[2];
    blades.value1[3] = blades0.value1[3];
    blades.value1[4] = blades0.value1[4];
    blades.value2[0] = -blades0.value2[0];
    blades.value2[1] = -blades0.value2[1];
    blades.value2[2] = -blades0.value2[2];
    blades.value2[3] = -blades0.value2[3];
    blades.value2[4] = -blades0.value2[4];
    blades.value2[5] = -blades0.value2[5];
    blades.value2[6] = -blades0.value2[6];
    blades.value2[7] = -blades0.value2[7];
    blades.value2[8] = -blades0.value2[8];
    blades.value2[9] = -blades0.value2[9];
    blades.value3[0] = -blades0.value3[0];
    blades.value3[1] = -blades0.value3[1];
    blades.value3[2] = -blades0.value3[2];
    blades.value3[3] = -blades0.value3[3];
    blades.value3[4] = -blades0.value3[4];
    blades.value3[5] = -blades0.value3[5];
    blades.value3[6] = -blades0.value3[6];
    blades.value3[7] = -blades0.value3[7];
    blades.value3[8] = -blades0.value3[8];
    blades.value3[9] = -blades0.value3[9];
    blades.value4[0] = blades0.value4[0];
    blades.value4[1] = blades0.value4[1];
    blades.value4[2] = blades0.value4[2];
    blades.value4[3] = blades0.value4[3];
    blades.value4[4] = blades0.value4[4];
    blades.value5[0] = blades0.value5[0];
    return blades;
}


static PyMultivectorObject *binary_dense1_product(PyMultivectorObject *data0, PyMultivectorObject *data1,ProductType ptype){
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)data0->data;
    gen1_DenseMultivector *pdense1 = (gen1_DenseMultivector*)data1->data;
    gen1_DenseMultivector *pdense  = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pdense0 || !pdense1 || !pdense || !out){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = gen1_dense_geometricproduct(*pdense0,*pdense1);
            break;
        case ProductType_inner:
            *pdense = gen1_dense_innerproduct(*pdense0,*pdense1);
            break;
        case ProductType_outer:
            *pdense = gen1_dense_outerproduct(*pdense0,*pdense1);
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades1_product(PyMultivectorObject *data0, PyMultivectorObject *data1,ProductType ptype){
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)data0->data;
    gen1_BladesMultivector *pblades1 = (gen1_BladesMultivector*)data1->data;
    gen1_BladesMultivector *pblades  = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pblades0 || !pblades1 || !pblades || !out){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = gen1_blades_geometricproduct(*pblades0,*pblades1);
            break;
        case ProductType_inner:
            *pblades = gen1_blades_innerproduct(*pblades0,*pblades1);
            break;
        case ProductType_outer:
            *pblades = gen1_blades_outerproduct(*pblades0,*pblades1);
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *ternary_dense1_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2,ProductType ptype){
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)data0->data;
    gen1_DenseMultivector *pdense1 = (gen1_DenseMultivector*)data1->data;
    gen1_DenseMultivector *pdense2 = (gen1_DenseMultivector*)data2->data;
    gen1_DenseMultivector *pdense  = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pdense0 || !pdense1 || !pdense2 || !pdense || !out){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = gen1_dense_geometricproduct(*pdense0,*pdense1);
            *pdense = gen1_dense_geometricproduct(*pdense,*pdense2);
            break;
        case ProductType_inner:
            *pdense = gen1_dense_innerproduct(*pdense0,*pdense1);
            *pdense = gen1_dense_innerproduct(*pdense,*pdense2);
            break;
        case ProductType_outer:
            *pdense = gen1_dense_outerproduct(*pdense0,*pdense1);
            *pdense = gen1_dense_outerproduct(*pdense,*pdense2);
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *ternary_blades1_product(PyMultivectorObject *data0, PyMultivectorObject *data1, PyMultivectorObject *data2,ProductType ptype){
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)data0->data;
    gen1_BladesMultivector *pblades1 = (gen1_BladesMultivector*)data1->data;
    gen1_BladesMultivector *pblades2 = (gen1_BladesMultivector*)data2->data;
    gen1_BladesMultivector *pblades  = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pblades0 || !pblades1 || !pblades2 || !pblades || !out){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = gen1_blades_geometricproduct(*pblades0,*pblades1);
            *pblades = gen1_blades_geometricproduct(*pblades,*pblades2);
            break;
        case ProductType_inner:
            *pblades = gen1_blades_innerproduct(*pblades0,*pblades1);
            *pblades = gen1_blades_innerproduct(*pblades,*pblades2);
            break;
        case ProductType_outer:
            *pblades = gen1_blades_outerproduct(*pblades0,*pblades1);
            *pblades = gen1_blades_outerproduct(*pblades,*pblades2);
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }

    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *unary_dense1_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;
    gen1_DenseMultivector dense = {{0}};
    gen1_DenseMultivector *pdense;
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)self->data;

    for(Py_ssize_t i = 0; i < size; i++){
        gen1densegradeprojectfunc gradeproject =
                    gen1denseproject.gradeproject[grades[i]];
        if(gradeproject != NULL)
            gradeproject(pdense0,&dense);
        else
            return NULL; // raise not implemented error
    }
    // should allocate all the necessary memory and also set the reference count to 1
    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    *pdense = dense;
    out->data = (void*)pdense;
    PyMem_RawFree(grades);
    return out;
}

static PyMultivectorObject *unary_blades1_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;

    gen1_BladesMultivector blades = {{0},{0},{0},{0},{0},{0},};
    gen1_BladesMultivector *pblades;
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)self->data;

    for(Py_ssize_t i = 0; i < size; i++){
        gen1bladesgradeprojectfunc gradeproject =
                    gen1bladesproject.gradeproject[grades[i]];
        if(gradeproject != NULL)
            gradeproject(pblades0,&blades);
        else
            return NULL; // raise not implemented error
    }
    // should allocate all the necessary memory and also set the reference count to 1
    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    *pblades = blades;
    out->data = (void*)pblades;
    return out;
}

static PyMultivectorObject* atomic_dense1_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorObject *out = new_multivector(data,-1);
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)PyMem_RawMalloc(size*sizeof(gen1_DenseMultivector));
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense0);
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    for(Py_ssize_t i = 0; i < size; i++)
        pdense0[i] = *((gen1_DenseMultivector*)data[i].data);

    *pdense = gen1_dense_atomicadd(pdense0,size);
    out->data = (void*)pdense;
    return out;
}
static PyMultivectorObject* atomic_blades1_add(PyMultivectorObject *data, Py_ssize_t size){
    PyMultivectorObject *out = new_multivector(data,-1);
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)PyMem_RawMalloc(size*sizeof(gen1_BladesMultivector));
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades0);
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    for(Py_ssize_t i = 0; i < size; i++)
        pblades0[i] = *((gen1_BladesMultivector*)data[i].data);

    *pblades = gen1_blades_atomicadd(pblades0,size);
    out->data = (void*)pblades;
    return out;
}

static PyMultivectorObject* binary_dense1_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    PyMultivectorObject *out = new_multivector(data0,-1);
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)data0->data;
    gen1_DenseMultivector *pdense1 = (gen1_DenseMultivector*)data1->data;
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    if(!out || !pdense0 || !pdense1 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen1_dense_add(*pdense0,*pdense1,sign);
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject* binary_blades1_add(PyMultivectorObject *data0, PyMultivectorObject *data1, int sign){
    PyMultivectorObject *out = new_multivector(data0,-1);
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)data0->data;
    gen1_BladesMultivector *pblades1 = (gen1_BladesMultivector*)data1->data;
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    if(!out || !pblades0 || !pblades1 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen1_blades_add(*pblades0,*pblades1,sign);
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}


static PyMultivectorObject* atomic_dense1_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    if(size < 2) return NULL;
    PyMultivectorObject *out = new_multivector(data,-1);
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    gen1_DenseMultivector dense;
    if(!out  || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    switch(ptype){
        case ProductType_geometric:
            dense = gen1_dense_geometricproduct(
                      *((gen1_DenseMultivector*)data[0].data),
                      *((gen1_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen1_dense_geometricproduct(
                          dense,
                          *((gen1_DenseMultivector*)data[i].data));
            }
            break;
        case ProductType_inner:
            dense = gen1_dense_innerproduct(
                      *((gen1_DenseMultivector*)data[0].data),
                      *((gen1_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen1_dense_innerproduct(
                          dense,
                          *((gen1_DenseMultivector*)data[i].data));
            }
            break;
        case ProductType_outer:
            dense = gen1_dense_outerproduct(
                      *((gen1_DenseMultivector*)data[0].data),
                      *((gen1_DenseMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                dense = gen1_dense_outerproduct(
                          dense,
                          *((gen1_DenseMultivector*)data[i].data));
            }
            break;
        default:
            PyMem_RawFree(pdense);
            free_multivector(out);
            return NULL;
    }
    *pdense = dense;
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject* atomic_blades1_product(PyMultivectorObject *data, Py_ssize_t size, ProductType ptype){
    if(size < 2) return NULL;
    PyMultivectorObject *out = new_multivector(data,-1);
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    gen1_BladesMultivector blades;
    if(!out  || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    switch(ptype){
        case ProductType_geometric:
            blades = gen1_blades_geometricproduct(
                      *((gen1_BladesMultivector*)data[0].data),
                      *((gen1_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen1_blades_geometricproduct(
                          blades,
                          *((gen1_BladesMultivector*)data[i].data));
            }
            break;
        case ProductType_inner:
            blades = gen1_blades_innerproduct(
                      *((gen1_BladesMultivector*)data[0].data),
                      *((gen1_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen1_blades_innerproduct(
                          blades,
                          *((gen1_BladesMultivector*)data[i].data));
            }
            break;
        case ProductType_outer:
            blades = gen1_blades_outerproduct(
                      *((gen1_BladesMultivector*)data[0].data),
                      *((gen1_BladesMultivector*)data[1].data));
            for(Py_ssize_t i = 2; i < size; i++){
                blades = gen1_blades_outerproduct(
                          blades,
                          *((gen1_BladesMultivector*)data[i].data));
            }
            break;
        default:
            PyMem_RawFree(pblades);
            free_multivector(out);
            return NULL;
    }
    *pblades = blades;
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_dense1_scalarproduct(PyMultivectorObject *self, ga_float value){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)self->data;
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen1_dense_scalarproduct(*pdense0,value); // multiply by a scalar
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades1_scalarproduct(PyMultivectorObject *self, ga_float value){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)self->data;
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen1_blades_scalarproduct(*pblades0,value); // multiply by a scalar
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *binary_dense1_scalaradd(PyMultivectorObject *self, ga_float value, int sign){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)self->data;
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen1_dense_scalaradd(*pdense0,value,sign); // add a scalar
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *binary_blades1_scalaradd(PyMultivectorObject *self, ga_float value, int sign){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)self->data;
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen1_blades_scalaradd(*pblades0,value,sign); // add a scalar
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}

static PyMultivectorObject *unary_dense1_reverse(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_DenseMultivector *pdense0 = (gen1_DenseMultivector*)self->data;
    gen1_DenseMultivector *pdense = (gen1_DenseMultivector*)PyMem_RawMalloc(sizeof(gen1_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen1_dense_reverse(*pdense0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *unary_blades1_reverse(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen1_BladesMultivector *pblades0 = (gen1_BladesMultivector*)self->data;
    gen1_BladesMultivector *pblades = (gen1_BladesMultivector*)PyMem_RawMalloc(sizeof(gen1_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen1_blades_reverse(*pblades0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}



static const PyMultivectorMath_Funcs dense0_math_funcs = {
    .atomic_add = (gaatomicfunc)atomic_dense0_add,
    .atomic_product = (gaatomicprodfunc) atomic_dense0_product,
    .add = (gaaddfunc) binary_dense0_add,
    .product = (gaprodfunc) binary_dense0_product,
    .grade_project = (gaunarygradefunc) unary_dense0_gradeproject,
    .scalar_product = (gascalarfunc) binary_dense0_scalarproduct,
    .scalar_add = (gascalaraddfunc) binary_dense0_scalaradd,
    .reverse = (gaunaryfunc) unary_dense0_reverse,
    .ternary_product = (gaternaryprodfunc) ternary_dense0_product,
};

static const PyMultivectorMath_Funcs blades0_math_funcs = {
    .atomic_add = (gaatomicfunc)atomic_blades0_add,
    .atomic_product = (gaatomicprodfunc) atomic_blades0_product,
    .add = (gaaddfunc) binary_blades0_add,
    .product = (gaprodfunc) binary_blades0_product,
    .grade_project = (gaunarygradefunc) unary_blades0_gradeproject,
    .scalar_product = (gascalarfunc) binary_blades0_scalarproduct,
    .scalar_add = (gascalaraddfunc) binary_blades0_scalaradd,
    .reverse = (gaunaryfunc) unary_blades0_reverse,
    .ternary_product = (gaternaryprodfunc) ternary_blades0_product,
};


static const PyMultivectorData_Funcs dense0_data_funcs = {
  .iter_next = (gaiternextfunc) dense0_iternext,
  .iter_init = (gaiterinitfunc) dense0_iterinit,
  .init = (gainitfunc) dense0_init,
};

static const PyMultivectorData_Funcs blades0_data_funcs = {
  .iter_next = (gaiternextfunc) blades0_iternext,
  .iter_init = (gaiterinitfunc) blades0_iterinit,
  .init = (gainitfunc) blades0_init,
};


static const PyMultivectorSubType dense0_subtype = {
    .math_funcs = dense0_math_funcs,
    .data_funcs = dense0_data_funcs,
    .name = "3DVGA",
    .type_name = "dense",
    .generated = 1,
    .metric = {1,1,1,},
    .metric_size = 3,
    .ntype = 3,
};

static const PyMultivectorSubType blades0_subtype = {
    .math_funcs = blades0_math_funcs,
    .data_funcs = blades0_data_funcs,
    .name = "3DVGA",
    .type_name = "blades",
    .generated = 1,
    .metric = {1,1,1,},
    .metric_size = 3,
    .ntype = 5,
};


static const PyMultivectorMath_Funcs dense1_math_funcs = {
    .atomic_add = (gaatomicfunc)atomic_dense1_add,
    .atomic_product = (gaatomicprodfunc) atomic_dense1_product,
    .add = (gaaddfunc) binary_dense1_add,
    .product = (gaprodfunc) binary_dense1_product,
    .grade_project = (gaunarygradefunc) unary_dense1_gradeproject,
    .scalar_product = (gascalarfunc) binary_dense1_scalarproduct,
    .scalar_add = (gascalaraddfunc) binary_dense1_scalaradd,
    .reverse = (gaunaryfunc) unary_dense1_reverse,
    .ternary_product = (gaternaryprodfunc) ternary_dense1_product,
};

static const PyMultivectorMath_Funcs blades1_math_funcs = {
    .atomic_add = (gaatomicfunc)atomic_blades1_add,
    .atomic_product = (gaatomicprodfunc) atomic_blades1_product,
    .add = (gaaddfunc) binary_blades1_add,
    .product = (gaprodfunc) binary_blades1_product,
    .grade_project = (gaunarygradefunc) unary_blades1_gradeproject,
    .scalar_product = (gascalarfunc) binary_blades1_scalarproduct,
    .scalar_add = (gascalaraddfunc) binary_blades1_scalaradd,
    .reverse = (gaunaryfunc) unary_blades1_reverse,
    .ternary_product = (gaternaryprodfunc) ternary_blades1_product,
};


static const PyMultivectorData_Funcs dense1_data_funcs = {
  .iter_next = (gaiternextfunc) dense1_iternext,
  .iter_init = (gaiterinitfunc) dense1_iterinit,
  .init = (gainitfunc) dense1_init,
};

static const PyMultivectorData_Funcs blades1_data_funcs = {
  .iter_next = (gaiternextfunc) blades1_iternext,
  .iter_init = (gaiterinitfunc) blades1_iterinit,
  .init = (gainitfunc) blades1_init,
};


static const PyMultivectorSubType dense1_subtype = {
    .math_funcs = dense1_math_funcs,
    .data_funcs = dense1_data_funcs,
    .name = "3DCGA",
    .type_name = "dense",
    .generated = 1,
    .metric = {1,1,1,1,-1,},
    .metric_size = 5,
    .ntype = 4,
};

static const PyMultivectorSubType blades1_subtype = {
    .math_funcs = blades1_math_funcs,
    .data_funcs = blades1_data_funcs,
    .name = "3DCGA",
    .type_name = "blades",
    .generated = 1,
    .metric = {1,1,1,1,-1,},
    .metric_size = 5,
    .ntype = 6,
};


PyMultivectorSubType gen_subtypes_array[4] = {
  dense0_subtype,
  blades0_subtype,
  dense1_subtype,
  blades1_subtype,
};
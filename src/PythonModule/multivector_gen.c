#include <Python.h>
#include "gasparse.h"
/* #include "multivector.h" */
#include "multivector_gen.h"

#define blades0zero {{0},{0},{0},{0}}


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

static int gen0_dense_gradeproject(gen0_DenseMultivector *dense, gen0_DenseMultivector *dense0, int *grades, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        gen0densegradeprojectfunc gradeproject =
                            gen0denseproject.gradeproject[grades[i]];
        if(gradeproject)
            gradeproject(dense0,dense);
        else return -1;
    }
    return 0;
}
static int gen0_blades_gradeproject(gen0_BladesMultivector *blades, gen0_BladesMultivector *blades0, int *grades, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        gen0bladesgradeprojectfunc gradeproject =
                            gen0bladesproject.gradeproject[grades[i]];
        if(gradeproject)
            gradeproject(blades0,blades);
        else return -1;
    }
    return 0;
}

static gen0_DenseMultivector dense0_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen0_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        if(bitmap[i] >= 8){
            return dense; // raise error
        }
        dense.value[bitmap[i]] += value[i]; // repeated blades are added to the same value
    }
    return dense;
}

static gen0_BladesMultivector blades0_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
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

 
static void* dense0_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen0_DenseMultivector *dense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    *dense = dense0_init_(bitmap,value,size,ga);
    return (void*)dense;
}

 
static void* blades0_init(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    gen0_BladesMultivector *blades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    *blades = blades0_init_(bitmap,value,size,ga);
    return (void*)blades;
}




static PyMultivectorIter dense0_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data = data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type.ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.size = 1;
    iter.niters = 8;
    iter.next = data->type.data_funcs.iter_next;
    iter.type_name = data->type.type_name;
    return iter;
}

static PyMultivectorIter blades0_iterinit(PyMultivectorObject *data){
    PyMultivectorIter iter;
    iter.data= data->data;
    iter.bitmap = -1;
    iter.value = 0;
    iter.type = data->type.ntype;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(2*sizeof(Py_ssize_t));
    iter.index[0] = 0;
    iter.index[1] = 0;
    iter.size = 2;
    iter.niters = 8;
    iter.next = data->type.data_funcs.iter_next;
    iter.type_name = data->type.type_name;
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
            return 1;
        case 1:
            iter->value = blades->value1[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap1[iter->index[1]];
            iter->grade = 1;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            return 1;
        case 2:
            iter->value = blades->value2[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap2[iter->index[1]];
            iter->grade = 2;
            if(++iter->index[1] >= 3){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            return 1;
        case 3:
            iter->value = blades->value3[iter->index[1]];
            iter->bitmap = gen0_gradesbitmap.bitmap3[iter->index[1]];
            iter->grade = 3;
            if(++iter->index[1] >= 1){
                iter->index[1] = 0;
                (*iter->index)++;
            }
            return 1;
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

// grade projection of the product of two multivectors <ab>_r
#define GEN0_DENSE_GRADE0GEOMETRICPRODUCT(dense,dense0,dense1){\
    dense.value[0] =\
    +dense0.value[0]*dense1.value[0]\
    +dense0.value[1]*dense1.value[1]\
    +dense0.value[2]*dense1.value[2]\
    -dense0.value[3]*dense1.value[3]\
    +dense0.value[4]*dense1.value[4]\
    -dense0.value[5]*dense1.value[5]\
    -dense0.value[6]*dense1.value[6]\
    -dense0.value[7]*dense1.value[7]\
;\
}
#define GEN0_DENSE_GRADE1GEOMETRICPRODUCT(dense,dense0,dense1){\
    dense.value[1] =\
    +dense0.value[0]*dense1.value[1]\
    +dense0.value[1]*dense1.value[0]\
    -dense0.value[2]*dense1.value[3]\
    +dense0.value[3]*dense1.value[2]\
    -dense0.value[4]*dense1.value[5]\
    +dense0.value[5]*dense1.value[4]\
    -dense0.value[6]*dense1.value[7]\
    -dense0.value[7]*dense1.value[6]\
;\
    dense.value[2] =\
    +dense0.value[0]*dense1.value[2]\
    +dense0.value[1]*dense1.value[3]\
    +dense0.value[2]*dense1.value[0]\
    -dense0.value[3]*dense1.value[1]\
    -dense0.value[4]*dense1.value[6]\
    +dense0.value[5]*dense1.value[7]\
    +dense0.value[6]*dense1.value[4]\
    +dense0.value[7]*dense1.value[5]\
;\
    dense.value[4] =\
    +dense0.value[0]*dense1.value[4]\
    +dense0.value[1]*dense1.value[5]\
    +dense0.value[2]*dense1.value[6]\
    -dense0.value[3]*dense1.value[7]\
    +dense0.value[4]*dense1.value[0]\
    -dense0.value[5]*dense1.value[1]\
    -dense0.value[6]*dense1.value[2]\
    -dense0.value[7]*dense1.value[3]\
;\
}
#define GEN0_DENSE_GRADE2GEOMETRICPRODUCT(dense,dense0,dense1){\
    dense.value[3] =\
    +dense0.value[0]*dense1.value[3]\
    +dense0.value[1]*dense1.value[2]\
    -dense0.value[2]*dense1.value[1]\
    +dense0.value[3]*dense1.value[0]\
    +dense0.value[4]*dense1.value[7]\
    -dense0.value[5]*dense1.value[6]\
    +dense0.value[6]*dense1.value[5]\
    +dense0.value[7]*dense1.value[4]\
;\
    dense.value[5] =\
    +dense0.value[0]*dense1.value[5]\
    +dense0.value[1]*dense1.value[4]\
    -dense0.value[2]*dense1.value[7]\
    +dense0.value[3]*dense1.value[6]\
    -dense0.value[4]*dense1.value[1]\
    +dense0.value[5]*dense1.value[0]\
    -dense0.value[6]*dense1.value[3]\
    -dense0.value[7]*dense1.value[2]\
;\
    dense.value[6] =\
    +dense0.value[0]*dense1.value[6]\
    +dense0.value[1]*dense1.value[7]\
    +dense0.value[2]*dense1.value[4]\
    -dense0.value[3]*dense1.value[5]\
    -dense0.value[4]*dense1.value[2]\
    +dense0.value[5]*dense1.value[3]\
    +dense0.value[6]*dense1.value[0]\
    +dense0.value[7]*dense1.value[1]\
;\
}
#define GEN0_DENSE_GRADE3GEOMETRICPRODUCT(dense,dense0,dense1){\
    dense.value[7] =\
    +dense0.value[0]*dense1.value[7]\
    +dense0.value[1]*dense1.value[6]\
    -dense0.value[2]*dense1.value[5]\
    +dense0.value[3]*dense1.value[4]\
    +dense0.value[4]*dense1.value[3]\
    -dense0.value[5]*dense1.value[2]\
    +dense0.value[6]*dense1.value[1]\
    +dense0.value[7]*dense1.value[0]\
;\
}

static gen0_DenseMultivector gen0_dense_gradegeometricproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1, int *grades, Py_ssize_t size){
    gen0_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_DENSE_GRADE0GEOMETRICPRODUCT(dense,dense0,dense1);
                break;
            case 1:
                GEN0_DENSE_GRADE1GEOMETRICPRODUCT(dense,dense0,dense1);
                break;
            case 2:
                GEN0_DENSE_GRADE2GEOMETRICPRODUCT(dense,dense0,dense1);
                break;
            case 3:
                GEN0_DENSE_GRADE3GEOMETRICPRODUCT(dense,dense0,dense1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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

// grade projection of the product of two multivectors <ab>_r
#define GEN0_DENSE_GRADE0INNERPRODUCT(dense,dense0,dense1){\
    dense.value[0] =\
    +dense0.value[1]*dense1.value[1]\
    +dense0.value[2]*dense1.value[2]\
    -dense0.value[3]*dense1.value[3]\
    +dense0.value[4]*dense1.value[4]\
    -dense0.value[5]*dense1.value[5]\
    -dense0.value[6]*dense1.value[6]\
    -dense0.value[7]*dense1.value[7]\
;\
}
#define GEN0_DENSE_GRADE1INNERPRODUCT(dense,dense0,dense1){\
    dense.value[1] =\
    -dense0.value[2]*dense1.value[3]\
    +dense0.value[3]*dense1.value[2]\
    -dense0.value[4]*dense1.value[5]\
    +dense0.value[5]*dense1.value[4]\
    -dense0.value[6]*dense1.value[7]\
    -dense0.value[7]*dense1.value[6]\
;\
    dense.value[2] =\
    +dense0.value[1]*dense1.value[3]\
    -dense0.value[3]*dense1.value[1]\
    -dense0.value[4]*dense1.value[6]\
    +dense0.value[5]*dense1.value[7]\
    +dense0.value[6]*dense1.value[4]\
    +dense0.value[7]*dense1.value[5]\
;\
    dense.value[4] =\
    +dense0.value[1]*dense1.value[5]\
    +dense0.value[2]*dense1.value[6]\
    -dense0.value[3]*dense1.value[7]\
    -dense0.value[5]*dense1.value[1]\
    -dense0.value[6]*dense1.value[2]\
    -dense0.value[7]*dense1.value[3]\
;\
}
#define GEN0_DENSE_GRADE2INNERPRODUCT(dense,dense0,dense1){\
    dense.value[3] =\
    +dense0.value[4]*dense1.value[7]\
    +dense0.value[7]*dense1.value[4]\
;\
    dense.value[5] =\
    -dense0.value[2]*dense1.value[7]\
    -dense0.value[7]*dense1.value[2]\
;\
    dense.value[6] =\
    +dense0.value[1]*dense1.value[7]\
    +dense0.value[7]*dense1.value[1]\
;\
}
#define GEN0_DENSE_GRADE3INNERPRODUCT(dense,dense0,dense1){\
}

static gen0_DenseMultivector gen0_dense_gradeinnerproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1, int *grades, Py_ssize_t size){
    gen0_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_DENSE_GRADE0INNERPRODUCT(dense,dense0,dense1);
                break;
            case 1:
                GEN0_DENSE_GRADE1INNERPRODUCT(dense,dense0,dense1);
                break;
            case 2:
                GEN0_DENSE_GRADE2INNERPRODUCT(dense,dense0,dense1);
                break;
            case 3:
                GEN0_DENSE_GRADE3INNERPRODUCT(dense,dense0,dense1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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

// grade projection of the product of two multivectors <ab>_r
#define GEN0_DENSE_GRADE0OUTERPRODUCT(dense,dense0,dense1){\
    dense.value[0] =\
    +dense0.value[0]*dense1.value[0]\
;\
}
#define GEN0_DENSE_GRADE1OUTERPRODUCT(dense,dense0,dense1){\
    dense.value[1] =\
    +dense0.value[0]*dense1.value[1]\
    +dense0.value[1]*dense1.value[0]\
;\
    dense.value[2] =\
    +dense0.value[0]*dense1.value[2]\
    +dense0.value[2]*dense1.value[0]\
;\
    dense.value[4] =\
    +dense0.value[0]*dense1.value[4]\
    +dense0.value[4]*dense1.value[0]\
;\
}
#define GEN0_DENSE_GRADE2OUTERPRODUCT(dense,dense0,dense1){\
    dense.value[3] =\
    +dense0.value[0]*dense1.value[3]\
    +dense0.value[1]*dense1.value[2]\
    -dense0.value[2]*dense1.value[1]\
    +dense0.value[3]*dense1.value[0]\
;\
    dense.value[5] =\
    +dense0.value[0]*dense1.value[5]\
    +dense0.value[1]*dense1.value[4]\
    -dense0.value[4]*dense1.value[1]\
    +dense0.value[5]*dense1.value[0]\
;\
    dense.value[6] =\
    +dense0.value[0]*dense1.value[6]\
    +dense0.value[2]*dense1.value[4]\
    -dense0.value[4]*dense1.value[2]\
    +dense0.value[6]*dense1.value[0]\
;\
}
#define GEN0_DENSE_GRADE3OUTERPRODUCT(dense,dense0,dense1){\
    dense.value[7] =\
    +dense0.value[0]*dense1.value[7]\
    +dense0.value[1]*dense1.value[6]\
    -dense0.value[2]*dense1.value[5]\
    +dense0.value[3]*dense1.value[4]\
    +dense0.value[4]*dense1.value[3]\
    -dense0.value[5]*dense1.value[2]\
    +dense0.value[6]*dense1.value[1]\
    +dense0.value[7]*dense1.value[0]\
;\
}

static gen0_DenseMultivector gen0_dense_gradeouterproduct(gen0_DenseMultivector dense0, gen0_DenseMultivector dense1, int *grades, Py_ssize_t size){
    gen0_DenseMultivector dense = {{0}};
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_DENSE_GRADE0OUTERPRODUCT(dense,dense0,dense1);
                break;
            case 1:
                GEN0_DENSE_GRADE1OUTERPRODUCT(dense,dense0,dense1);
                break;
            case 2:
                GEN0_DENSE_GRADE2OUTERPRODUCT(dense,dense0,dense1);
                break;
            case 3:
                GEN0_DENSE_GRADE3OUTERPRODUCT(dense,dense0,dense1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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

static gen0_DenseMultivector gen0_dense_dual(gen0_DenseMultivector dense0){
    gen0_DenseMultivector dense = {{0}};

    dense.value[7] = -dense0.value[0];
    dense.value[6] = -dense0.value[1];
    dense.value[5] = dense0.value[2];
    dense.value[4] = dense0.value[3];
    dense.value[3] = -dense0.value[4];
    dense.value[2] = -dense0.value[5];
    dense.value[1] = dense0.value[6];
    dense.value[0] = dense0.value[7];

    return dense;
}

static gen0_DenseMultivector gen0_dense_undual(gen0_DenseMultivector dense0){
    gen0_DenseMultivector dense = {{0}};

    dense.value[7] = dense0.value[0];
    dense.value[6] = dense0.value[1];
    dense.value[5] = -dense0.value[2];
    dense.value[4] = -dense0.value[3];
    dense.value[3] = dense0.value[4];
    dense.value[2] = dense0.value[5];
    dense.value[1] = -dense0.value[6];
    dense.value[0] = -dense0.value[7];

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

#define GEN0_BLADES_GRADE0GEOMETRICPRODUCT(blades,blades0,blades1){\
    blades.value0[0] =\
    +blades0.value0[0]*blades1.value0[0]\
    +blades0.value1[0]*blades1.value1[0]\
    +blades0.value1[1]*blades1.value1[1]\
    -blades0.value2[0]*blades1.value2[0]\
    +blades0.value1[2]*blades1.value1[2]\
    -blades0.value2[1]*blades1.value2[1]\
    -blades0.value2[2]*blades1.value2[2]\
    -blades0.value3[0]*blades1.value3[0]\
;\
}
#define GEN0_BLADES_GRADE1GEOMETRICPRODUCT(blades,blades0,blades1){\
    blades.value1[0] =\
    +blades0.value0[0]*blades1.value1[0]\
    +blades0.value1[0]*blades1.value0[0]\
    -blades0.value1[1]*blades1.value2[0]\
    +blades0.value2[0]*blades1.value1[1]\
    -blades0.value1[2]*blades1.value2[1]\
    +blades0.value2[1]*blades1.value1[2]\
    -blades0.value2[2]*blades1.value3[0]\
    -blades0.value3[0]*blades1.value2[2]\
;\
    blades.value1[1] =\
    +blades0.value0[0]*blades1.value1[1]\
    +blades0.value1[0]*blades1.value2[0]\
    +blades0.value1[1]*blades1.value0[0]\
    -blades0.value2[0]*blades1.value1[0]\
    -blades0.value1[2]*blades1.value2[2]\
    +blades0.value2[1]*blades1.value3[0]\
    +blades0.value2[2]*blades1.value1[2]\
    +blades0.value3[0]*blades1.value2[1]\
;\
    blades.value1[2] =\
    +blades0.value0[0]*blades1.value1[2]\
    +blades0.value1[0]*blades1.value2[1]\
    +blades0.value1[1]*blades1.value2[2]\
    -blades0.value2[0]*blades1.value3[0]\
    +blades0.value1[2]*blades1.value0[0]\
    -blades0.value2[1]*blades1.value1[0]\
    -blades0.value2[2]*blades1.value1[1]\
    -blades0.value3[0]*blades1.value2[0]\
;\
}
#define GEN0_BLADES_GRADE2GEOMETRICPRODUCT(blades,blades0,blades1){\
    blades.value2[0] =\
    +blades0.value0[0]*blades1.value2[0]\
    +blades0.value1[0]*blades1.value1[1]\
    -blades0.value1[1]*blades1.value1[0]\
    +blades0.value2[0]*blades1.value0[0]\
    +blades0.value1[2]*blades1.value3[0]\
    -blades0.value2[1]*blades1.value2[2]\
    +blades0.value2[2]*blades1.value2[1]\
    +blades0.value3[0]*blades1.value1[2]\
;\
    blades.value2[1] =\
    +blades0.value0[0]*blades1.value2[1]\
    +blades0.value1[0]*blades1.value1[2]\
    -blades0.value1[1]*blades1.value3[0]\
    +blades0.value2[0]*blades1.value2[2]\
    -blades0.value1[2]*blades1.value1[0]\
    +blades0.value2[1]*blades1.value0[0]\
    -blades0.value2[2]*blades1.value2[0]\
    -blades0.value3[0]*blades1.value1[1]\
;\
    blades.value2[2] =\
    +blades0.value0[0]*blades1.value2[2]\
    +blades0.value1[0]*blades1.value3[0]\
    +blades0.value1[1]*blades1.value1[2]\
    -blades0.value2[0]*blades1.value2[1]\
    -blades0.value1[2]*blades1.value1[1]\
    +blades0.value2[1]*blades1.value2[0]\
    +blades0.value2[2]*blades1.value0[0]\
    +blades0.value3[0]*blades1.value1[0]\
;\
}
#define GEN0_BLADES_GRADE3GEOMETRICPRODUCT(blades,blades0,blades1){\
    blades.value3[0] =\
    +blades0.value0[0]*blades1.value3[0]\
    +blades0.value1[0]*blades1.value2[2]\
    -blades0.value1[1]*blades1.value2[1]\
    +blades0.value2[0]*blades1.value1[2]\
    +blades0.value1[2]*blades1.value2[0]\
    -blades0.value2[1]*blades1.value1[1]\
    +blades0.value2[2]*blades1.value1[0]\
    +blades0.value3[0]*blades1.value0[0]\
;\
}

static gen0_BladesMultivector gen0_blades_gradegeometricproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1, int *grades, Py_ssize_t size){
    gen0_BladesMultivector blades = blades0zero;
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_BLADES_GRADE0GEOMETRICPRODUCT(blades,blades0,blades1);
                break;
            case 1:
                GEN0_BLADES_GRADE1GEOMETRICPRODUCT(blades,blades0,blades1);
                break;
            case 2:
                GEN0_BLADES_GRADE2GEOMETRICPRODUCT(blades,blades0,blades1);
                break;
            case 3:
                GEN0_BLADES_GRADE3GEOMETRICPRODUCT(blades,blades0,blades1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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

#define GEN0_BLADES_GRADE0INNERPRODUCT(blades,blades0,blades1){\
    blades.value0[0] =\
    +blades0.value1[0]*blades1.value1[0]\
    +blades0.value1[1]*blades1.value1[1]\
    -blades0.value2[0]*blades1.value2[0]\
    +blades0.value1[2]*blades1.value1[2]\
    -blades0.value2[1]*blades1.value2[1]\
    -blades0.value2[2]*blades1.value2[2]\
    -blades0.value3[0]*blades1.value3[0]\
;\
}
#define GEN0_BLADES_GRADE1INNERPRODUCT(blades,blades0,blades1){\
    blades.value1[0] =\
    -blades0.value1[1]*blades1.value2[0]\
    +blades0.value2[0]*blades1.value1[1]\
    -blades0.value1[2]*blades1.value2[1]\
    +blades0.value2[1]*blades1.value1[2]\
    -blades0.value2[2]*blades1.value3[0]\
    -blades0.value3[0]*blades1.value2[2]\
;\
    blades.value1[1] =\
    +blades0.value1[0]*blades1.value2[0]\
    -blades0.value2[0]*blades1.value1[0]\
    -blades0.value1[2]*blades1.value2[2]\
    +blades0.value2[1]*blades1.value3[0]\
    +blades0.value2[2]*blades1.value1[2]\
    +blades0.value3[0]*blades1.value2[1]\
;\
    blades.value1[2] =\
    +blades0.value1[0]*blades1.value2[1]\
    +blades0.value1[1]*blades1.value2[2]\
    -blades0.value2[0]*blades1.value3[0]\
    -blades0.value2[1]*blades1.value1[0]\
    -blades0.value2[2]*blades1.value1[1]\
    -blades0.value3[0]*blades1.value2[0]\
;\
}
#define GEN0_BLADES_GRADE2INNERPRODUCT(blades,blades0,blades1){\
    blades.value2[0] =\
    +blades0.value1[2]*blades1.value3[0]\
    +blades0.value3[0]*blades1.value1[2]\
;\
    blades.value2[1] =\
    -blades0.value1[1]*blades1.value3[0]\
    -blades0.value3[0]*blades1.value1[1]\
;\
    blades.value2[2] =\
    +blades0.value1[0]*blades1.value3[0]\
    +blades0.value3[0]*blades1.value1[0]\
;\
}
#define GEN0_BLADES_GRADE3INNERPRODUCT(blades,blades0,blades1){\
}

static gen0_BladesMultivector gen0_blades_gradeinnerproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1, int *grades, Py_ssize_t size){
    gen0_BladesMultivector blades = blades0zero;
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_BLADES_GRADE0INNERPRODUCT(blades,blades0,blades1);
                break;
            case 1:
                GEN0_BLADES_GRADE1INNERPRODUCT(blades,blades0,blades1);
                break;
            case 2:
                GEN0_BLADES_GRADE2INNERPRODUCT(blades,blades0,blades1);
                break;
            case 3:
                GEN0_BLADES_GRADE3INNERPRODUCT(blades,blades0,blades1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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

#define GEN0_BLADES_GRADE0OUTERPRODUCT(blades,blades0,blades1){\
    blades.value0[0] =\
    +blades0.value0[0]*blades1.value0[0]\
;\
}
#define GEN0_BLADES_GRADE1OUTERPRODUCT(blades,blades0,blades1){\
    blades.value1[0] =\
    +blades0.value0[0]*blades1.value1[0]\
    +blades0.value1[0]*blades1.value0[0]\
;\
    blades.value1[1] =\
    +blades0.value0[0]*blades1.value1[1]\
    +blades0.value1[1]*blades1.value0[0]\
;\
    blades.value1[2] =\
    +blades0.value0[0]*blades1.value1[2]\
    +blades0.value1[2]*blades1.value0[0]\
;\
}
#define GEN0_BLADES_GRADE2OUTERPRODUCT(blades,blades0,blades1){\
    blades.value2[0] =\
    +blades0.value0[0]*blades1.value2[0]\
    +blades0.value1[0]*blades1.value1[1]\
    -blades0.value1[1]*blades1.value1[0]\
    +blades0.value2[0]*blades1.value0[0]\
;\
    blades.value2[1] =\
    +blades0.value0[0]*blades1.value2[1]\
    +blades0.value1[0]*blades1.value1[2]\
    -blades0.value1[2]*blades1.value1[0]\
    +blades0.value2[1]*blades1.value0[0]\
;\
    blades.value2[2] =\
    +blades0.value0[0]*blades1.value2[2]\
    +blades0.value1[1]*blades1.value1[2]\
    -blades0.value1[2]*blades1.value1[1]\
    +blades0.value2[2]*blades1.value0[0]\
;\
}
#define GEN0_BLADES_GRADE3OUTERPRODUCT(blades,blades0,blades1){\
    blades.value3[0] =\
    +blades0.value0[0]*blades1.value3[0]\
    +blades0.value1[0]*blades1.value2[2]\
    -blades0.value1[1]*blades1.value2[1]\
    +blades0.value2[0]*blades1.value1[2]\
    +blades0.value1[2]*blades1.value2[0]\
    -blades0.value2[1]*blades1.value1[1]\
    +blades0.value2[2]*blades1.value1[0]\
    +blades0.value3[0]*blades1.value0[0]\
;\
}

static gen0_BladesMultivector gen0_blades_gradeouterproduct(gen0_BladesMultivector blades0, gen0_BladesMultivector blades1, int *grades, Py_ssize_t size){
    gen0_BladesMultivector blades = blades0zero;
    for(Py_ssize_t i = 0; i < size; i++){
        switch(grades[i]){
            case 0:
                GEN0_BLADES_GRADE0OUTERPRODUCT(blades,blades0,blades1);
                break;
            case 1:
                GEN0_BLADES_GRADE1OUTERPRODUCT(blades,blades0,blades1);
                break;
            case 2:
                GEN0_BLADES_GRADE2OUTERPRODUCT(blades,blades0,blades1);
                break;
            case 3:
                GEN0_BLADES_GRADE3OUTERPRODUCT(blades,blades0,blades1);
                break;
            default:// do nothing for different values
                break;
        }
    }
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


static gen0_BladesMultivector gen0_blades_dual(gen0_BladesMultivector blades0){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    blades.value3[0] = -blades0.value0[0];
    blades.value2[2] = -blades0.value1[0];
    blades.value2[1] =  blades0.value1[1];
    blades.value1[2] =  blades0.value2[0];
    blades.value2[0] = -blades0.value1[2];
    blades.value1[1] = -blades0.value2[1];
    blades.value1[0] =  blades0.value2[2];
    blades.value0[0] =  blades0.value3[0];
    return blades;
}

static gen0_BladesMultivector gen0_blades_undual(gen0_BladesMultivector blades0){
    gen0_BladesMultivector blades = {{0},{0},{0},{0},};

    blades.value3[0] =  blades0.value0[0];
    blades.value2[2] =  blades0.value1[0];
    blades.value2[1] = -blades0.value1[1];
    blades.value1[2] = -blades0.value2[0];
    blades.value2[0] =  blades0.value1[2];
    blades.value1[1] =  blades0.value2[1];
    blades.value1[0] = -blades0.value2[2];
    blades.value0[0] = -blades0.value3[0];
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


static PyMultivectorObject *binary_dense0_gradeproduct(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType ptype, GradeProjectMap gpmap){
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)data0->data;
    gen0_DenseMultivector *pdense1 = (gen0_DenseMultivector*)data1->data;
    gen0_DenseMultivector *pdense  = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));

    gen0_DenseMultivector projdense0 = {{0}};
    gen0_DenseMultivector projdense1 = {{0}};

    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pdense0 || !pdense1 || !pdense || !out){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise error
    }
    // grade projection of the input
    if(gen0_dense_gradeproject(&projdense0,pdense0,gpmap.grades0,gpmap.size0) == -1) return NULL;
    if(gen0_dense_gradeproject(&projdense1,pdense1,gpmap.grades0,gpmap.size0) == -1) return NULL;


    switch(ptype){
        case ProductType_geometric:
            *pdense = gen0_dense_gradegeometricproduct(projdense0,projdense1,gpmap.grades,gpmap.size);
            break;
        case ProductType_inner:
            *pdense = gen0_dense_gradeinnerproduct(projdense0,projdense1,gpmap.grades,gpmap.size);
            break;
        case ProductType_outer:
            *pdense = gen0_dense_gradeouterproduct(projdense0,projdense1,gpmap.grades,gpmap.size);
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
static PyMultivectorObject *binary_blades0_gradeproduct(PyMultivectorObject *data0, PyMultivectorObject *data1, ProductType ptype, GradeProjectMap gpmap){
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)data0->data;
    gen0_BladesMultivector *pblades1 = (gen0_BladesMultivector*)data1->data;
    gen0_BladesMultivector *pblades  = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));

    gen0_BladesMultivector projblades0 =  blades0zero;
    gen0_BladesMultivector projblades1 =  blades0zero;

    PyMultivectorObject *out = new_multivector(data0,-1);
    if(!pblades0 || !pblades1 || !pblades || !out){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise error
    }
    // grade projection of the input
    if(gen0_blades_gradeproject(&projblades0,pblades0,gpmap.grades0,gpmap.size0) == -1) return NULL;
    if(gen0_blades_gradeproject(&projblades1,pblades1,gpmap.grades0,gpmap.size0) == -1) return NULL;


    switch(ptype){
        case ProductType_geometric:
            *pblades = gen0_blades_gradegeometricproduct(projblades0,projblades1,gpmap.grades,gpmap.size);
            break;
        case ProductType_inner:
            *pblades = gen0_blades_gradeinnerproduct(projblades0,projblades1,gpmap.grades,gpmap.size);
            break;
        case ProductType_outer:
            *pblades = gen0_blades_gradeouterproduct(projblades0,projblades1,gpmap.grades,gpmap.size);
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

    if(gen0_dense_gradeproject(&dense,pdense0,grades,size) == -1) return NULL;

    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    *pdense = dense;
    out->data = (void*)pdense;
    PyMem_RawFree(grades);
    return out;
}
static PyMultivectorObject *unary_blades0_gradeproject(PyMultivectorObject *self, int *grades, Py_ssize_t size){
    PyMultivectorObject *out = NULL;
    gen0_BladesMultivector blades =  blades0zero;
    gen0_BladesMultivector *pblades;
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;

    if(gen0_blades_gradeproject(&blades,pblades0,grades,size) == -1) return NULL;

    out = new_multivector(self,-1); // pass -1 to inherit type of self
    pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    *pblades = blades;
    out->data = (void*)pblades;
    PyMem_RawFree(grades);
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
static PyMultivectorObject *unary_dense0_dual(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_dual(*pdense0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pdense;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *unary_dense0_undual(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_DenseMultivector *pdense0 = (gen0_DenseMultivector*)self->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!out || !pdense0 || !pdense){
        PyMem_RawFree(pdense);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pdense = gen0_dense_undual(*pdense0); // revert the order of the basis vectors of the multivector
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
static PyMultivectorObject *unary_blades0_dual(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_dual(*pblades0); // revert the order of the basis vectors of the multivector
    out->data = (void*)pblades;
    Py_SET_REFCNT(out,1);
    return out;
}
static PyMultivectorObject *unary_blades0_undual(PyMultivectorObject *self){
    PyMultivectorObject *out = new_multivector(self,-1);
    gen0_BladesMultivector *pblades0 = (gen0_BladesMultivector*)self->data;
    gen0_BladesMultivector *pblades = (gen0_BladesMultivector*)PyMem_RawMalloc(sizeof(gen0_BladesMultivector));
    if(!out || !pblades0 || !pblades){
        PyMem_RawFree(pblades);
        free_multivector(out);
        return NULL; // raise memory error
    }
    *pblades = gen0_blades_undual(*pblades0); // revert the order of the basis vectors of the multivector
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
    .dual = (gaunaryfunc) unary_dense0_dual,
    .undual = (gaunaryfunc) unary_dense0_undual,
    .ternary_product = (gaternaryprodfunc) ternary_dense0_product,
    .graded_product = (gabinarygradefunc) binary_dense0_gradeproduct,
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
    .dual = (gaunaryfunc) unary_blades0_dual,
    .undual = (gaunaryfunc) unary_blades0_undual,
    .ternary_product = (gaternaryprodfunc) ternary_blades0_product,
    .graded_product = (gabinarygradefunc) binary_blades0_gradeproduct,
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
    .type_name = "dense0",
    .generated = 1,
    .metric = {1,1,1,},
    .msize = 3,
    .ntype = 3,
    .asize = 8,
};

static const PyMultivectorSubType blades0_subtype = {
    .math_funcs = blades0_math_funcs,
    .data_funcs = blades0_data_funcs,
    .name = "3DVGA",
    .type_name = "blades0",
    .generated = 1,
    .metric = {1,1,1,},
    .msize = 3,
    .ntype = 4,
    .asize = 8,
};


PyMultivectorSubType gen_subtypes_array[2] = {
  dense0_subtype,
  blades0_subtype,
};
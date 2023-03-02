#ifndef GRADE_SPARSE_H_
#define GRADE_SPARSE_H_
#include "sparse.h"

typedef struct grade_info{
    size_t *grade_size;// number of basis vectors by grade
    size_t max_grade;
}grade_info;

typedef struct blades{
    sparse *data;
    unsigned int *grades;
    size_t size;
}blades;

typedef struct graded_multivectors{
    blades *data;
    graded_map m;
    size_t size;
    float precision;
    grade_info ginfo;
}graded_multivectors;


blades initialize_blades(grade_info ginfo);

#endif // GRADE_SPARSE_H_

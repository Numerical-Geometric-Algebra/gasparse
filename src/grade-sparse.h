#ifndef GRADE_SPARSE_H_
#define GRADE_SPARSE_H_
#include "sparse.h"

typedef struct blades{
    sparse *data;
    unsigned int *grades;
}blades;

typedef struct graded_multivectors{
    blades *data;
    map m;
    size_t size;
    float precision;
}graded_multivectors;

#endif // GRADE_SPARSE_H_

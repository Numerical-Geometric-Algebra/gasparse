#ifndef GRADE_SPARSE_H_
#define GRADE_SPARSE_H_
#include "sparse.h"

typedef struct blades{
    sparse *data;
    unsigned int *grade;
    size_t size;
}blades;

typedef struct graded_multivectors{
    blades *data;
    map m;
    grade_map gm;
    size_t size;
    float precision;
}graded_multivectors;


blades initialize_blades_empty(size_t);
blades initialize_blades(unsigned int *, size_t);
blades graded_product(graded_multivectors);
unsigned int *initialize_grade_size(grade_map);
blades sparse_to_graded(sparse,grade_map);
void free_blades(blades y);
#endif // GRADE_SPARSE_H_

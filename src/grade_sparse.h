#ifndef GRADE_SPARSE_H_
#define GRADE_SPARSE_H_
#include "sparse.h"

typedef struct blades{
    sparse *data;
    unsigned int *grade;
    size_t size;
}blades;

typedef struct graded_multivectors{
    blades a;
    blades b;
    map m;
    grade_map gm;
    float precision;
}graded_multivectors;


blades initialize_blades_empty(size_t);
blades initialize_blades(unsigned int *,size_t);
blades graded_product(graded_multivectors);
blades graded_product_(blades,blades,map,grade_map,float);
blades graded_scalar_multiply(float,blades);
unsigned int *initialize_grade_size(grade_map);
blades sparse_to_graded(sparse,grade_map);
void free_blades(blades);
void graded_remove_small(blades,float,unsigned int*);
blades grade_dense_to_grade_sparse(blades,unsigned int *);
blades graded_general_product(graded_multivectors,project_map);
blades graded_general_product_(blades,blades,map,grade_map,project_map,float);
blades grade_sparse_grade_project(blades,unsigned int *,size_t,size_t);

blades graded_atomic_add_add_(blades **,size_t,grade_map,float);
blades graded_add_add_(blades,blades,grade_map,float);
blades graded_atomic_add_append(blades **,size_t);
blades graded_add_append(blades,blades);

blades graded_add_add(graded_multivectors);

#endif // GRADE_SPARSE_H_

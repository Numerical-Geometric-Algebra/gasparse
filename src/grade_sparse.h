#ifndef GRADE_SPARSE_H_
#define GRADE_SPARSE_H_
#include "sparse.h"

typedef struct blades{
    sparse *data;
    size_t *grade;
    size_t size;
}blades;

typedef struct graded_multivectors{
    blades a;
    blades b;
    map m;
    grade_map gm;
    float precision;
}graded_multivectors;

typedef struct blades_extra{
    map *m;
    size_t size;
    grade_map gm;
    float precision;
}blades_extra;


blades initialize_blades_empty(size_t);
blades initialize_blades(size_t *,size_t);
blades graded_product(graded_multivectors);
blades graded_scalar_multiply(float,blades);
size_t *initialize_grade_size(grade_map);
blades sparse_to_graded(sparse,grade_map);
void graded_remove_small(blades,float,size_t*);
blades grade_dense_to_grade_sparse(blades,size_t *);
blades graded_general_product(graded_multivectors,project_map);
blades grade_sparse_grade_project(blades,size_t *,size_t,size_t);

blades graded_atomic_add_append(blades **,size_t);
blades graded_add_append(blades,blades);

blades blades_product_(blades,blades,map,grade_map,float);
void *blades_operator_product__(void*,void*,void*,size_t);
blades blades_general_product_(blades,blades,map,grade_map,project_map,float);
blades blades_atomic_add_add_(blades **,size_t,grade_map,float);
blades blades_add_add_(blades,blades,grade_map,float);
void blades_free_(blades);

blades graded_add_add(graded_multivectors);

// functions to use in einsum
void blades_assign__(void*,void*);
void blades_init__(void*,size_t);
void *blades_product__(void*,void*,void*);
void *blades_add_add__(void*,void*,void*);
void *blades_atomic_add__(void*,size_t,void*);
void blades_free__(void*,size_t);
#endif // GRADE_SPARSE_H_

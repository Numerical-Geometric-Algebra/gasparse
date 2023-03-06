#ifndef MAIN_H_
#define MAIN_H_
#include "cayley.h"
#include "sparse.h"
#include "grade_sparse.h"
#include "dense.h"

void time_generator(void);
void print_generator_table(unsigned int,unsigned int,unsigned int);

void print_map(map);
void print_dense(dense,int,char*);
void print_sparse(sparse,char*);
void print_blades(blades,char*);
void print_all_types(sparse,dense,blades);

void compute_product(void);
void compute_graded_product(void);
void grade_project(void);
void test_all_types(void);

dense project_dense_product(dense_multivectors,project_map,dense_grade_map);
sparse project_sparse_product(sparse_multivectors,project_map,dense_grade_map);
blades project_blades_product(graded_multivectors,project_map);
#endif // MAIN_H_

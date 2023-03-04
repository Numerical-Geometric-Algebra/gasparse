#include "main.h"
#include "cayley.h"
#include "sparse.h"
#include "grade_sparse.h"
#include "dense.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_METRIC_SIZE 32

int main(){

    /* time_generator(); */
    /* print_generator_table(3,2,1); */
    /* compute_product(); */
    /* printf("\n"); */
    /* grade_map m = bitmap_grade_map(32); */
    /* for(size_t i = 0; i < m.size; i++){ */
    /*     printf("%zu,%d,%d\n",i,m.position[i],m.grade[i]); */
    /* } */
    /* printf("\n"); */
    /* for(size_t i = 0; i < m.max_grade; i++){ */
    /*     printf("%zu,%d\n",i,m.grade_size[i]); */
    /* } */

    compute_graded_product();
    /* grade_project(); */
}


void time_generator(void){
    clock_t start, end;
    double cpu_time_used;
    size_t p,q,r;

    // For metric with more than 14 basis vector the cpu runs out of memory
    for(size_t n = 1; n <= 14; n++){
        p = n; q = 0; r = 0;
        start = clock();
        map m = cayley_table(p,q,r);
        end = clock();
        cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
        printf("(p,q,r)=(%zu,%zu,%zu) | time=%f\n",p,q,r,cpu_time_used);
        free_map(m);
    }
}


void print_generator_table(unsigned int p,unsigned int q,unsigned int r){

    map m = cayley_table(p,q,r);
    print_map(m);
    free_map(m);

}

void grade_project(void){
    unsigned int p=4, q=1, r=1;

    sparse_multivectors mvs;
    graded_multivectors gmvs;
    dense_multivectors dmvs;
    float precision = 1e-12;

    int size = 5;
    float value1[] = {1,2,3,4,6};
    float value0[] = {3,4,7,3,9};

    int bitmap1[] = {1,3,7,15,31};
    int bitmap0[] = {1,3,7,15,31};

    // Generate the cayley table
    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    mvs.m = m;
    mvs.size = 2;
    mvs.precision = precision;
    mvs.data = (sparse*)malloc(mvs.size*sizeof(sparse));
    mvs.data[0].bitmap = bitmap0;
    mvs.data[1].bitmap = bitmap1;
    mvs.data[0].value = value0;
    mvs.data[1].value = value1;
    mvs.data[0].size = size;
    mvs.data[1].size = size;

    /* sparse y = geometric_product(mvs); */
    blades x0 = sparse_to_graded(mvs.data[0],gm);
    blades x1 = sparse_to_graded(mvs.data[1],gm);
    gmvs.size = 2;
    gmvs.data = (blades*)malloc(gmvs.size*sizeof(blades));
    gmvs.data[0] = x0;
    gmvs.data[1] = x1;
    gmvs.gm = gm;
    gmvs.m = m;
    gmvs.precision = precision;

    blades y = graded_product(gmvs);

    map inv_m = invert_map(m);
    dense dense_x0 = sparse_to_dense(mvs.data[0],m.size);
    dense dense_x1 = sparse_to_dense(mvs.data[1],m.size);
    dmvs.size = gmvs.size;
    dmvs.data = (dense*)malloc(dmvs.size*sizeof(dense));
    dmvs.data[0] = dense_x0;
    dmvs.data[1] = dense_x1;
    dmvs.m = inv_m;
    dense dense_y = inverse_dense_product(dmvs);
    unsigned int project_grade[3] = {2,3,5};

    dense proj_mv = dense_grade_project(dense_x0,gm.grade,project_grade,gm.size,3);
    sparse proj_mv_sparse = sparse_grade_project(mvs.data[0],gm.grade,project_grade,gm.size,3);

    for(unsigned int i = 0; i < dense_x0.size; i++){
        printf("%d,%.0f\n",i,dense_x0.value[i]);
    }
    printf("\n");

    printf("grade projection of dense:\n");
    for(unsigned int i = 0; i < proj_mv.size; i++){
        printf("%d,%.0f\n",i,proj_mv.value[i]);
    }
    printf("\n");
    printf("grade projection of sparse\n");
    for(unsigned int i = 0; i < proj_mv_sparse.size; i++){
        printf("%d,%.0f\n",proj_mv_sparse.bitmap[i],proj_mv_sparse.value[i]);
    }


    free_blades(x0);
    free_blades(x1);
    free_blades(y);
    free_grade_map(gm);
    free(mvs.data);
    free(gmvs.data);
    free_map(m);
    free_map(inv_m);
    free(dense_x0.value);
    free(dense_x1.value);
    free(dense_y.value);
    free(proj_mv.value);
    free(dmvs.data);
    free_sparse(proj_mv_sparse);

}


void compute_graded_product(void){
    unsigned int p=4, q=1, r=1;

    sparse_multivectors mvs;
    graded_multivectors gmvs;
    dense_multivectors dmvs;
    float precision = 1e-12;

    int size = 6;
    float value1[] = {1,2,3,4,4,6};
    float value0[] = {3,4,7,3,8,3};

    int bitmap1[] = {1,2,7,9,6,31};
    int bitmap0[] = {0,1,3,7,15,31};

    // Generate the cayley table
    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    mvs.m = m;
    mvs.size = 2;
    mvs.precision = precision;
    mvs.data = (sparse*)malloc(mvs.size*sizeof(sparse));
    mvs.data[0].bitmap = bitmap0;
    mvs.data[1].bitmap = bitmap1;
    mvs.data[0].value = value0;
    mvs.data[1].value = value1;
    mvs.data[0].size = size;
    mvs.data[1].size = size;

    /* sparse y = geometric_product(mvs); */
    blades x0 = sparse_to_graded(mvs.data[0],gm);
    blades x1 = sparse_to_graded(mvs.data[1],gm);
    gmvs.size = 2;
    gmvs.data = (blades*)malloc(gmvs.size*sizeof(blades));
    gmvs.data[0] = x0;
    gmvs.data[1] = x1;
    gmvs.gm = gm;
    gmvs.m = m;
    gmvs.precision = precision;

    unsigned int project_grade[4] = {3,2,5,0};
    blades proj_x0 = grade_sparse_grade_project(x0,project_grade,m.size,4);


    printf("\nprojected graded\n");
    for(size_t i = 0; i < proj_x0.size; i++){
        for(size_t j = 0; j < proj_x0.data[i].size; j++){
            printf("(%d,%d,%.0f)\n",proj_x0.grade[i],proj_x0.data[i].bitmap[j],proj_x0.data[i].value[j]);
        }
        printf("\n");
    }

    blades y = graded_product(gmvs);

    /* for(size_t i = 0; i < y.size; i++){ */
    /*     for(size_t j = 0; j < y.data[i].size; j++){ */
    /*         printf("(%d,%d,%.0f)\n",y.grade[i],y.data[i].bitmap[j],y.data[i].value[j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

    map inv_m = invert_map(m);
    dense dense_x0 = sparse_to_dense(mvs.data[0],m.size);
    dense dense_x1 = sparse_to_dense(mvs.data[1],m.size);
    dmvs.size = gmvs.size;
    dmvs.data = (dense*)malloc(dmvs.size*sizeof(dense));
    dmvs.data[0] = dense_x0;
    dmvs.data[1] = dense_x1;
    dmvs.m = inv_m;
    dense dense_y = inverse_dense_product(dmvs);

    /* printf("\n"); */
    /* for(size_t i = 0; i < dense_y.size; i++){ */
    /*     printf("(%zu,%f)\n",i,dense_y.value[i]); */
    /* } */

    free_blades(x0);
    free_blades(x1);
    free_blades(y);
    free_blades(proj_x0);
    free_grade_map(gm);
    free(mvs.data);
    free(gmvs.data);
    free_map(m);
    free_map(inv_m);
    free(dense_x0.value);
    free(dense_x1.value);
    free(dense_y.value);
    free(dmvs.data);

}


void compute_product(void){
    unsigned int p=4, q=1, r=1;

    sparse_multivectors mvs;
    float precision = 1e-12;

    int size = 2;
    float value1[] = {1,2,3,4};
    float value0[] = {1,2,7,3};

    int bitmap1[] = {1,2,7,9};
    int bitmap0[] = {1,2,0,35};

    // Generate the cayley table
    map m = cayley_table(p,q,r);
    mvs.m = m;
    mvs.size = 2;
    mvs.precision = precision;
    mvs.data = (sparse*)malloc(mvs.size*sizeof(sparse));
    mvs.data[0].bitmap = (int*)malloc(size*sizeof(int));
    mvs.data[1].bitmap = (int*)malloc(size*sizeof(int));
    mvs.data[0].value = (float*)malloc(size*sizeof(float));
    mvs.data[1].value = (float*)malloc(size*sizeof(float));
    mvs.data[0].size = size;
    mvs.data[1].size = size;

    for(int i = 0; i < size; i++){
        mvs.data[0].bitmap[i] = bitmap0[i];
        mvs.data[1].bitmap[i] = bitmap1[i];
        mvs.data[0].value[i] = value0[i];
        mvs.data[1].value[i] = value1[i];
    }

    sparse y = sparse_product(mvs);

    for (unsigned int i=0; i < y.size; i++) {
        printf("%d,%.0f\n",y.bitmap[i],y.value[i]);
    }

    free(mvs.data[0].bitmap);
    free(mvs.data[1].bitmap);
    free(mvs.data[0].value);
    free(mvs.data[1].value);
    free(mvs.data);
    free(y.bitmap);
    free(y.value);
    free_map(m);

}

#include "main.h"
#include "cayley.h"
#include "sparse.h"
#include "grade_sparse.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_METRIC_SIZE 32

int main(){

    /* time_generator(); */
    /* print_generator_table(3,2,1); */
    compute_product();
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

void compute_graded_product(void){
    unsigned int p=4, q=1, r=1;

    sparse_multivectors mvs;
    graded_multivectors gmvs;
    float precision = 1e-12;

    int size = 2;
    float value1[] = {1,2,3,4};
    float value0[] = {1,2,7,3};

    int bitmap1[] = {1,2,7,9};
    int bitmap0[] = {1,2,0,35};

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

    for(size_t i = 0; i < y.size; i++){
        for(size_t j = 0; j < y.data[i].size; j++){
            printf("(%d,%d,%.0f)\n",y.grade[i],y.data[i].bitmap[j],y.data[i].value[j]);
        }
        printf("\n");
    }
    free_blades(x0);
    free_blades(x1);
    free_blades(y);
    free_grade_map(gm);
    free(mvs.data);
    free(gmvs.data);
    free_map(m);

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

    sparse y = geometric_product(mvs);

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

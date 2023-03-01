#include "main.h"
#include "cayley.h"
#include "sparse.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_METRIC_SIZE 32

int main(){

    /* time_generator(); */
    /* print_generator_table(3,2,1); */
    compute_product();
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

void compute_product(void){
    unsigned int p=4, q=1, r=1;

    sparse_multivectors mvs;
    float precision = 1e-12;

    int size = 2;
    float value1[] = {1,2,3,4};
    float value0[] = {5,9,7,3};

    int bitmap1[] = {32,4,7,9};
    int bitmap0[] = {33,2,0,33};

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

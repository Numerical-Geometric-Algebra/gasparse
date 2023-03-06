#include "main.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_METRIC_SIZE 32

int main(){
    test_all_types();
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


void test_all_types(void){
    unsigned int p = 4, q = 1, r = 1;

    float precision = 1e-12;
    int size = 6;
    int mvs_size = 2;
    size_t l_size = 2;
    size_t r_size = 3;
    size_t k_size = 1;

    float value_l[] = {1,4,6,8,2,3};
    int bitmap_l[] = {1,4,5,3,7,6};

    float value_r[] = {2,7,33,2,4,9};
    int bitmap_r[] = {3,12,2,32,5,27};

    unsigned int l_grade[] = {1,3};
    unsigned int r_grade[] = {2,5};
    unsigned int k_grade[] = {1};

    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    sparse data[] = {{bitmap_l,value_l,size},
                     {bitmap_r,value_r,size}};
    sparse_multivectors mvs = {data,m,mvs_size,precision};

    blades xl = sparse_to_graded(data[0],gm);
    blades xr = sparse_to_graded(data[1],gm);
    dense dense_xl = sparse_to_dense(mvs.data[0],m.size);
    dense dense_xr = sparse_to_dense(mvs.data[1],m.size);

    blades data_blades[] = {xl,xr};
    dense data_dense[] = {dense_xl,dense_xr};

    graded_multivectors gmvs = {data_blades,m,gm,size,precision};
    dense_multivectors dmvs = {data_dense,m,size};

    project_map pm = {l_grade,l_size,r_grade,r_size,k_grade,k_size};
    dense_grade_map dgm = {gm.max_grade,gm.size,gm.grade};

    sparse sparse_y = project_sparse_product(mvs,pm,dgm);
    blades blades_y = project_blades_product(gmvs,pm);
    dense dense_y = project_dense_product(dmvs,pm,dgm);

    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);

    sparse_y = sparse_general_product(mvs,pm,dgm);
    blades_y = graded_general_product(gmvs,pm);
    dense_y = project_dense_product(dmvs,pm,dgm);

    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);


    free(dense_xl.value);
    free(dense_xr.value);
    free_blades(xl);
    free_blades(xr);
    free_map(m);
    free_grade_map(gm);
}


void print_all_types(sparse sparse_y, dense dense_y, blades blades_y){
    char tab[] = "\t";
    printf("Dense:\n");
    print_dense(dense_y,0,tab);
    printf("Sparse:\n");
    print_sparse(sparse_y,tab);
    printf("Blades:\n");
    print_blades(blades_y,tab);
}

void print_dense(dense y, int non_zero,char *s){
    for(size_t i = 0; i < y.size; i++ ){
        if(non_zero || y.value[i] != 0){
            printf("%s(%zu,%.0f)\n",s,i,y.value[i]);
        }
    }
}

void print_sparse(sparse y,char *s){
    for(size_t i = 0; i < y.size; i++ )
        printf("%s(%d,%.0f)\n",s,y.bitmap[i],y.value[i]);
}

void print_blades(blades y,char *s){
    for(size_t i = 0; i < y.size; i++){
        printf("%sGrade %d:\n",s, y.grade[i]);
        for(size_t j = 0; j < y.data[i].size; j++){
            printf("%s\t(%d,%.0f)\n",s,y.data[i].bitmap[j],y.data[i].value[j]);
        }
        printf("\n");
    }
}

dense project_dense_product(dense_multivectors mvs, project_map pm, dense_grade_map dgm){
    dense y,proj_y;
    dense xl = dense_grade_project(mvs.data[0],pm.l,pm.l_size,dgm);
    dense xr = dense_grade_project(mvs.data[1],pm.r,pm.r_size,dgm);
    dense data[] = {xl,xr};
    mvs.data = data;
    y = dense_product(mvs);
    proj_y = dense_grade_project(y,pm.k,pm.k_size,dgm);

    free(xl.value);
    free(xr.value);
    free(y.value);

    return proj_y;
}

sparse project_sparse_product(sparse_multivectors mvs, project_map pm, dense_grade_map dgm){
    sparse y,proj_y;
    sparse xl = sparse_grade_project(mvs.data[0],pm.l,pm.l_size,dgm);
    sparse xr = sparse_grade_project(mvs.data[1],pm.r,pm.r_size,dgm);
    sparse data[] = {xl,xr};
    mvs.data = data;
    y = sparse_product(mvs);
    proj_y = sparse_grade_project(y,pm.k,pm.k_size,dgm);

    free_sparse(xl);
    free_sparse(xr);
    free_sparse(y);

    return proj_y;
}

blades project_blades_product(graded_multivectors mvs, project_map pm){
    blades y,proj_y;
    blades xl = grade_sparse_grade_project(mvs.data[0],pm.l,pm.l_size,mvs.gm.size);
    blades xr = grade_sparse_grade_project(mvs.data[1],pm.r,pm.r_size,mvs.gm.size);
    blades data[] = {xl,xr};
    mvs.data = data;
    y = graded_product(mvs);
    proj_y = grade_sparse_grade_project(y,pm.k,pm.k_size,mvs.gm.size);

    free_blades(xl);
    free_blades(xr);
    free_blades(y);

    return proj_y;
}


void print_map(map m){
     for(size_t i = 0; i < m.size; i++ ){
        for(size_t j = 0; j < m.size; j++ ){
            printf(" (%zu,%zu) | ", i, j);
            printf("(%d,%d)\n", m.sign[i][j], m.bitmap[i][j]);
        }
        printf("\n");
    }
}

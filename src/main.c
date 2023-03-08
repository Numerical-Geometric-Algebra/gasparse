#include "main.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_METRIC_SIZE 32

int main(){
    /* test_product_all_types(); */
    /* test_sum_all_types(); */
    /* test_scalar_multiply_all_types(); */
    test_parse();
    /* test_matrix_mult(); */
}


void test_parse(void){
    /* labels l = parse_subscripts("abbcbc",6,6); */
    /* labels l = parse_subscripts("ab...bc",7,6); */
    /* free(l.op_labels); */
    size_t ndims[] = {2,2,2};
    char args[] = "ik,kj->ij";
    symbols s = parse_all(args,strlen(args),ndims,3);

    for(size_t i = 0; i < s.size_; i++){
        for(size_t j = 0; j< s.size[i]; j++){
            printf("%d,",s.subscripts[i][j]);
        }printf("\n");
    }

    free_symbols(s);
}

void time_generator(void){
    clock_t start, end;
    double cpu_time_used;
    size_t p,q,r;

    // For metric with more than 14 basis vectors the cpu runs out of memory
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

void test_scalar_multiply_all_types(void){
    unsigned int p = 4, q = 1, r = 1;

    float precision = 1e-12;
    int size = 3;
    float scalar = 8.763;

    float value_l[] = {1,4,6,8,2,3};
    int bitmap_l[] =  {1,4,0,3,7,6};

    float value_r[] = {-1,7, 33,2,4,9};
    int bitmap_r[] =  { 1,12,4, 32,5,27};

    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    sparse a = {bitmap_l,value_l,size},b = {bitmap_r,value_r,size};

    dense_grade_map dgm = {gm.max_grade,gm.size,gm.grade};
    sparse_multivectors mvs = {a,b,m,precision,dgm};

    blades xl = sparse_to_graded(a,gm);
    blades xr = sparse_to_graded(b,gm);
    dense dense_xl = sparse_to_dense(mvs.a,m.size);
    dense dense_xr = sparse_to_dense(mvs.b,m.size);

    graded_multivectors gmvs = {xl,xr,m,gm,precision};
    dense_multivectors dmvs = {dense_xl,dense_xr,m,dgm};

    sparse sparse_y = sparse_scalar_multiply(scalar,mvs.a);
    blades blades_y = graded_scalar_multiply(scalar,gmvs.a);
    dense dense_y = dense_scalar_multiply(scalar,dmvs.a);

    printf("scalar multiply:\n");
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


void test_sum_all_types(void){
    unsigned int p = 4, q = 1, r = 1;

    float precision = 1e-12;
    int size = 3;

    float value_l[] = {1,4,6,8,2,3};
    int bitmap_l[] =  {1,4,0,3,7,6};

    float value_r[] = {-1,7, 33,2,4,9};
    int bitmap_r[] =  { 1,12,4, 32,5,27};

    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    sparse a = {bitmap_l,value_l,size},b = {bitmap_r,value_r,size};

    dense_grade_map dgm = {gm.max_grade,gm.size,gm.grade};
    sparse_multivectors mvs = {a,b,m,precision,dgm};

    blades xl = sparse_to_graded(a,gm);
    blades xr = sparse_to_graded(b,gm);
    dense dense_xl = sparse_to_dense(mvs.a,m.size);
    dense dense_xr = sparse_to_dense(mvs.b,m.size);

    graded_multivectors gmvs = {xl,xr,m,gm,precision};
    dense_multivectors dmvs = {dense_xl,dense_xr,m,dgm};

    sparse sparse_y = sparse_add_add(mvs);
    blades blades_y = graded_add_add(gmvs);
    dense dense_y = dense_add(dmvs.a,dmvs.b);

    printf("add_add:\n");
    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);

    sparse_y = sparse_add_append(mvs.a,mvs.b);
    blades_y = graded_add_append(gmvs.a,gmvs.b);
    dense_y = dense_add(dmvs.a,dmvs.b);

    printf("add_append:\n");
    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);

    sparse **s_mv = (sparse**)malloc(2*sizeof(sparse*));
    s_mv[0] = &mvs.a;
    s_mv[1] = &mvs.b;

    blades **b_mv = (blades**)malloc(2*sizeof(blades*));
    b_mv[0] = &gmvs.a;
    b_mv[1] = &gmvs.b;

    dense **d_mv = (dense**)malloc(2*sizeof(dense*));
    d_mv[0] = &dmvs.a;
    d_mv[1] = &dmvs.b;

    sparse_y = sparse_atomic_add_append(s_mv,2);
    blades_y = graded_atomic_add_append(b_mv,2);
    dense_y = dense_atomic_add(d_mv,2);

    printf("atomic_add_append:\n");
    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);

    sparse_y = sparse_atomic_add_add_(s_mv,2,m.size,precision);
    blades_y = graded_atomic_add_add_(b_mv,2,gm,precision);
    dense_y = dense_atomic_add(d_mv,2);

    printf("atomic_add_add:\n");
    print_all_types(sparse_y,dense_y,blades_y);


    free(s_mv);
    free(b_mv);
    free(d_mv);
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


void test_matrix_mult(void){
    unsigned int p = 4, q = 1, r = 1;

    float precision = 1e-12;
    int size = 1;

    float value_11[] = {1,4,6,8,2,3};
    int bitmap_11[] = {0,4,7,3,7,6};

    float value_12[] = {1,4,8,8,2,3};
    int bitmap_12[] = {1,4,5,3,7,6};

    float value_21[] = {1,4,6,8,21,3};
    int bitmap_21[] = {2,4,5,3,7,6};

    float value_22[] = {1,4,6,8,2,3};
    int bitmap_22[] = {3,4,5,3,32,6};

    float value_1[] = {1,7,33,2,4,9};
    int bitmap_1[] = {0,12,2,32,12,27};

    float value_2[] = {1,12,33,2,4,9};
    int bitmap_2[] = {1,5,2,32,17,27};

    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);
    sparse a[] = {{bitmap_11,value_11,size},{bitmap_12,value_12,size},
                  {bitmap_21,value_21,size},{bitmap_22,value_22,size}};
    sparse b[] = {{bitmap_1,value_1,size},{bitmap_2,value_2,size}};


    blades x[4]; // 2x2 matrix
    blades y[2]; // 2 vector
    for(size_t i = 0; i < 4; i++)
        x[i] = sparse_to_graded(a[i],gm);

    for(size_t i = 0; i < 2; i++)
        y[i] = sparse_to_graded(b[i],gm);

    blades **data = (blades**)malloc(2*sizeof(blades*));
    size_t **shapes = (size_t**)malloc(2*sizeof(size_t*));
    data[0] = x; data[1] = y;


    shapes[0] = (size_t*)malloc(2*sizeof(size_t));
    shapes[0][0] = 2;
    shapes[0][1] = 2;

    graded_tensor_multivectors tmvs = {data,shapes,NULL,1,m,gm,precision};

    graded_tensor out = vector_matrix_mult(tmvs);

    for(size_t i = 0; i < 2; i++)
        print_blades(out.data[i],"");


    for(size_t i = 0; i < 4; i++)
        free_blades(x[i]);

    free(shapes[0]);
    for(size_t i = 0; i < 2; i++){
        free_blades(y[i]);
        free_blades(out.data[i]);
    }
    free(data);
    free(out.data);
    free(out.shapes);
    free(shapes);
    free_map(m);
    free_grade_map(gm);

}


void test_product_all_types(void){
    unsigned int p = 4, q = 1, r = 1;

    float precision = 1e-12;
    int size = 6;
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
    sparse a = {bitmap_l,value_l,size},b = {bitmap_r,value_r,size};

    dense_grade_map dgm = {gm.max_grade,gm.size,gm.grade};
    sparse_multivectors mvs = {a,b,m,precision,dgm};

    blades xl = sparse_to_graded(a,gm);
    blades xr = sparse_to_graded(b,gm);
    dense dense_xl = sparse_to_dense(mvs.a,m.size);
    dense dense_xr = sparse_to_dense(mvs.b,m.size);

    graded_multivectors gmvs = {xl,xr,m,gm,precision};
    dense_multivectors dmvs = {dense_xl,dense_xr,m,dgm};

    project_map pm = {l_grade,l_size,r_grade,r_size,k_grade,k_size};

    sparse sparse_y = project_sparse_product(mvs,pm);
    blades blades_y = project_blades_product(gmvs,pm);
    dense dense_y = project_dense_product(dmvs,pm);

    print_all_types(sparse_y,dense_y,blades_y);

    free_blades(blades_y);
    free_sparse(sparse_y);
    free(dense_y.value);

    sparse_y = sparse_general_product(mvs,pm);
    blades_y = graded_general_product(gmvs,pm);
    dense_y = dense_general_product(dmvs,pm);

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
            printf("%s(%zu,%f)\n",s,i,y.value[i]);
        }
    }
}

void print_sparse(sparse y,char *s){
    for(size_t i = 0; i < y.size; i++ )
        printf("%s(%d,%f)\n",s,y.bitmap[i],y.value[i]);
}

void print_blades(blades y,char *s){
    for(size_t i = 0; i < y.size; i++){
        printf("%sGrade %d:\n",s, y.grade[i]);
        for(size_t j = 0; j < y.data[i].size; j++){
            printf("%s\t(%d,%f)\n",s,y.data[i].bitmap[j],y.data[i].value[j]);
        }
        printf("\n");
    }
}

dense project_dense_product(dense_multivectors mvs, project_map pm){
    dense y,proj_y;
    dense xl = dense_grade_project(mvs.a,pm.l,pm.l_size,mvs.dgm);
    dense xr = dense_grade_project(mvs.b,pm.r,pm.r_size,mvs.dgm);
    mvs.a = xl;
    mvs.b = xr;
    y = dense_product(mvs);
    proj_y = dense_grade_project(y,pm.k,pm.k_size,mvs.dgm);

    free(xl.value);
    free(xr.value);
    free(y.value);

    return proj_y;
}

sparse project_sparse_product(sparse_multivectors mvs, project_map pm){
    dense_grade_map dgm = mvs.dgm;
    sparse y,proj_y;
    sparse xl = sparse_grade_project(mvs.a,pm.l,pm.l_size,dgm);
    sparse xr = sparse_grade_project(mvs.b,pm.r,pm.r_size,dgm);
    mvs.a = xl;
    mvs.b = xr;
    y = sparse_product(mvs);
    proj_y = sparse_grade_project(y,pm.k,pm.k_size,dgm);

    free_sparse(xl);
    free_sparse(xr);
    free_sparse(y);

    return proj_y;
}

blades project_blades_product(graded_multivectors mvs, project_map pm){
    blades y,proj_y;
    blades xl = grade_sparse_grade_project(mvs.a,pm.l,pm.l_size,mvs.gm.size);
    blades xr = grade_sparse_grade_project(mvs.b,pm.r,pm.r_size,mvs.gm.size);
    mvs.a = xl;
    mvs.b = xr;
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

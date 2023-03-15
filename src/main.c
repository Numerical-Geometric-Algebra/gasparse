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
    /* test_parse(); */
    /* test_matrix_mult(); */
    test_general_einsum();
    /* test_einsum(); */
    /* test_parser_expression(); */

}

void test_parser_expression(void){
    /* char expression[] = "a"; */
    char expression[100];
    printf("expression: ");
    scanf("%s",expression);
    parse_expression(expression,strlen(expression));
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


sparse *gen_random_tensor(size_t size, int bitmap_max, size_t n_values){
    sparse *mvs = (sparse*)malloc(size*sizeof(sparse));
    for(size_t j = 0; j < size; j++){
        float *value = (float*)malloc(n_values*sizeof(float));
        int *bitmap = (int*)malloc(n_values*sizeof(int));
        for(size_t i = 0; i < n_values; i++){
            bitmap[i] = rand() % bitmap_max;
            value[i] = (float)rand()/(float)(RAND_MAX);
        }
        mvs[j].value = value;
        mvs[j].bitmap = bitmap;
        mvs[j].size = n_values;
    }
    return mvs;
}

sparse **gen_random_tensors(size_t **shapes, size_t *shape_size, size_t size, int vec_size, size_t n_values, size_t **data_size){
    int bitmap_max = 1 << vec_size;
    sparse **data = (sparse**)malloc(size*sizeof(sparse*));
    for(size_t i = 0; i < size; i++){
        int size_i = 1;
        for(size_t j = 0; j < shape_size[i]; j++){ // determine the size of the data
            size_i *= shapes[i][j];
        }
        (*data_size)[i] = size_i;
        data[i] = gen_random_tensor(size_i,bitmap_max,n_values);
    }
    return data;
}

void free_sparse_tensors(sparse **data, size_t *data_size, size_t size){
    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < data_size[i]; j++){
            free_sparse(data[i][j]);
        }free(data[i]);
    }free(data);
}

void free_graded_tensors(blades **data, size_t *data_size, size_t size){
    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < data_size[i]; j++){
            free_blades(data[i][j]);
        }free(data[i]);
    }free(data);
}

void free_graded_tensor_(blades *data, size_t data_size){
    for(size_t j = 0; j < data_size; j++){
        free_blades(data[j]);
    }
    free(data);
}

blades **sparse_to_graded_tensors(sparse **data, size_t *data_size, size_t size, grade_map gm){
    blades **graded_data = (blades**)malloc(size*sizeof(blades*));
    for(size_t i = 0; i < size; i++){
        graded_data[i] = (blades*)malloc(data_size[i]*sizeof(blades));
        for(size_t j = 0; j < data_size[i]; j++){
            graded_data[i][j] = sparse_to_graded(data[i][j],gm);
        }
    }
    return graded_data;
}

void test_general_einsum(void){
    unsigned int p = 4, q = 1, r = 1;
    float precision = 1e-12;
    size_t shape0[] = {2,3};
    size_t shape1[] = {3,2};
    size_t shape2[] = {3,3};

    size_t tensor_size = 3;
    size_t sparse_size = 1;
    size_t **shapes = (size_t**)malloc(tensor_size*sizeof(size_t*));
    size_t *data_size = (size_t*)malloc(tensor_size*sizeof(size_t*));

    shapes[0] = shape0;
    shapes[1] = shape1;
    shapes[2] = shape2;

    size_t shape_size[] = {2,2,2}; // just 2D tensors

    map m = cayley_table(p,q,r);
    grade_map gm = bitmap_grade_map(m.size);

    /* time_t t; */
    /* srand((unsigned) time(&t)); // set random seed */
    srand(75843);
    sparse **data = gen_random_tensors(shapes,shape_size,tensor_size,p+q+r,sparse_size,&data_size);
    blades **graded_data = sparse_to_graded_tensors(data,data_size,tensor_size,gm);

    size_t ndims[] = {2,2,2,2};
    /* char args[] = "ik,kj->ij"; */
    char args[] = "ik,kl,lj->ij";
    symbols s = parse_all(args,strlen(args),ndims,tensor_size+1);
    /* graded_tensor_multivectors tmvs = {graded_data,shapes,shape_size,data_size,tensor_size,m,gm,precision}; */

    printf("in left tensor:\n");
    for(size_t i = 0; i < data_size[0]; i++){
        printf("index %zu:\n",i);
        print_blades(graded_data[0][i],"\t");
    }

    printf("in right tensor:\n");
    for(size_t i = 0; i < data_size[1]; i++){
        printf("index %zu:\n",i);
        print_blades(graded_data[1][i],"\t");
    }

    tensor out_tensor;
    operator_functions opf =
        { graded_atomic_add__,
          graded_add_add__,
          graded_product__,
          graded_init__,
          graded_assign__,
          graded_free__};

    graded_extra gextra = {m,gm,precision};
    tensor_multivectors gtmvs =
        {(void**)graded_data,shapes,shape_size,data_size,tensor_size,sizeof(blades)};
    int error = main_einsum(gtmvs,(void*)&gextra,opf,s,&out_tensor);

    printf("error: %d\n",error);
    blades *graded_out_data = out_tensor.data;

    printf("out tensor:\n");
    for(size_t i = 0; i < out_tensor.data_size; i++){
        printf("index %zu:\n",i);
        print_blades(graded_out_data[i],"\t");
    }


    free_graded_tensor_(out_tensor.data,out_tensor.data_size);
    free(out_tensor.shapes);
    free_sparse_tensors(data,data_size,tensor_size); // free the data allocated to the input tensors
    free_graded_tensors(graded_data,data_size,tensor_size); // free the data allocated to the input converted tensors
    free(shapes);
    free(data_size);
    free_symbols(s);
    free_grade_map(gm);
    free_map(m);

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

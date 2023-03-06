#include "sparse.h"

// returns true if abs(v) < p
int comp_abs(float v,float p){
    float r = (v < 0) ? -v : v;
    return r < p;
}

sparse sparse_copy(sparse mv){
    sparse mv_copy;
    mv_copy.bitmap = (int*)malloc(mv.size*sizeof(int));
    mv_copy.value = (float*)malloc(mv.size*sizeof(float));
    mv_copy.size = mv.size;
    for(unsigned int i = 0; i < mv.size; i++){
        mv_copy.bitmap[i] = mv.bitmap[i];
        mv_copy.value[i] = mv.value[i];
    }
    return mv_copy;
}


sparse sparse_grade_project(
    sparse mv,
    unsigned int *project_grade,
    size_t project_size,
    dense_grade_map dgm){

    // computes a bool array where each index indicates the selected grade
    unsigned int *g = get_grade_bool(project_grade,project_size,dgm.grade_size);
    sparse projected_mv;

    unsigned int proj_size = 0;
    for(unsigned int i = 0; i < mv.size; i++)
        if(g[dgm.grade[mv.bitmap[i]]])
            proj_size++;

    projected_mv = initialize_sparse(proj_size--);

    // copies the values of the selected grades
    for(unsigned int i = 0; i < mv.size; i++){
        if(g[dgm.grade[mv.bitmap[i]]]){
            projected_mv.value[proj_size] = mv.value[i];
            projected_mv.bitmap[proj_size] = mv.bitmap[i];
            proj_size--;
            if(proj_size < 0)
                break;
        }
    }

    free(g);
    return projected_mv;
}

sparse sparse_general_product(sparse_multivectors mvs, project_map pm, dense_grade_map dgm){
    if(mvs.size < 2){
        // return error or something
    }
    unsigned int *ga = get_grade_bool(pm.l,pm.l_size,dgm.max_grade+1);
    unsigned int *gb = get_grade_bool(pm.r,pm.r_size,dgm.max_grade+1);
    unsigned int *gy = get_grade_bool(pm.k,pm.k_size,dgm.max_grade+1);
    sparse a = mvs.data[0];
    sparse b = mvs.data[1];
    map m = mvs.m;
    unsigned int m_size = mvs.m.size;

    int a_size = a.size;
    int b_size = b.size;

    // Allocate memory for a dense y
    sparse dense_y = initialize_sparse(m_size);
    sparse sparse_y;
    unsigned int sparse_size = 0;
    /* unsigned int k; */

    for(int i = 0; i < a_size; i++){
        if(!ga[dgm.grade[a.bitmap[i]]]) continue; // skip grade
        for(int j = 0; j < b_size; j++){
            if(!gb[dgm.grade[b.bitmap[j]]]) continue; // skip grade
            int sign = m.sign[a.bitmap[i]][b.bitmap[j]];
            // skip product if sign is null
            if(sign == 0)
                continue;
            unsigned int bitmap = m.bitmap[a.bitmap[i]][b.bitmap[j]];
            if(!gy[dgm.grade[bitmap]]) continue; // skip grade
            float value = a.value[i]*b.value[j];

            // write bitmap once to memory
            if(dense_y.bitmap[bitmap] == -1){
                dense_y.bitmap[bitmap] = bitmap;
                sparse_size++;// increment size of sparse
            }
            dense_y.value[bitmap] += value*sign;
        }
    }

    sparse_remove_small(dense_y,mvs.precision,&sparse_size);
    sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    free_sparse(dense_y);
    free(ga);
    free(gb);
    free(gy);
    return sparse_y;
}

// computes the geometric product between two sparse multivectors
sparse sparse_product(sparse_multivectors mvs) {
    if(mvs.size < 2){
        // return error or something
    }
    sparse a = mvs.data[0];
    sparse b = mvs.data[1];
    map m = mvs.m;
    unsigned int m_size = mvs.m.size;

    int a_size = a.size;
    int b_size = b.size;

    // Allocate memory for a dense y
    sparse dense_y = initialize_sparse(m_size);
    sparse sparse_y;
    unsigned int sparse_size = 0;
    /* unsigned int k; */

    for(int i = 0; i < a_size; i++){
        for(int j = 0; j < b_size; j++){
            int sign = m.sign[a.bitmap[i]][b.bitmap[j]];
            // skip product if sign is null
            if(sign == 0)
                continue;
            unsigned int bitmap = m.bitmap[a.bitmap[i]][b.bitmap[j]];
            float value = a.value[i]*b.value[j];

            // write bitmap once to memory
            if(dense_y.bitmap[bitmap] == -1){
                dense_y.bitmap[bitmap] = bitmap;
                sparse_size++;// increment size of sparse
            }
            dense_y.value[bitmap] += value*sign;
        }
    }

    sparse_remove_small(dense_y,mvs.precision,&sparse_size);
    sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    free_sparse(dense_y);
    return sparse_y;
}



sparse sparse_dense_to_sparse_sparse(sparse dense_y, unsigned int size){
    unsigned int k = 0;
    sparse sparse_y = initialize_sparse(size);
    for(unsigned int i = 0; i < dense_y.size; i++){
        if(k < size){
            if(dense_y.bitmap[i] != -1){
                sparse_y.bitmap[k] = dense_y.bitmap[i];
                sparse_y.value[k] = dense_y.value[i];
                k++;
            }
        }
    }
    sparse_y.size = size;
    return sparse_y;
}

void sparse_remove_small(sparse y, float precision,unsigned int *size){
     // Remove if value is too small
    for(unsigned int i = 0; i < y.size; i++){
        // Check if value was set
        if(y.bitmap[i] > -1){
            // Check if value is too small
            if(comp_abs(y.value[i],precision)){
                y.bitmap[i] = -1;
                (*size)--;
            }
        }
    }
}

sparse initialize_sparse(unsigned int size){
    sparse y;
    y.bitmap = (int*)malloc(size*sizeof(int));
    y.value = (float*)malloc(size*sizeof(float));
    y.size = size;
    for(unsigned int i = 0; i < size; i++){
        y.bitmap[i] = -1;
        y.value[i] = 0;
    }
    return y;
}

void free_sparse(sparse mv){
    free(mv.bitmap);
    free(mv.value);
}

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

sparse sparse_general_product(sparse_multivectors mvs, project_map pm){
    dense_grade_map dgm = mvs.dgm;
    sparse a = mvs.a;
    sparse b = mvs.b;
    map m = mvs.m;
    float precision = mvs.precision;
    return sparse_general_product_(a,b,m,pm,dgm,precision);
}


sparse sparse_general_product_(sparse a, sparse b, map m, project_map pm, dense_grade_map dgm, float precision){

    unsigned int *ga = get_grade_bool(pm.l,pm.l_size,dgm.max_grade+1);
    unsigned int *gb = get_grade_bool(pm.r,pm.r_size,dgm.max_grade+1);
    unsigned int *gy = get_grade_bool(pm.k,pm.k_size,dgm.max_grade+1);

    unsigned int m_size = m.size;

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

    sparse_remove_small(dense_y,precision,&sparse_size);
    sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    free_sparse(dense_y);
    free(ga);
    free(gb);
    free(gy);
    return sparse_y;
}

sparse sparse_scalar_multiply(float scalar, sparse b){
    sparse y = initialize_sparse(b.size);
    for(unsigned int i = 0; i < b.size; i++){
        y.value[i] = b.value[i]*scalar;
        y.bitmap[i] = b.bitmap[i];
    }
    return y;
}

sparse sparse_add_append(sparse a, sparse b){
    sparse y = initialize_sparse(a.size+b.size);
    for(size_t i = 0; i < a.size; i++){
        y.value[i] = a.value[i];
        y.bitmap[i] = a.bitmap[i];
    }

    for(size_t i = 0; i < b.size; i++){
        y.value[i+a.size] = b.value[i];
        y.bitmap[i+a.size] = b.bitmap[i];
    }

    return y;
}


// appends a bunch of multivectors together
sparse sparse_atomic_add_append(sparse **mv, size_t mv_size){
    size_t size_y = 0;
    for(size_t i =0; i < mv_size; i++)
        size_y += mv[i]->size;
    sparse y = initialize_sparse(size_y);

    for(size_t j = 0; j < mv[0]->size; j++){
        y.value[j] = mv[0]->value[j];
        y.bitmap[j] = mv[0]->bitmap[j];
    }

    size_t p = 0;
    for(size_t i = 1; i < mv_size; i++){
        p += mv[i-1]->size;
        for(size_t j = 0; j < mv[i]->size; j++){
            y.value[j+p] = mv[i]->value[j];
            y.bitmap[j+p] = mv[i]->bitmap[j];
        }
    }

    return y;
}


sparse sparse_add_add(sparse_multivectors mvs){
    sparse a = mvs.a;
    sparse b = mvs.b;
    float precision = mvs.precision;
    map m = mvs.m;
    return sparse_add_add_(a,b,m.size,precision);
}

sparse sparse_add_add_(sparse a, sparse b, unsigned int size, float precision){
    sparse dense_y = initialize_sparse(size);
    sparse sparse_y;
    unsigned int sparse_size = 0;
    for(size_t i = 0; i < a.size; i++){
        if(dense_y.bitmap[a.bitmap[i]]==-1){
            dense_y.bitmap[a.bitmap[i]] = a.bitmap[i];
            sparse_size++;
        }
        dense_y.value[a.bitmap[i]] += a.value[i];
    }
    for(size_t i = 0; i < b.size; i++){
        if(dense_y.bitmap[b.bitmap[i]]==-1){
            dense_y.bitmap[b.bitmap[i]] = b.bitmap[i];
            sparse_size++;
        }
        dense_y.value[b.bitmap[i]] += b.value[i];
    }
    sparse_remove_small(dense_y,precision,&sparse_size);
    sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    free_sparse(dense_y);
    return sparse_y;
}


// adds a bunch of multivectors together
sparse sparse_atomic_add_add_(sparse **mv, size_t mv_size, size_t size, float precision){
    sparse dense_y = initialize_sparse(size);
    sparse sparse_y;
    /* printf("dense_y size: %d",dense_y.size); */
    unsigned int sparse_size = 0;
    for(size_t i = 0; i < mv_size; i++){
        for(size_t j = 0; j < mv[i]->size; j++){
            int bitmap = mv[i]->bitmap[j];
            if(dense_y.bitmap[bitmap] == -1){
                dense_y.bitmap[bitmap] = bitmap;
                sparse_size++;
            }
            dense_y.value[bitmap] += mv[i]->value[j];
        }
    }

    sparse_remove_small(dense_y,precision,&sparse_size);
    sparse_y = sparse_dense_to_sparse_sparse(dense_y,sparse_size);

    free_sparse(dense_y);
    return sparse_y;
}



// computes the geometric product between two sparse multivectors
sparse sparse_product(sparse_multivectors mvs) {
    sparse a = mvs.a;
    sparse b = mvs.b;
    map m = mvs.m;
    float precision = mvs.precision;
    return sparse_product_(a,b,m,precision);
}


sparse sparse_product_(sparse a, sparse b, map m, float precision) {
    unsigned int m_size = m.size;

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

    sparse_remove_small(dense_y,precision,&sparse_size);
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

#include "dense.h"

dense dense_grade_project(
    dense mv,
    unsigned int *project_grade,
    size_t project_size,
    dense_grade_map dgm){

    // computes a bool array which each index indicates the selected grade
    unsigned int *g = get_grade_bool(project_grade,project_size,dgm.grade_size);
    dense projected_mv = initialize_dense(mv.size);

    // copies the values of the selected grades
    for(unsigned int i = 0; i < mv.size; i++)
        if(g[dgm.grade[i]])
            projected_mv.value[i] = mv.value[i];

    free(g);
    return projected_mv;
}

dense dense_general_product(dense_multivectors mvs, project_map pm){
    dense a = mvs.a;
    dense b = mvs.b;
    map m = mvs.m;
    dense_grade_map dgm = mvs.dgm;

    return dense_general_product_(a,b,m,pm,dgm);
}

dense dense_scalar_multiply(float scalar, dense b){
    dense y = initialize_dense(b.size);
    for(unsigned int i = 0; i < b.size; i++){
        y.value[i] = b.value[i]*scalar;
    }
    return y;
}

dense dense_add(dense a, dense b){
    dense y = initialize_dense(a.size);
    for(size_t i = 0; i < a.size; i++){
        y.value[i] += a.value[i] + b.value[i];
    }
    return y;
}

dense dense_atomic_add(dense *mv, size_t size){
    dense y = initialize_dense(mv[0].size);
    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < mv[i].size; j++){
            y.value[j] += mv[i].value[j];
        }
    }
    return y;
}

dense dense_general_product_(dense a, dense b, map m, project_map pm, dense_grade_map dgm){

    unsigned int *ga = get_grade_bool(pm.l,pm.l_size,dgm.max_grade+1);
    unsigned int *gb = get_grade_bool(pm.r,pm.r_size,dgm.max_grade+1);
    unsigned int *gy = get_grade_bool(pm.k,pm.k_size,dgm.max_grade+1);

    unsigned int m_size = m.size;

    int a_size = a.size;
    int b_size = b.size;

    // Allocate memory for a dense y
    dense dense_y = initialize_dense(m_size);

    for(int i = 0; i < a_size; i++){
        if(!ga[dgm.grade[i]]) continue;
        for(int j = 0; j < b_size; j++){
            if(!gb[dgm.grade[j]]) continue;
            int sign = m.sign[i][j];
            // skip product if sign is null
            if(sign == 0)
                continue;
            unsigned int bitmap = m.bitmap[i][j];
            if(!gy[dgm.grade[bitmap]]) continue;
            float value = a.value[i]*b.value[j];

            dense_y.value[bitmap] += value*sign;
        }
    }

    free(ga);
    free(gb);
    free(gy);
    return dense_y;
}



// computes the geometric product between two dense multivectors
dense dense_product(dense_multivectors mvs) {
    dense a = mvs.a;
    dense b = mvs.b;
    map m = mvs.m;
    return dense_product_(a,b,m);
}


dense dense_product_(dense a,dense b, map m) {
    unsigned int m_size = m.size;

    int a_size = a.size;
    int b_size = b.size;

    // Allocate memory for a dense y
    dense dense_y = initialize_dense(m_size);

    for(int i = 0; i < a_size; i++){
        for(int j = 0; j < b_size; j++){
            int sign = m.sign[i][j];
            // skip product if sign is null
            if(sign == 0)
                continue;
            unsigned int bitmap = m.bitmap[i][j];
            float value = a.value[i]*b.value[j];

            dense_y.value[bitmap] += value*sign;
        }
    }

    return dense_y;
}


dense sparse_to_dense(sparse mv, unsigned int size){
    dense dense_mv = initialize_dense(size);

    for(unsigned int i = 0; i < mv.size; i++)
        dense_mv.value[mv.bitmap[i]] = mv.value[i];
    return dense_mv;
}

// computes the geometric product using the dense representation
dense inverse_dense_product(dense_multivectors mvs){
    dense a = mvs.a;
    dense b = mvs.b;
    map m = mvs.m;
    unsigned int m_size = mvs.m.size;

    // Allocate memory for a dense y
    dense dense_y = initialize_dense(m_size);

    // iterate over each basis blade of the output multivector y
    for(size_t i = 0; i < m_size; i++){
        for(size_t j = 0; j < m_size; j++){
            int sign = m.sign[i][j];
            // skip product if sign is null
            if(sign == 0)
                continue;
            unsigned int bitmap = m.bitmap[i][j];
            float value = a.value[j]*b.value[bitmap];
            dense_y.value[i] += value*sign;
        }
    }

    return dense_y;
}


dense initialize_dense(unsigned int size){
    dense y;
    y.value = (float*)malloc(size*sizeof(float));
    y.size = size;
    for(unsigned int i = 0; i < size; i++){
        y.value[i] = 0;
    }
    return y;
}


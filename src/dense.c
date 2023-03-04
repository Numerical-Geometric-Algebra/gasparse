#include "dense.h"

dense dense_grade_project(
    dense mv,
    unsigned int *grade,
    unsigned int *project_grade,
    size_t grade_size,
    size_t project_size){

    // computes a bool array which each index indicates the selected grade
    unsigned int *g = get_grade_bool(project_grade,project_size,grade_size);
    dense projected_mv = initialize_dense(mv.size);

    // copies the values of the selected grades
    for(unsigned int i = 0; i < mv.size; i++)
        if(g[grade[i]])
            projected_mv.value[i] = mv.value[i];

    free(g);
    return projected_mv;
}


// computes the geometric product between two dense multivectors
dense dense_product(dense_multivectors mvs) {
    if(mvs.size < 2){
        // return error or something
    }
    dense a = mvs.data[0];
    dense b = mvs.data[1];
    map m = mvs.m;
    unsigned int m_size = mvs.m.size;

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
    if(mvs.size < 2){
        // return error or something
    }
    dense a = mvs.data[0];
    dense b = mvs.data[1];
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


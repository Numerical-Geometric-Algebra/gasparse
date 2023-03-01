#include "sparse.h"

// returns true if abs(v) < p
int comp_abs(float v,float p){
    float r = (v < 0) ? -v : v;
    return r < p;
}


sparse geometric_product(sparse_multivectors mvs) {
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
    unsigned int k;

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

    // Remove if value is too small
    for(unsigned int i = 0; i < m_size; i++){
        // Check if value was set
        if(dense_y.bitmap[i] > -1){
            // Check if value is too small
            if(comp_abs(dense_y.value[i],mvs.precision)){
                dense_y.bitmap[i] = -1;
                sparse_size--;
            }
        }
    }

    k = 0;
    sparse_y = initialize_sparse(sparse_size);
    for(unsigned int i = 0; i < m_size; i++){
        if(k < sparse_size){
            if(dense_y.bitmap[i] != -1){
                sparse_y.bitmap[k] = dense_y.bitmap[i];
                sparse_y.value[k] = dense_y.value[i];
                k++;
            }
        }
    }
    sparse_y.size = sparse_size;
    free(dense_y.bitmap);
    free(dense_y.value);

    return sparse_y;
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


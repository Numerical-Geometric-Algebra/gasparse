#include "grade-sparse.h"

/* sparse_multivector grade_selection(sparse_multivector a, unsigned int grade){ */
/*     if(grade < a.graded_mvector.grades_size){ */
/*         a.comps = a.graded_mvector.graded_comps[grade]->comps; */
/*         a.graded_mvector.graded_comps = NULL; */
/*     } */
/*     return a; */
/* } */


blades *graded_product(graded_multivectors mvs){
    if(mvs.size < 2){
        // return error
    }
    blades a = mvs.data[0];
    blades b = mvs.data[1];
    graded_map m = mvs.m;
    size_t m_size = mvs.size;
    grade_info ginfo = mvs.ginfo;

    blades dense_y = initialize_blades(ginfo);
    unsigned int *grade_size = (unsigned int *)malloc(ginfo.max_grade*sizeof(unsigned int));

    for(size_t i = 0; i < a.size; i++){
        for(size_t j = 0; j < b.size; j++){ // iterate over grades
            sparse mv_a = a.data[i];
            sparse mv_b = b.data[j];
            for(size_t k = 0; k < a.data[i].size; k++){
                for(size_t l = 0; l < b.data[j].size; l++){ // iterate over basis vectors
                    int sign = m.sign[mv_a.bitmap[k]][mv_b.bitmap[l]];
                    if(sign == 0)
                        continue;
                    unsigned int grade = m.grade[mv_a.bitmap[k]][mv_b.bitmap[l]];
                    unsigned int bitmap = m.bitmap[mv_a.bitmap[k]][mv_b.bitmap[l]];
                    float value = mv_a.value[k]*mv_b.value[l];

                    // write bitmap once to memory
                    if(dense_y.data[grade].bitmap[bitmap-grade] == -1){
                        dense_y.bitmap[bitmap] = bitmap;
                        sparse_size++;// increment size of sparse
                    }
                    dense_y.value[bitmap] += value*sign;
                }
            }
        }
    }
}


set_dense_value(blades y, int bitmap, int grade, grade_info ginfo){

}

blades initialize_blades(grade_info ginfo){
    blades y;
    y.data = (sparse*)malloc(ginfo.max_grade*sizeof(sparse));
    y.grades = (unsigned int*)malloc(ginfo.max_grade*sizeof(unsigned int));
    y.size = ginfo.max_grade;
    for(size_t i = 0; i < ginfo.max_grade; i++){
        y.data[i] = initialize_sparse(ginfo.grade_size[i]);
        y.grades[i] = i;
    }
    return y;
}


sparse_multivector geometric_product(multivectors data) {
    if(data.mvs_size < 2){
        // return error or something
    }
    sparse_multivector a = data.mvs[0];
    sparse_multivector b = data.mvs[1];
    algebra_map **map = data.map;
    unsigned int algebra_size = data.algebra_size;
    unsigned int grades_size = data.grades_size;

    int a_size = a.comps.size;
    int b_size = b.comps.size;

    sparse_multivector y;
    components comps;
    graded_multivector g_mv = initialize_graded_mv(grades_size);

    // allocate memory as if a dense multivector is to be computed
    comps.sing_comp = initialize_comps(algebra_size);
    y.comps.size = 0;
    y.graded_mvector.size = 0;

    for(int i = 0; i < a_size; i++){
        unsigned int n = a.comps.sing_comp[i]->index;

        for(int j = 0; j < b_size; j++){
            unsigned int m = b.comps.sing_comp[j]->index;
            algebra_map map_nm = map[n][m];

            if(map_nm.sign != 0){
                // Allocate memory if not already
                if(comps.sing_comp[map_nm.index] == NULL){
                    y.comps.size++;
                    comps.sing_comp[map_nm.index] = (single_component*)malloc(sizeof(single_component));
                    comps.sing_comp[map_nm.index]->index = map_nm.index;
                    comps.sing_comp[map_nm.index]->grade = map_nm.grade;
                    graded_components *comps = g_mv.graded_comps[map_nm.grade];
                    if(comps == NULL){
                        y.graded_mvector.size++;
                        comps = (graded_components*)malloc(sizeof(graded_components));
                        comps->grade = map_nm.grade;
                        comps->comps.size = 1;
                        g_mv.graded_comps[map_nm.grade] = comps;
                    }else{
                        comps->comps.size++;
                    }
                }
                // Compute the product
                comps.sing_comp[map_nm.index]->value
                    += map_nm.sign*a.comps.sing_comp[i]->value*b.comps.sing_comp[j]->value;
            }
        }
    }

    y.comps.sing_comp = populate_comps(comps.sing_comp,y.comps.size,algebra_size);
    y.graded_mvector = populate_graded_mv(g_mv,comps.sing_comp,y.graded_mvector.size,algebra_size);
    free(comps.sing_comp);
    free(g_mv.graded_comps);
    return y;
}

graded_multivector populate_graded_mv(
    graded_multivector dense_mv,
    single_component **comp,
    unsigned int size,
    unsigned int algebra_size){

    graded_multivector sparse_mv;
    sparse_mv.graded_comps = (graded_components**)malloc(size*sizeof(graded_components*));

    unsigned int j = 0;
    unsigned int* k = (unsigned int*)malloc(dense_mv.size*sizeof(unsigned int));
    for(unsigned int i = 0; i < dense_mv.size; i++){
        if(dense_mv.graded_comps[i] != NULL){
            if(j < size){
                sparse_mv.graded_comps[j] = dense_mv.graded_comps[i];
                j++;
            }
            k[i] = dense_mv.graded_comps[i]->comps.size;
        }else{
            k[i] = 0;
        }
    }
    // writes the components for each grade
    for(unsigned int i = 0; i < algebra_size; i++){
        if(comp[i] != NULL){
            dense_mv.graded_comps[comp[i]->grade]->comps.sing_comp[--k[comp[i]->grade]] = comp[i];
        }
    }

    return sparse_mv;
}


graded_multivector initialize_graded_mv(unsigned int grades_size){
    graded_multivector g_mv;
    g_mv.size = grades_size;
    g_mv.graded_comps = (graded_components**)malloc(grades_size*sizeof(graded_components*));
    for(unsigned int i = 0; i < grades_size; i++){
        g_mv.graded_comps[i] = NULL;
    }
    return g_mv;
}


// Initializes all components to NULL
single_component **initialize_comps(int algebra_size){
    single_component **comps;
    comps = (single_component**)malloc(algebra_size*sizeof(single_component*));
    for(int i = 0; i < algebra_size; i++)
        comps[i] = NULL;
    return comps;
}

/* Populates sparse multivectors with only the non-zero components */
single_component **populate_comps(single_component **full_comps, int size, int algebra_size){
    single_component **comps = initialize_comps(size);

    int j = 0;
    for(int i = 0; i < algebra_size; i++){
        if(full_comps[i] != NULL){
            if(j < size){
                comps[j] = full_comps[i];
                j++;
            }
        }
    }
    return comps;
}

algebra_map geo_prod(unsigned int a,unsigned int b,int *m, size_t m_size){
    unsigned int bitmap = a ^ b;
    int sign = sign_reorder(a,b)*sign_metric(a,b,m,m_size);
    unsigned int grade = grade_compute(bitmap);
    algebra_map map = {sign,bitmap,grade};
    return map;
}


algebra_map **cayley_table(int *m,size_t m_size,size_t *mv_size){
    size_t n = *mv_size = (size_t)pow((double)2,(double)m_size);

    algebra_map **map = (algebra_map**)malloc(n*sizeof(algebra_map*));
    for(size_t i = 0; i < n; i++ ){
        map[i] = malloc(n*sizeof(algebra_map));
        for(size_t j = 0; j < n; j++){
            map[i][j] = geo_prod(i,j,m,m_size);
        }
    }
    return map;
}

#include "grade_sparse.h"

blades grade_sparse_grade_project(blades mv,
    unsigned int *project_grade,
    size_t grade_size,
    size_t project_size){

    // computes a bool array which each index indicates the selected grade
    unsigned int *g = get_grade_bool(project_grade,project_size,grade_size);
    blades projected_mv;

    size_t projected_size = 0;

    for(unsigned int i = 0; i < mv.size; i++)
        if(g[mv.grade[i]])
            projected_size++;

    projected_mv = initialize_blades_empty(projected_size--);

    for(unsigned int i = 0; i < mv.size; i++){
        if(g[mv.grade[i]]){
            projected_mv.data[projected_size] = sparse_copy(mv.data[i]);
            projected_mv.grade[projected_size--] = mv.grade[i];
        }
    }
    free(g);
    return projected_mv;
}


// takes a sparse representation of a multivector and converts it to grade sparse
blades sparse_to_graded(sparse mv, grade_map gm){
    blades sparse_mv;
    blades dense_mv;
    unsigned int *grade_size = initialize_grade_size(gm);
    size_t n_grades = 0;
    size_t max_grade = 0;
    // determine how many components for each grade
    for(size_t i = 0; i < mv.size; i++){
        unsigned int grade = gm.grade[mv.bitmap[i]];
        grade_size[grade]++;
        if(grade > max_grade)
            max_grade = grade;
    }
    // counts the number of different grades
    for(size_t i = 0; i <= max_grade; i++)
        if(grade_size[i] > 0)
            n_grades++;

    dense_mv = initialize_blades_empty(max_grade+1);
    sparse_mv = initialize_blades_empty(n_grades);

    // convert sparse to grade dense
    for(size_t i = 0; i < mv.size; i++){
        unsigned int grade = gm.grade[mv.bitmap[i]];
        if(dense_mv.data[grade].bitmap == NULL)
            dense_mv.data[grade] = initialize_sparse(grade_size[grade]);

        grade_size[grade]--;
        dense_mv.data[grade].bitmap[grade_size[grade]] = mv.bitmap[i];
        dense_mv.data[grade].value[grade_size[grade]] = mv.value[i];
    }

    n_grades = 0;
    for(size_t i = 0; i <= max_grade; i++){
        if(dense_mv.data[i].bitmap != NULL){
            sparse_mv.data[n_grades] = dense_mv.data[i];
            sparse_mv.grade[n_grades] = dense_mv.grade[i];
            n_grades++;
        }
    }
    free(dense_mv.data);
    free(dense_mv.grade);
    free(grade_size);
    return sparse_mv;
}


blades graded_product(graded_multivectors mvs){
    if(mvs.size < 2){
        // return error
    }
    blades a = mvs.data[0];
    blades b = mvs.data[1];
    map m = mvs.m;
    grade_map gm = mvs.gm;
    blades dense_y = initialize_blades_empty(gm.max_grade + 1);
    blades sparse_y;
    unsigned int n_grades = 0;
    unsigned int *grade_size = initialize_grade_size(gm);

    // iterate over grades
    for(size_t i = 0; i < a.size; i++){
        for(size_t j = 0; j < b.size; j++){
            sparse mv_a = a.data[i];
            sparse mv_b = b.data[j];
            // iterate over basis vectors
            for(size_t k = 0; k < a.data[i].size; k++){
                for(size_t l = 0; l < b.data[j].size; l++){
                    int sign = m.sign[mv_a.bitmap[k]][mv_b.bitmap[l]];
                    if(sign == 0) continue;
                    unsigned int bitmap = m.bitmap[mv_a.bitmap[k]][mv_b.bitmap[l]];
                    unsigned int grade = gm.grade[bitmap];
                    unsigned int position = gm.position[bitmap];
                    float value = mv_a.value[k]*mv_b.value[l];

                    if(dense_y.data[grade].bitmap == NULL){
                        dense_y.data[grade] = initialize_sparse(gm.grade_size[grade]);
                        n_grades++;
                    }

                    // write bitmap once to memory
                    if(dense_y.data[grade].bitmap[position] == -1){
                        dense_y.data[grade].bitmap[position] = bitmap;
                        dense_y.data[grade].value[position] = 0;
                        grade_size[grade]++;
                    }
                    // compute the geometric product
                    dense_y.data[grade].value[position] += sign*value;
                }
            }
        }
    }

    // Remove if value is too small
    for(size_t i = 0; i <= gm.max_grade; i++){
        // Check if grade was set
        if(dense_y.data[i].bitmap != NULL){
            for(size_t j = 0; j < gm.grade_size[i]; j++){
                // Check if value was set
                if(dense_y.data[i].bitmap[j] != -1){
                    // compare value with precision
                    if(comp_abs(dense_y.data[i].value[j],mvs.precision)){
                        dense_y.data[i].bitmap[j] = -1;
                        grade_size[i]--;

                        // remove grade if all elements have small value
                        if(grade_size[i] <= 0){
                            free(dense_y.data[i].bitmap);
                            free(dense_y.data[i].value);
                            dense_y.data[i].bitmap = NULL;
                            dense_y.data[i].value = NULL;
                            break;
                        }
                    }
                }
            }
        }
    }

    sparse_y = initialize_blades(grade_size,gm.max_grade + 1);
    n_grades = 0;
    for(size_t i = 0; i <= gm.max_grade; i++){
        if(dense_y.data[i].bitmap != NULL){
            sparse_y.grade[n_grades] = dense_y.grade[i];
            for(size_t j = 0; j < gm.grade_size[i]; j++){
                if(dense_y.data[i].bitmap[j] != -1){
                    grade_size[i]--;
                    int *bitmap = sparse_y.data[n_grades].bitmap;
                    float *value = sparse_y.data[n_grades].value;
                    bitmap[grade_size[i]] = dense_y.data[i].bitmap[j];
                    value[grade_size[i]] = dense_y.data[i].value[j];
                }
            }
            n_grades++;
        }
    }
    free(grade_size);
    free_blades(dense_y);

    return sparse_y;
}


unsigned int *initialize_grade_size(grade_map gm){
    unsigned int *grade_size = (unsigned int *)malloc((gm.max_grade+1)*sizeof(unsigned int));
    for(size_t i = 0; i <= gm.max_grade; i++)
        grade_size[i] = 0;
    return grade_size;
}

blades initialize_blades_empty(size_t n_grades){// allocate the necessary memory for each grade
    blades y;
    y.data = (sparse*)malloc(n_grades*sizeof(sparse));
    y.grade = (unsigned int*)malloc(n_grades*sizeof(unsigned int));
    y.size = n_grades;
    for(size_t i = 0; i < n_grades; i++){
        y.data[i].bitmap = NULL;
        y.data[i].value = NULL;
        y.grade[i] = i;
    }
    return y;
}

blades initialize_blades(unsigned int *grade_size, size_t n_grades){// allocate the necessary memory for each grade
    blades y;
    size_t m_grades = 0;
    // determine the total number of nonempty grades
    for(size_t i = 0; i < n_grades; i++)
        if(grade_size[i] > 0)
            m_grades++;

    y.data = (sparse*)malloc(m_grades*sizeof(sparse));
    y.grade = (unsigned int*)malloc(m_grades*sizeof(unsigned int));
    y.size = m_grades;
    unsigned int j = 0;
    for(size_t i = 0; i < n_grades; i++){
        if(grade_size[i] > 0 ){
            y.data[j] = initialize_sparse(grade_size[i]);
            y.grade[j] = -1;
            j++;
        }
    }
    return y;
}


void free_blades(blades y){
    for(size_t i = 0; i < y.size; i++){
        if(y.data[i].bitmap != NULL)
            free(y.data[i].bitmap);

        if(y.data[i].value != NULL)
            free(y.data[i].value);
    }
    free(y.data);
    free(y.grade);
}

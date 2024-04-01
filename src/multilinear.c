#include <Python.h>
#include "types.h"
#include "multilinear.h"
#include "multivector_types.h"
#include "common.h"

SparseMultivector sparse_remove_relative_small(SparseMultivector x, ga_float percentage){
    ga_float x_max = 0;
    Py_ssize_t size = 0;
    for(Py_ssize_t i = 0; i < x.size; i++){
        if(x.bitmap[i] != -1) size++;
        if(x_max < fabs(x.value[i]))
            x_max = fabs(x.value[i]);
    }

    // Compare with the maximum
    for(Py_ssize_t i = 0; i < x.size; i++){
        if(fabs(x.value[i]) < x_max*percentage){
            x.bitmap[i] = -1;
            size--;
        }
    }

    SparseMultivector sparse = init_sparse_empty(size);
    if(sparse.size == -1)
        return sparse;

    Py_ssize_t j = 0;
    for(Py_ssize_t i = 0; i < x.size; i++){
        if(x.bitmap[i] != -1){
            sparse.bitmap[j] = x.bitmap[i];
            sparse.value[j] = x.value[i];
            j++;
        }
    }
    return sparse;

}

static int ternary_sparse_product(void *out, void *data0, void *data1, void *data2, PyAlgebraObject *ga, ProductType ptype){
    
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse2 = (SparseMultivector*)data2;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign = m.sign[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(!sign) continue;
            bitmap = m.bitmap[sparse0->bitmap[i]][sparse1->bitmap[j]];
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0->value[i]*sparse1->value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    *sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return 1;
}


// Takes the sandwich product of data0 with data1: mult(data1,data0,data1)
static int ternary_sparse_sandwich_product(void *out, void *data0, void *data1, PyAlgebraObject *ga, ProductType ptype){
    
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse = (SparseMultivector*)out;
    CliffordMap m = ga->product[ptype];
    
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return 0;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap0, bitmap1;
    int sign0, sign1;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        for(Py_ssize_t j = 0; j < sparse1->size; j++){
            sign0 = m.sign[sparse1->bitmap[j]][sparse0->bitmap[i]];
            if(!sign0) continue;
            bitmap0 = m.bitmap[sparse1->bitmap[j]][sparse0->bitmap[i]];
            for(Py_ssize_t k = 0; k < sparse1->size; k++){
                sign1 = m.sign[bitmap0][sparse1->bitmap[j]];
                if(!sign1) continue;
                bitmap1 = m.bitmap[bitmap0][sparse1->bitmap[j]];
                if(dense.bitmap[bitmap1] == -1) dense.bitmap[bitmap1] = bitmap1, size++;
                dense.value[bitmap1] += sparse0->value[i]*sparse1->value[j]*sign1;
            }
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    *sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return 1;
}

// the reverse argument reverses the last element data2
static int specialized_ternary_blades_product0(void *out, void *data0, void *data1, void *data2, ga_float scalar, int *grades_out, Py_ssize_t gsize, PyMultilinearObject *ml, ProductType ptype1, ProductType ptype2, int reverse){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades1 = (BladesMultivector*)data1;
    BladesMultivector *blades2 = (BladesMultivector*)data2;
    BladesMultivector *blades = (BladesMultivector*)out;
    PyAlgebraObject *ga = ml->GA;
    GradeMap gm = ga->gm;
    GradeTable gt = ga->gt;

    // Specialized ternary product when the multivectors are of unique grade
    if(blades0->size == 1 && blades1->size == 1 && blades2->size == 1 && gsize == 1){
        SparseMultivector sparse = init_sparse_empty(gm.grade_size[*grades_out]); // Initialize multivectors of unique grade
        if(sparse.size == -1) return 0;
        
        Py_ssize_t grade0 = *blades0->grade;
        Py_ssize_t grade1 = *blades1->grade;
        Py_ssize_t grade2 = *blades2->grade;

        int gsign;
        if(reverse) gsign = (*blades2->grade & 2) ? -1 : 1; // Sign for reversing the blades
        else gsign = 1;
        
        PositionMap tmap = ml->ternary_product[ptype1][ptype2][*blades0->grade][*blades1->grade][*blades2->grade];
        
        SparseMultivector sparse0 = *blades0->data;
        SparseMultivector sparse1 = *blades1->data;
        SparseMultivector sparse2 = *blades2->data;
                    
        for(Py_ssize_t i0 = 0; i0 < sparse0.size; i0++){
            Py_ssize_t position0 = gm.position[sparse0.bitmap[i0]];
            MapBase **tmap0 = tmap[position0];
            for(Py_ssize_t i1 = 0; i1 < sparse1.size; i1++){
                ga_float value0 = sparse0.value[i0]*sparse1.value[i1];
                Py_ssize_t position1 = gm.position[sparse1.bitmap[i1]];
                MapBase *tmap1 = tmap0[position1];
                for(Py_ssize_t i2 = 0; i2 < sparse2.size; i2++){
                    Py_ssize_t position2 = gm.position[sparse2.bitmap[i2]];
                    MapBase tmap2 = tmap1[position2];
                    if(*grades_out != tmap2.grade || !tmap2.sign) continue;
                    sparse.bitmap[tmap2.position] = gt.bitmaps[tmap2.grade][tmap2.position];
                    sparse.value[tmap2.position] += scalar*gsign*tmap2.sign*value0*sparse2.value[i2];
                }
            }
        }

        blades->data = (SparseMultivector*)PyMem_RawMalloc(sizeof(SparseMultivector));
        blades->grade = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
        blades->size = 1;

        *blades->data = sparse_remove_relative_small(sparse,ga->precision);

        return 1;
    }


    Py_ssize_t *g = get_grade_bool(grades_out, gsize, gm.max_grade + 1);
    SparseMultivector *sparse = (SparseMultivector*)PyMem_RawMalloc((gm.max_grade + 1)*sizeof(SparseMultivector));
    for(Py_ssize_t i = 0; i < gm.max_grade + 1; i++){
        sparse[i] = init_sparse_empty(gm.grade_size[i]); // Initialize multivectors of unique grade
        if(sparse[i].size == -1) return 0;
    }
    int gsign;

    for(Py_ssize_t j0 = 0; j0 < blades0->size; j0++){
        SparseMultivector sparse0 = blades0->data[j0];
        for(Py_ssize_t j1 = 0; j1 < blades1->size; j1++){
            SparseMultivector sparse1 = blades1->data[j1];
            for(Py_ssize_t j2 = 0; j2 < blades2->size; j2++){
                SparseMultivector sparse2 = blades2->data[j2];
                PositionMap tmap = ml->ternary_product[ptype1][ptype2][blades0->grade[j0]][blades1->grade[j1]][blades2->grade[j2]];
                if(reverse) gsign = (blades2->grade[j2] & 2) ? -1 : 1; // Sign for reversing the multivectors
                else gsign = 1;
                for(Py_ssize_t i0 = 0; i0 < sparse0.size; i0++){
                    Py_ssize_t position0 = gm.position[sparse0.bitmap[i0]];
                    MapBase **tmap0 = tmap[position0];
                    for(Py_ssize_t i1 = 0; i1 < sparse1.size; i1++){
                        ga_float value0 = sparse0.value[i0]*sparse1.value[i1];
                        Py_ssize_t position1 = gm.position[sparse1.bitmap[i1]];
                        MapBase *tmap1 = tmap0[position1];
                        for(Py_ssize_t i2 = 0; i2 < sparse2.size; i2++){
                            Py_ssize_t position2 = gm.position[sparse2.bitmap[i2]];
                            MapBase tmap2 = tmap1[position2];
                            if(!g[tmap2.grade] || !tmap2.sign) continue;
                            sparse[tmap2.grade].bitmap[tmap2.position] = gt.bitmaps[tmap2.grade][tmap2.position];
                            sparse[tmap2.grade].value[tmap2.position] += scalar*gsign*tmap2.sign*value0*sparse2.value[i2];
                        }
                    }
                }
            }
        }
    }

    blades->data = (SparseMultivector*)PyMem_RawMalloc(gsize*sizeof(SparseMultivector));
    blades->grade = (Py_ssize_t*)PyMem_RawMalloc(gsize*sizeof(Py_ssize_t));
    blades->size = gsize;

    Py_ssize_t j = 0;
    for(Py_ssize_t i = 0; i < gm.max_grade + 1; i++){
        if(g[i]){
            SparseMultivector temp = sparse_remove_relative_small(sparse[i],ga->precision);
            
            blades->data[j].bitmap = temp.bitmap;
            blades->data[j].value = temp.value;
            blades->grade[j] = i;
            
        }else{
            PyMem_RawFree(sparse[i].bitmap);
            PyMem_RawFree(sparse[i].value);
        }
    }
    PyMem_RawFree(sparse);

    return 1;
}


static int specialized_ternary_gradeddense_product0(void *out, void *data0, void *data1, void *data2, ga_float scalar, int *grades_out, Py_ssize_t gsize, PyMultilinearObject *ml, ProductType ptype1, ProductType ptype2, int flags){
    GrdDenseMv *mv0 = (GrdDenseMv*)data0;
    GrdDenseMv *mv1 = (GrdDenseMv*)data1;
    GrdDenseMv *mv2 = (GrdDenseMv*)data2;
    GrdDenseMv *mv = (GrdDenseMv*)out;

    PyAlgebraObject *ga = ml->GA;
    GradeMap gm = ga->gm;
    GradeTable gt = ga->gt;
    Py_ssize_t *g = get_grade_bool(grades_out, gsize, gm.max_grade + 1);
    SparseTernaryMap map = ml->sparse_ternary_product[ptype1][ptype2];
    GrdDenseMv dense = {.values = (ga_float**)PyMem_RawMalloc(gsize)};
    SparseMap smap;

    // Specialized when all multivectors are of unique grade
    if(gsize == 1 && mv0->size == 1 && mv1->size == 1 && mv2->size == 1){
        ga_float *values_out = (ga_float*)PyMem_RawMalloc(gm.grade_size[*grades_out]);

        *dense.values = values_out;
        *dense.grades = *grades_out;

        ga_float *values0 = *mv0->values;
        ga_float *values1 = *mv1->values;
        ga_float *values2 = *mv2->values;
        
        int gsign = 1;
        if(flags & OpFlag_reverse0) gsign *= (*mv0->grades & 2) ? -1 : 1; // Sign for reversing the blades
        if(flags & OpFlag_reverse1) gsign *= (*mv1->grades & 2) ? -1 : 1; // Sign for reversing the blades
        if(flags & OpFlag_reverse2) gsign *= (*mv2->grades & 2) ? -1 : 1; // Sign for reversing the blades

        smap = map[*mv0->grades][*mv1->grades][*mv2->grades][*dense.grades];
                    
        for(Py_ssize_t j4 = 0; j4 < smap.size; j4++){
            SparseMapBase base = smap.base[j4];
            values_out[base.position[3]] += gsign*scalar*values0[base.position[0]]*values1[base.position[1]]*values2[base.position[2]]*base.sign;
        }
        *mv = dense;
        return 1;
    }

    // Generalized for mixed grade multivectors
    int gsign;
    for(Py_ssize_t j3 = 0; j3 < gsize; j3++){
        ga_float *values_out = (ga_float*)PyMem_RawMalloc(gm.grade_size[grades_out[j3]]);
        dense.values[j3] = values_out;
        dense.grades[j3] = grades_out[j3];
        for(Py_ssize_t j0 = 0; j0 < mv0->size; j0++){
            ga_float *values0 = mv0->values[j0];
            for(Py_ssize_t j1 = 0; j1 < mv1->size; j1++){
                ga_float *values1 = mv1->values[j1];
                for(Py_ssize_t j2 = 0; j2 < mv2->size; j2++){
                    gsign = 1;
                    if(flags & OpFlag_reverse0) gsign *= (mv0->grades[j0] & 2) ? -1 : 1; // Sign for reversing the blades
                    if(flags & OpFlag_reverse1) gsign *= (mv1->grades[j1] & 2) ? -1 : 1; // Sign for reversing the blades
                    if(flags & OpFlag_reverse2) gsign *= (mv2->grades[j2] & 2) ? -1 : 1; // Sign for reversing the blades
                    
                    ga_float *values2 = mv2->values[j2];
                    smap = map[mv0->grades[j0]][mv1->grades[j1]][mv2->grades[j2]][dense.grades[j3]];
                    
                    for(Py_ssize_t j4 = 0; j4 < smap.size; j4++){
                        SparseMapBase base = smap.base[j4];
                        values_out[base.position[3]] += scalar*gsign*values0[base.position[0]]*values1[base.position[1]]*values2[base.position[2]]*base.sign;
                    }
                }
            }
        }
    }

    *mv = dense;

    return 1;
}

static int apply_multilinear_operator(void *out, void *data0, PyOperatorObject *op){
    
    if(op->flags & OpFlag_gradepreserving){ // Check if the grade preserving flag is set

    }
    return 1;
}

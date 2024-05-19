#include <Python.h>
#include "types.h"
#include "common.h"
#include "multivector_types.h"

// macro to append bitmap and values to the end of the graph
#define GRAPH_APPEND_NEXT(bitmap_,graph,addr,value_,size) \
    if(!addr[bitmap_]){        \
        graph->next = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement));\
        graph = graph->next;\
        graph->bitmap = bitmap_;\
        graph->value = (value_);\
        addr[bitmap_] = graph;\
        graph->next = NULL;\
        (size)++;\
    } else\
        addr[bitmap_]->value += (value_);\

static Py_ssize_t *init_grade_size(PyAlgebraObject *ga){
    Py_ssize_t *gsize = (Py_ssize_t*)PyMem_RawMalloc((MAX_GRADE(ga)+1)*sizeof(Py_ssize_t));
    if(!gsize){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory for grade size array");
        return NULL;
    }
    for(Py_ssize_t i = 0; i <= MAX_GRADE(ga); i++)
        gsize[i] = 0;
    return gsize;
}


static BladesMultivector sparse_dense_to_blades_sparse(SparseMultivector dense, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    Py_ssize_t ssize = 0, grade = -1;
    Py_ssize_t *gsize = init_grade_size(ga);
    Py_ssize_t *gindex = init_grade_size(ga);
    if(!gsize || !gindex){
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        return sparse;
    }
    int bitmap;
    for(Py_ssize_t i = 0; i < dense.size; i++){
        if(dense.bitmap[i] == -1) continue;
        grade = GRADE(dense.bitmap[i]);
        if(!gsize[grade]) gindex[grade] = ssize++; // first time incrementing
        gsize[grade]++;
    }

    if(!ssize){
        sparse.data = NULL;
        sparse.grade = NULL;
        sparse.size = 0;
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        return sparse;
    }

    sparse.data = (SparseMultivector*)PyMem_RawMalloc(ssize*sizeof(SparseMultivector));
    sparse.grade =  (Py_ssize_t*)PyMem_RawMalloc(ssize*sizeof(Py_ssize_t));
    if(!sparse.data || !sparse.grade){
        PyMem_RawFree(gsize);
        PyMem_RawFree(gindex);
        sparse.size = -1;
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory");
        return sparse;
    }
    sparse.size = ssize;

    // initialize each grade
    for(Py_ssize_t i = 0; i <= MAX_GRADE(ga); i++){ // iterate over grades
        if(!gsize[i]) continue;
        sparse.data[gindex[i]] = init_sparse_empty(gsize[i]);
        sparse.grade[gindex[i]] = i;
    }

    for(Py_ssize_t i = 0; i < dense.size; i++){
        bitmap = dense.bitmap[i];
        if(bitmap == -1) continue;
        grade = GRADE(bitmap); gsize[grade]--;
        sparse.data[gindex[grade]].bitmap[gsize[grade]] = bitmap;
        sparse.data[gindex[grade]].value[gsize[grade]] = dense.value[i];
    }

    PyMem_RawFree(gsize);
    PyMem_RawFree(gindex);
    return sparse;
}



static BladesMultivector blades_init_(int *bitmap, ga_float *value, Py_ssize_t size, PyAlgebraObject *ga){
    if(!size){
        BladesMultivector blades = {.size = 0,.grade = NULL, .data = NULL};
        return blades;
    }
    SparseMultivector ssparse = {.bitmap = bitmap, .value = value, .size = size};
    return sparse_dense_to_blades_sparse(ssparse,ga);
}

static int cast_to_blades(PyMultivectorIter *from, void *to, PyAlgebraObject *GA){
    BladesMultivector *pblades = (BladesMultivector*)to;
    
    if(!from || !pblades){
        return 0;
    }

    SparseMultivector sparse = {.size = from->niters, .value = NULL, .bitmap = NULL};
    sparse.value = (ga_float*)PyMem_RawMalloc(from->niters*sizeof(ga_float));
    sparse.bitmap = (int*)PyMem_RawMalloc(from->niters*sizeof(int));
    Py_ssize_t i = 0;
    while(from->next(from)){
        sparse.value[i] = from->value;
        sparse.bitmap[i] = from->bitmap;
        i++;
    }
    *pblades = sparse_dense_to_blades_sparse(sparse,GA);
    sparse_free_(sparse);
    
    return 1;
}


static SparseMultivector graph_to_sparse_multivector(BasisElement *graph, Py_ssize_t size){
    SparseMultivector sparse = alloc_sparse(size);
    BasisElement *prev;

    Py_ssize_t i = 0;
    while(graph){
        sparse.value[i] = graph->value;
        sparse.bitmap[i] = graph->bitmap;
        i++;
        prev = graph;
        graph = graph->next;
        PyMem_RawFree(prev);
    }
    return sparse;
}

static void graph_free(BasisElement *graph){
    BasisElement *prev;

    while(graph){
        prev = graph;
        graph = graph->next;
        PyMem_RawFree(prev);
    }
}

static BasisElement* graph_remove_rel_small(BasisElement *graph, Py_ssize_t *size, ga_float percentage){
    ga_float max = 0;
    BasisElement *head = graph;
    BasisElement *prev = NULL;

    while(graph){
        if(max < fabsl(graph->value))
            max = fabsl(graph->value);

        graph = graph->next;
    }
    graph = head;

    while(graph){
        if(fabsl(graph->value) < max*percentage){ // is the value relatively small
            // the previous does not change prev <- prev
            if(prev){
                prev->next = graph->next; // skip the element
                PyMem_RawFree(graph);
                graph = prev->next;
            } else{ // The head of the graph
                head = graph->next;
                PyMem_RawFree(graph);
                graph = head;
            }
            
            (*size)--;
        }else{
            prev = graph; // remember the previous
            graph = graph->next; // next element in the graph
        }
        
    }
    return head;

}


static SparseMultivector binary_sparse_geometricproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_geometricproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }
    graph->next = NULL;
    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph 
    PyMem_RawFree(addr);
    PyMem_RawFree(head);
    return sparse;
}


static SparseMultivector ternary_sparse_geometricproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head1 = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head1 of the graph1
    BasisElement *graph1 = head1;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }

    memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
    
    size = 0;
    graph = head->next;
    while(graph){
        for(Py_ssize_t i = 0; i < sparse2.size; i++){
            if(!(sign = m.sign[graph->bitmap][sparse2.bitmap[i]])) continue;
            bitmap = graph->bitmap ^ sparse2.bitmap[i];
            
            GRAPH_APPEND_NEXT(bitmap,graph1,addr,graph->value*sparse2.value[i]*sign,size)
        }

        graph = graph->next;
    }
    

    graph1->next = NULL;
    graph1 = graph_remove_rel_small(head1->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph1,size); // also frees memory for the graph
    
    graph_free(head);
    PyMem_RawFree(addr);
    PyMem_RawFree(head1);

    return sparse;
}

static SparseMultivector ternary_sparse_geometricproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_outerproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_outerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }
    graph->next = NULL;
    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph 
    PyMem_RawFree(addr);
    PyMem_RawFree(head);
    return sparse;
}


static SparseMultivector ternary_sparse_outerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head1 = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head1 of the graph1
    BasisElement *graph1 = head1;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }

    memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
    
    size = 0;
    graph = head->next;
    while(graph){
        for(Py_ssize_t i = 0; i < sparse2.size; i++){
            if(!(sign = m.sign[graph->bitmap][sparse2.bitmap[i]])) continue;
            bitmap = graph->bitmap ^ sparse2.bitmap[i];
            if(GRADE(graph->bitmap)+GRADE(sparse2.bitmap[i])!=GRADE(bitmap)) continue;
            GRAPH_APPEND_NEXT(bitmap,graph1,addr,graph->value*sparse2.value[i]*sign,size)
        }

        graph = graph->next;
    }
    

    graph1->next = NULL;
    graph1 = graph_remove_rel_small(head1->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph1,size); // also frees memory for the graph
    
    graph_free(head);
    PyMem_RawFree(addr);
    PyMem_RawFree(head1);

    return sparse;
}

static SparseMultivector ternary_sparse_outerproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(GRADE(sparse0.bitmap[i])+GRADE(sparse1.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            if(GRADE(dense1.bitmap[i])+GRADE(sparse2.bitmap[j])!=GRADE(bitmap)) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_innerproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(labs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_innerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(labs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }
    graph->next = NULL;
    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph 
    PyMem_RawFree(addr);
    PyMem_RawFree(head);
    return sparse;
}


static SparseMultivector ternary_sparse_innerproduct_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head1 = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head1 of the graph1
    BasisElement *graph1 = head1;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            if(!(sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]])) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(labs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign,size)
        }
    }

    memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
    
    size = 0;
    graph = head->next;
    while(graph){
        for(Py_ssize_t i = 0; i < sparse2.size; i++){
            if(!(sign = m.sign[graph->bitmap][sparse2.bitmap[i]])) continue;
            bitmap = graph->bitmap ^ sparse2.bitmap[i];
            if(labs((_grade0=GRADE(graph->bitmap))-(_grade1=GRADE(sparse2.bitmap[i])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            GRAPH_APPEND_NEXT(bitmap,graph1,addr,graph->value*sparse2.value[i]*sign,size)
        }

        graph = graph->next;
    }
    

    graph1->next = NULL;
    graph1 = graph_remove_rel_small(head1->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph1,size); // also frees memory for the graph
    
    graph_free(head);
    PyMem_RawFree(addr);
    PyMem_RawFree(head1);

    return sparse;
}

static SparseMultivector ternary_sparse_innerproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, SparseMultivector sparse2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    if(dense0.size == -1) return sparse;
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;

    // dense0 = product(sparse0,sparse1)
    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            sign = m.sign[sparse0.bitmap[i]][sparse1.bitmap[j]];
            if(!sign) continue;
            bitmap = sparse0.bitmap[i] ^ sparse1.bitmap[j];
            if(labs((_grade0=GRADE(sparse0.bitmap[i]))-(_grade1=GRADE(sparse1.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += sparse0.value[i]*sparse1.value[j]*sign;
        }
    }
    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    // dense1 = copy(dense0)
    // dense0 = reset(dense0)
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.value[size] = dense0.value[i];
            dense1.bitmap[size] = dense0.bitmap[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    // dense0 = product(dense1,sparse2)
    size = 0;
    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < sparse2.size; j++){
            sign = m.sign[dense1.bitmap[i]][sparse2.bitmap[j]];
            if(!sign) continue;
            bitmap = dense1.bitmap[i] ^ sparse2.bitmap[j];
            if(labs((_grade0=GRADE(dense1.bitmap[i]))-(_grade1=GRADE(sparse2.bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
            dense0.value[bitmap] += dense1.value[i]*sparse2.value[j]*sign;
        }
    }

    sparse_remove_small(dense0,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense0,size);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}


static SparseMultivector binary_sparse_regressiveproduct0_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    SparseMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    Py_ssize_t _grade0;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap,inner_bitmap;
    int undualsign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    Py_ssize_t pss = ga->asize - 1;
    Py_ssize_t l,r; int lsign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        l = pss ^ sparse0.bitmap[i];
        _grade0 = GRADE(l);
        lsign = undualsign*dm.sign[sparse0.bitmap[i]];
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            r = pss ^ sparse1.bitmap[j];
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += sparse0.value[i]*sparse1.value[j]*m.sign[l][r]*dm.sign[sparse1.bitmap[j]]*lsign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}

static SparseMultivector binary_sparse_regressiveproduct_(SparseMultivector sparse0, SparseMultivector sparse1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    Py_ssize_t _grade0;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap,inner_bitmap;
    int undualsign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    Py_ssize_t pss = ga->asize - 1;
    Py_ssize_t l,r; int lsign,sign;

    for(Py_ssize_t i = 0; i < sparse0.size; i++){
        l = pss ^ sparse0.bitmap[i];
        _grade0 = GRADE(l);
        lsign = undualsign*dm.sign[sparse0.bitmap[i]];
        for(Py_ssize_t j = 0; j < sparse1.size; j++){
            r = pss ^ sparse1.bitmap[j];
            if(!(sign = m.sign[l][r])) continue;
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0.value[i]*sparse1.value[j]*sign*dm.sign[sparse1.bitmap[j]]*lsign,size)
        }
    }

    graph->next = NULL;
    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph
    PyMem_RawFree(addr);
    PyMem_RawFree(head);
    return sparse;
}

static int unary_sparse_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    Py_ssize_t *g = get_grade_bool(grades,grade_size,MAX_GRADE(ga) + 1);
    if(!g)
        return 0;

    int size = 0;
    for(Py_ssize_t i = 0; i < sparse0->size; i++)
        if(g[GRADE(sparse0->bitmap[i])])
            size++;

    *sparse = init_sparse_empty(size--);
    if(sparse->size == -1){
        PyMem_RawFree(g);
        return 0;
    }

    // copies the values of the selected grades
    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        if(g[GRADE(sparse0->bitmap[i])]){
            sparse->value[size] = sparse0->value[i];
            sparse->bitmap[size] = sparse0->bitmap[i];
            size--;
            if(size < 0)
                break;
        }
    }

    PyMem_RawFree(g);
    return 1;
}



static int unary_sparse_reverse(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    *sparse = init_sparse_empty(sparse0->size);

    if(sparse->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        int sign = (GRADE(sparse0->bitmap[i]) & 2) ? -1 : 1;
        sparse->value[i] = sign*sparse0->value[i];
        sparse->bitmap[i] = sparse0->bitmap[i];
    }

    return 1;
}

static int unary_sparse_dual(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    *sparse = init_sparse_empty(sparse0->size);
    if(sparse->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        Py_ssize_t bitmap = sparse0->bitmap[i];
        sparse->value[i] = dm.sign[bitmap]*sparse0->value[i];
        sparse->bitmap[i] = pss ^ bitmap;
    }

    return 1;
}

static int unary_sparse_undual(void *out, void *data0, PyAlgebraObject *ga){
    SparseMultivector *sparse0 = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;
    *sparse = init_sparse_empty(sparse0->size);
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    
    if(sparse->size == -1)
        return 0;
    int sign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        Py_ssize_t bitmap = sparse0->bitmap[i];
        sparse->value[i] = sign*dm.sign[bitmap]*sparse0->value[i];
        sparse->bitmap[i] = pss ^ bitmap;
    }

    return 1;
}

static int binary_sparse_add(void *out, void *data0, void *data1, PyAlgebraObject *ga, int sign){
    SparseMultivector *sparse0 = (SparseMultivector*)data0; 
    SparseMultivector *sparse1 = (SparseMultivector*)data1;
    SparseMultivector *sparse = (SparseMultivector*)out;
    
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(ga->product->size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    Py_ssize_t bitmap;
    Py_ssize_t size = 0;



    for(Py_ssize_t i = 0; i < sparse0->size; i++){
        bitmap = sparse0->bitmap[i];
        GRAPH_APPEND_NEXT(bitmap,graph,addr,sparse0->value[i],size)
    }

    for(Py_ssize_t i = 0; i < sparse1->size; i++){
        bitmap = sparse1->bitmap[i];
        GRAPH_APPEND_NEXT(bitmap,graph,addr,sign*sparse1->value[i],size)
    }

    graph->next = NULL;
    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    *sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph 
    PyMem_RawFree(addr);
    PyMem_RawFree(head);

    return 1;
}

static int binary_blades_add(void *out, void *data0, void *data1,  PyAlgebraObject *ga, int sign){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades1 = (BladesMultivector*)data1;
    BladesMultivector *sparse = (BladesMultivector*)out;

    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1)
        return 0;

    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga->precision;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        sub = blades0->data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sub.value[j];
        }
    }

    for(Py_ssize_t i = 0; i < blades1->size; i++){
        sub = blades1->data[i];
        for(Py_ssize_t j = 0; j < sub.size; j++){
            bitmap = sub.bitmap[j];
            dense.bitmap[bitmap] = bitmap;
            dense.value[bitmap] += sign*sub.value[j];
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    *sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse->size == -1){
        sparse_free_(dense);
        return 0;
    }

    sparse_free_(dense);
    return 1;
}



static BladesMultivector binary_blades_geometricproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_geometricproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.bitmap[size] = dense0.bitmap[i];
            dense1.value[size] = dense0.value[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_outerproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(grade0+grade1!=GRADE(bitmap)) continue;
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_outerproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    Py_ssize_t grade0,grade1,grade2;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(grade0+grade1!=GRADE(bitmap)) continue;
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.bitmap[size] = dense0.bitmap[i];
            dense1.value[size] = dense0.value[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            grade2 = blades2.grade[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                if(GRADE(dense1.bitmap[i])+grade2!=GRADE(bitmap)) continue;
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_innerproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;

    Py_ssize_t bitmap;
    int sign;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(labs(grade0-grade1)!=GRADE(bitmap)||!grade0||!grade1) continue;
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}

static BladesMultivector ternary_blades_innerproduct_(BladesMultivector blades0, BladesMultivector blades1, BladesMultivector blades2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    ga_float precision = ga->precision;
    Py_ssize_t size = 0;
    Py_ssize_t bitmap;
    int sign;
    Py_ssize_t grade0,grade1,grade2;
    SparseMultivector ssparse0, ssparse1, ssparse2;

    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense0 = init_sparse_empty(m.size);
    SparseMultivector dense1;
    Py_ssize_t _grade0;

    if(dense0.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = blades1.grade[j];
            for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                for(Py_ssize_t l = 0; l < ssparse0.size; l++){
                    sign = m.sign[ssparse0.bitmap[l]][ssparse1.bitmap[k]];
                    if(!sign) continue;
                    bitmap = ssparse0.bitmap[l] ^ ssparse1.bitmap[k];
                    if(labs(grade0-grade1)!=GRADE(bitmap)||!grade0||!grade1) continue;
                    if(dense0.bitmap[bitmap] == -1) dense0.bitmap[bitmap] = bitmap, size++;
                    dense0.value[bitmap] += ssparse0.value[l]*ssparse1.value[k]*sign;
                }
            }
        }
    }

    dense1 = init_sparse_empty(size--);
    if(dense1.size == -1){
        sparse_free_(dense0);
        return sparse;
    }
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        if(dense0.bitmap[i] != -1 && size >= 0){
            dense1.bitmap[size] = dense0.bitmap[i];
            dense1.value[size] = dense0.value[i];
            size--;
        }
        dense0.bitmap[i] = -1;
        dense0.value[i] = 0;
    }

    for(Py_ssize_t i = 0; i < dense1.size; i++){
        for(Py_ssize_t j = 0; j < blades2.size; j++){
            ssparse2 = blades2.data[j];
            grade2 = blades2.grade[j];
            for(Py_ssize_t k = 0; k < ssparse2.size; k++){
                sign = m.sign[dense1.bitmap[i]][ssparse2.bitmap[k]];
                if(!sign) continue;
                bitmap = dense1.bitmap[i] ^ ssparse2.bitmap[k];
                if(labs((_grade0=GRADE(dense1.bitmap[i]))-grade2)!=GRADE(bitmap)||!_grade0||!grade2) continue;
                dense0.bitmap[bitmap] = bitmap;
                dense0.value[bitmap] += dense1.value[i]*ssparse2.value[k]*sign;
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense0.size; i++)
        if(dense0.bitmap[i] != -1 && comp_abs(dense0.value[i],precision))
            dense0.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense0,ga);
    if(sparse.size == -1){
        sparse_free_(dense0);
        sparse_free_(dense1);
        return sparse;
    }

    sparse_free_(dense0);
    sparse_free_(dense1);
    return sparse;
}

static BladesMultivector binary_blades_regressiveproduct_(BladesMultivector blades0, BladesMultivector blades1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize - 1;
    ga_float precision = ga->precision;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    Py_ssize_t l,r; int lsign;
    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    Py_ssize_t bitmap,inner_bitmap;
    SparseMultivector ssparse0, ssparse1;
    Py_ssize_t grade0,grade1;
    BladesMultivector sparse = {.size = -1};
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1)
        return sparse;

    for(Py_ssize_t i = 0; i < blades0.size; i++){
        ssparse0 = blades0.data[i];
        grade0 = max_grade - blades0.grade[i];
        for(Py_ssize_t j = 0; j < blades1.size; j++){
            ssparse1 = blades1.data[j];
            grade1 = max_grade - blades1.grade[j];
            for(Py_ssize_t n = 0; n < ssparse0.size; n++){
                l = pss ^ ssparse0.bitmap[n]; // product with the pseudoscalar
                lsign = undualsign*dm.sign[ssparse0.bitmap[n]];
                for(Py_ssize_t k = 0; k < ssparse1.size; k++){
                    r = pss ^ ssparse1.bitmap[k];
                    bitmap = pss^(inner_bitmap = l ^ r);
                    if(grade0 + grade1 != GRADE(inner_bitmap)) continue;
                    dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += ssparse0.value[n]*ssparse1.value[k]*m.sign[l][r]*dm.sign[ssparse1.bitmap[k]]*lsign;
                }
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    sparse = sparse_dense_to_blades_sparse(dense,ga);
    if(sparse.size == -1){
        sparse_free_(dense);
        return sparse;
    }

    sparse_free_(dense);
    return sparse;
}


static int unary_blades_dual(void *out, void *data0, PyAlgebraObject *ga){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;

    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    *blades = init_blades_empty(blades0->size);
    if(blades->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < blades0->size; i++){
        Py_ssize_t bsize = blades0->data[i].size;
        blades->data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades->data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades->data[i].bitmap || !blades->data[i].value){
            blades_free_(*blades);
            return 0;
        }
        blades->data[i].size = bsize;
        blades->grade[i] = max_grade - blades0->grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0->data[i].bitmap[j];
            blades->data[i].bitmap[j] = pss^bitmap;
            blades->data[i].value[j] = dm.sign[bitmap]*blades0->data[i].value[j];
        }
    }

    return 1;
}

static int unary_blades_undual(void *out, void *data0, PyAlgebraObject *ga){
    BladesMultivector *blades0 = (BladesMultivector*)data0;
    BladesMultivector *blades = (BladesMultivector*)out;
    DualMap dm = ga->dm;
    *blades = init_blades_empty(blades0->size);
    if(blades->size == -1)
        return 0;
    Py_ssize_t pss = ga->asize-1;
    Py_ssize_t max_grade = MAX_GRADE(ga);
    int sign = max_grade & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < blades0->size; i++){
        Py_ssize_t bsize = blades0->data[i].size;
        blades->data[i].bitmap = (int*)PyMem_RawMalloc(bsize*sizeof(int));
        blades->data[i].value = (ga_float*)PyMem_RawMalloc(bsize*sizeof(ga_float));
        if(!blades->data[i].bitmap || !blades->data[i].value){
            blades_free_(*blades);
            return 0;
        }
        blades->data[i].size = bsize;
        blades->grade[i] = max_grade - blades0->grade[i];
        for(Py_ssize_t j = 0; j < bsize; j++){
            Py_ssize_t bitmap = blades0->data[i].bitmap[j];
            blades->data[i].bitmap[j] = pss^bitmap;
            blades->data[i].value[j] = sign*dm.sign[bitmap]*blades0->data[i].value[j];
        }
    }

    return 1;
}



static DenseMultivector binary_dense_geometricproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_geometricproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_outerproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_outerproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    
    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(GRADE(i)+GRADE(j)!=GRADE(bitmap)) continue;
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_innerproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap;
    Py_ssize_t _grade0, _grade1;

    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    int sign;
    for(Py_ssize_t i = 0; i < dense0.size; i++){
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(labs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    return dense;
}

static DenseMultivector ternary_dense_innerproduct_(DenseMultivector dense0, DenseMultivector dense1, DenseMultivector dense2, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = {.size = -1};
    DenseMultivector temp = {.size = -1};
    dense = init_dense_empty(m.size);
    temp = init_dense_empty(m.size);
    if(temp.size == -1) return temp;
    if(dense.size == -1) return dense;
    int sign;
    Py_ssize_t bitmap;
    Py_ssize_t _grade0, _grade1;

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(labs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            temp.value[bitmap] += dense0.value[i]*dense1.value[j]*sign;
        }
    }

    for(Py_ssize_t i = 0; i < m.size; i++){
        for(Py_ssize_t j = 0; j < m.size; j++){
            sign = m.sign[i][j];
            if(!sign) continue;
            bitmap = i^j;
            if(labs((_grade0=GRADE(i))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            dense.value[bitmap] += temp.value[i]*dense2.value[j]*sign;
        }
    }
    dense_free_(temp);
    return dense;
}

static DenseMultivector binary_dense_regressiveproduct_(DenseMultivector dense0, DenseMultivector dense1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    DenseMultivector dense = {.size = -1};
    Py_ssize_t bitmap,inner_bitmap;
    Py_ssize_t _grade0;
    Py_ssize_t pss = ga->asize - 1;
    dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;

    Py_ssize_t l,r; int lsign;
    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar

    for(Py_ssize_t i = 0; i < dense0.size; i++){
        l = pss ^ i; // dual of the left multivector
        lsign = undualsign*dm.sign[i];
        _grade0 = GRADE(l);
        for(Py_ssize_t j = 0; j < dense1.size; j++){
            r = pss ^ j; // dual of the right multivector
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            dense.value[bitmap] += dense0.value[i]*dense1.value[j]*m.sign[l][r]*dm.sign[j]*lsign;
        }
    }

    return dense;
}

static int unary_dense_gradeproject(void *out, void *data0, PyAlgebraObject *ga, int *grades, Py_ssize_t grade_size){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;
    
    Py_ssize_t *g = get_grade_bool(grades,grade_size,MAX_GRADE(ga)+1);
    if(!g) return 0;
    *dense = init_dense_empty(dense0->size);
    for(Py_ssize_t i = 0; i < dense->size; i++)
        if(g[GRADE(i)])
            dense->value[i] = dense0->value[i];

    PyMem_RawFree(g);
    return 1;
}

static int unary_dense_reverse(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;

    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < dense0->size; i++){
        int sign = (GRADE(i) & 2) ? -1 : 1;
        dense->value[i] = sign*dense0->value[i];
    }

    return 1;
}

static int unary_dense_dual(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;

    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[pss^i] = dm.sign[i]*dense0->value[i];

    return 1;
}

static int unary_dense_undual(void *out, void *data0, PyAlgebraObject *ga){
    DenseMultivector *dense0 = (DenseMultivector*)data0;
    DenseMultivector *dense = (DenseMultivector*)out;

    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize-1;
    *dense = init_dense_empty(dense0->size);
    if(dense->size == -1)
        return 0;

    int sign = MAX_GRADE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    for(Py_ssize_t i = 0; i < dense0->size; i++)
        dense->value[pss^i] = sign*dm.sign[i]*dense0->value[i];

    return 1;
}




static SparseMultivector atomic_sparse_geometricproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head_out = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph_out = head_out;

    BasisElement *prev = NULL;

    graph->next = NULL;
    graph->bitmap = 0;
    graph->value = 1;

    graph_out->next = NULL;
    graph_out->bitmap = -1;
    graph_out->value = 0;
    

    int sign; Py_ssize_t bitmap,size=0;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        size = 0;
        while(graph){
            for(Py_ssize_t j = 0; j < data[i].size; j++){
                sign = m.sign[graph->bitmap][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = graph->bitmap ^ data[i].bitmap[j];
                
                GRAPH_APPEND_NEXT(bitmap,graph_out,addr,data[i].value[j]*graph->value*sign,size)
            }
            graph =  graph->next;
        }
        
        graph_free(prev); // release memory for graph needed in the previous computation 
        memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
        
        prev = head;
        graph = head = head_out->next;
        head_out->next = NULL;
        graph_out = head_out;
    }

    graph_out->next = NULL;
    graph_out = graph_remove_rel_small(head,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph_out,size); // also frees memory for the graph
    
    graph_free(prev);
    PyMem_RawFree(addr);
    PyMem_RawFree(head_out);

    return sparse;
}


static SparseMultivector atomic_sparse_geometricproduct0_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static SparseMultivector atomic_sparse_outerproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head_out = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph_out = head_out;

    BasisElement *prev = NULL;

    graph->next = NULL;
    graph->bitmap = 0;
    graph->value = 1;

    graph_out->next = NULL;
    graph_out->bitmap = -1;
    graph_out->value = 0;
    

    int sign; Py_ssize_t bitmap,size=0;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        size = 0;
        while(graph){
            for(Py_ssize_t j = 0; j < data[i].size; j++){
                sign = m.sign[graph->bitmap][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = graph->bitmap ^ data[i].bitmap[j];
                if(GRADE(graph->bitmap)+GRADE(data[i].bitmap[j])!=GRADE(bitmap)) continue;
                GRAPH_APPEND_NEXT(bitmap,graph_out,addr,data[i].value[j]*graph->value*sign,size)
            }
            graph =  graph->next;
        }
        
        graph_free(prev); // release memory for graph needed in the previous computation 
        memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
        
        prev = head;
        graph = head = head_out->next;
        head_out->next = NULL;
        graph_out = head_out;
    }

    graph_out->next = NULL;
    graph_out = graph_remove_rel_small(head,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph_out,size); // also frees memory for the graph
    
    graph_free(prev);
    PyMem_RawFree(addr);
    PyMem_RawFree(head_out);

    return sparse;
}


static SparseMultivector atomic_sparse_outerproduct0_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                if(GRADE(temp.bitmap[k])+GRADE(data[i].bitmap[j])!=GRADE(bitmap)) continue;
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static SparseMultivector atomic_sparse_innerproduct_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector sparse = {.size = -1};

    // initiallize the addresses to null: addr[i] <- NULL
    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(m.size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    BasisElement *head_out = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph_out = head_out;

    BasisElement *prev = NULL;

    graph->next = NULL;
    graph->bitmap = 0;
    graph->value = 1;

    graph_out->next = NULL;
    graph_out->bitmap = -1;
    graph_out->value = 0;
    Py_ssize_t _grade0, _grade1;


    int sign; Py_ssize_t bitmap,size=0;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        size = 0;
        while(graph){
            for(Py_ssize_t j = 0; j < data[i].size; j++){
                sign = m.sign[graph->bitmap][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = graph->bitmap ^ data[i].bitmap[j];
                if(labs((_grade0=GRADE(graph->bitmap))-(_grade1=GRADE(data[i].bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
                GRAPH_APPEND_NEXT(bitmap,graph_out,addr,data[i].value[j]*graph->value*sign,size)
            }
            graph =  graph->next;
        }
        
        graph_free(prev); // release memory for graph needed in the previous computation 
        memset(addr,0,m.size*sizeof(BasisElement*));// reset the address array
        
        prev = head;
        graph = head = head_out->next;
        head_out->next = NULL;
        graph_out = head_out;
    }

    graph_out->next = NULL;
    graph_out = graph_remove_rel_small(head,&size,ga->precision);
    sparse = graph_to_sparse_multivector(graph_out,size); // also frees memory for the graph
    
    graph_free(prev);
    PyMem_RawFree(addr);
    PyMem_RawFree(head_out);

    return sparse;
}


static SparseMultivector atomic_sparse_innerproduct0_(SparseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    // Allocate memory for a dense y
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1) {
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    Py_ssize_t _grade0, _grade1;

    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar

    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][data[i].bitmap[j]];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ data[i].bitmap[j];
                if(labs((_grade0=GRADE(temp.bitmap[k]))-(_grade1=GRADE(data[i].bitmap[j])))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    if(sparse.size == -1){
        sparse_free_(dense);
        sparse_free_(temp);
        return sparse;
    }

    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}


static int atomic_blades_add(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t size){
    BladesMultivector *data = (BladesMultivector*)data0;
    BladesMultivector *sparse = (BladesMultivector*)out;

    SparseMultivector dense = init_sparse_empty(ga->asize);
    if(dense.size == -1) return 0;
    SparseMultivector sub;
    Py_ssize_t bitmap;
    ga_float precision = ga->precision;

    for(Py_ssize_t k = 0; k < size; k++){
        for(Py_ssize_t i = 0; i < data[k].size; i++){
            sub = data[k].data[i];
            for(Py_ssize_t j = 0; j < sub.size; j++){
                bitmap = sub.bitmap[j];
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += sub.value[j];
            }
        }
    }

    // remove small values
    for(Py_ssize_t i = 0; i < dense.size; i++)
        if(dense.bitmap[i] != -1 && comp_abs(dense.value[i],precision))
            dense.bitmap[i] = -1;

    *sparse = sparse_dense_to_blades_sparse(dense,ga);
    sparse_free_(dense);
    return 1;
}

static int atomic_sparse_add(void *out, void *data0, PyAlgebraObject *ga, Py_ssize_t dsize){
    SparseMultivector *data = (SparseMultivector*)data0;
    SparseMultivector *sparse = (SparseMultivector*)out;

    BasisElement **addr = (BasisElement**)PyMem_RawCalloc(ga->product->size,sizeof(BasisElement*)); // A list of addresses for all basis elements
    BasisElement *head = (BasisElement*)PyMem_RawMalloc(sizeof(BasisElement)); // The head of the graph
    BasisElement *graph = head;

    graph->next = NULL;
    graph->bitmap = -1;
    graph->value = 0;

    Py_ssize_t bitmap;
    Py_ssize_t size = 0;

    for(Py_ssize_t j = 0; j < dsize; j++){
        for(Py_ssize_t i = 0; i < data[j].size; i++){
            bitmap = data[j].bitmap[i];
            GRAPH_APPEND_NEXT(bitmap,graph,addr,data[j].value[i],size)
        }
    }

    graph = graph_remove_rel_small(head->next,&size,ga->precision);
    *sparse = graph_to_sparse_multivector(graph,size); // also frees memory for the graph 
    PyMem_RawFree(addr);
    PyMem_RawFree(head);
    return 1;
}

static BladesMultivector atomic_blades_geometricproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        
                                        dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += temp.value[j]*sdata.value[l]*sign;
                }
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static BladesMultivector atomic_blades_outerproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    Py_ssize_t sgrade;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                sgrade = data[i].grade[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        if(GRADE(temp.bitmap[j])+sgrade!=GRADE(bitmap)) continue;
                                        dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += temp.value[j]*sdata.value[l]*sign;
                }
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static BladesMultivector atomic_blades_innerproduct_(BladesMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    BladesMultivector sparse = {.size = -1};
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return sparse;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return sparse;
    }
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t tsize = 1;
    Py_ssize_t sgrade;
    Py_ssize_t _grade0;
    int sign; int bitmap;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < tsize; j++){ // iterate over temp
            if(temp.bitmap[j] == -1) continue; // ignore if value not set
            for(Py_ssize_t k = 0; k < data[i].size; k++){ // iterate over grades
                SparseMultivector sdata = data[i].data[k];
                sgrade = data[i].grade[k];
                for(Py_ssize_t l = 0; l < sdata.size; l++){ // iterate over values and bitmaps of data[i]
                    sign = m.sign[temp.bitmap[j]][sdata.bitmap[l]];
                    if(!sign) continue;
                    bitmap = temp.bitmap[j] ^ sdata.bitmap[l];
                                        if(labs((_grade0=GRADE(temp.bitmap[j]))-sgrade)!=GRADE(bitmap)||!_grade0||!sgrade) continue;
                                        dense.bitmap[bitmap] = bitmap;
                    dense.value[bitmap] += temp.value[j]*sdata.value[l]*sign;
                }
            }
        }
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_blades_sparse(temp,ga);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}

static DenseMultivector atomic_dense_geometricproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    
    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}
static DenseMultivector atomic_dense_outerproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    
    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                if(GRADE(k)+GRADE(j)!=GRADE(bitmap)) continue;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}
static DenseMultivector atomic_dense_innerproduct_(DenseMultivector *data, Py_ssize_t dsize, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DenseMultivector dense = init_dense_empty(m.size);
    if(dense.size == -1) return dense;
    DenseMultivector temp = init_dense_empty(m.size);
    if(temp.size == -1) {
        dense_free_(dense);
        return temp;
    }

    *temp.value = 1; // initialize temp to unit scalar
    Py_ssize_t _grade0, _grade1;

    Py_ssize_t bitmap;
    int sign;
    for(Py_ssize_t i = 0; i < dsize; i++){ // iterate over multivectors
        for(Py_ssize_t j = 0; j < data[i].size; j++){
            for(Py_ssize_t k = 0; k < temp.size; k++){
                sign = m.sign[k][j];
                if(!sign) continue;
                bitmap = k ^ j;
                if(labs((_grade0=GRADE(k))-(_grade1=GRADE(j)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
                dense.value[bitmap] += temp.value[k]*data[i].value[j]*sign;
            }
        }
        // copy values
        for(Py_ssize_t l = 0; l < dense.size; l++){
            temp.value[l] = dense.value[l];
            dense.value[l] = 0;
        }
    }

    dense_free_(dense);
    return temp;
}

static SparseMultivector binary_mixed_regressiveproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    DualMap dm = ga->dm;
    Py_ssize_t pss = ga->asize - 1;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    int undualsign = METRIC_SIZE(ga) & 2 ? -1 : 1; // sign of reversing the pseudoscalar
    Py_ssize_t _grade0;
    SparseMultivector sparse;
    Py_ssize_t size = 0;
    Py_ssize_t l,r; int lsign;

    Py_ssize_t bitmap, inner_bitmap;
    while(iter0->next(iter0)){
        l = pss^iter0->bitmap;
        lsign = undualsign*dm.sign[iter0->bitmap];
        _grade0 = GRADE(l);
        while(iter1->next(iter1)){
            r = pss^iter1->bitmap;
            bitmap = pss^(inner_bitmap = l^r);
            if(_grade0 + GRADE(r) != GRADE(inner_bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*m.sign[l][r]*lsign*dm.sign[iter1->bitmap];
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}



static SparseMultivector binary_mixed_geometricproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    
    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_geometricproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*iter->value*sign;
            }
        }iter++;
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static SparseMultivector binary_mixed_outerproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    
    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            if(GRADE(iter0->bitmap)+GRADE(iter1->bitmap)!=GRADE(bitmap)) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_outerproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    
    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                if(GRADE(temp.bitmap[k])+GRADE(iter->bitmap)!=GRADE(bitmap)) continue;
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*iter->value*sign;
            }
        }iter++;
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
static SparseMultivector binary_mixed_innerproduct_(PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;

    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;

    SparseMultivector sparse;
    Py_ssize_t size = 0;
    Py_ssize_t _grade0, _grade1;

    int sign; Py_ssize_t bitmap;
    while(iter0->next(iter0)){
        while(iter1->next(iter1)){
            sign = m.sign[iter0->bitmap][iter1->bitmap];
            if(!sign) continue;
            bitmap = iter0->bitmap ^ iter1->bitmap;
            if(labs((_grade0=GRADE(iter0->bitmap))-(_grade1=GRADE(iter1->bitmap)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
            if(dense.bitmap[bitmap] == -1) dense.bitmap[bitmap] = bitmap, size++;
            dense.value[bitmap] += iter0->value*iter1->value*sign;
        }
    }

    sparse_remove_small(dense,ga->precision,&size);
    sparse = sparse_dense_to_sparse_sparse(dense,size);
    sparse_free_(dense);
    return sparse;
}


static SparseMultivector atomic_mixed_innerproduct_(PyMultivectorIter *iter, Py_ssize_t size, PyAlgebraObject *ga){
    CliffordMap m = *ga->product;
    SparseMultivector dense = init_sparse_empty(m.size);
    if(dense.size == -1) return dense;
    SparseMultivector temp = init_sparse_empty(m.size);
    if(temp.size == -1){
        sparse_free_(dense);
        return temp;
    }

    SparseMultivector sparse;
    Py_ssize_t tsize = 1;
    Py_ssize_t _grade0, _grade1;

    *temp.bitmap = 0; *temp.value = 1; // initialize temp to unit scalar
    int sign; Py_ssize_t bitmap;
    for(Py_ssize_t i = 0; i < size; i++){ // iterate over multivectors
        while(iter->next(iter)){
            for(Py_ssize_t k = 0; k < tsize; k++){
                if(temp.bitmap[k] == -1) continue;
                sign = m.sign[temp.bitmap[k]][iter->bitmap];
                if(!sign) continue;
                bitmap = temp.bitmap[k] ^ iter->bitmap;
                if(labs((_grade0=GRADE(temp.bitmap[k]))-(_grade1=GRADE(iter->bitmap)))!=GRADE(bitmap)||!_grade0||!_grade1) continue;
                dense.bitmap[bitmap] = bitmap;
                dense.value[bitmap] += temp.value[k]*iter->value*sign;
            }
        }iter++;
        tsize = 0;
        for(Py_ssize_t l = 0; l < dense.size; l++){
            if(dense.bitmap[l] != -1){
                temp.bitmap[tsize] = dense.bitmap[l];
                temp.value[tsize] = dense.value[l];
                tsize++;
            }
            dense.bitmap[l] = -1;
            dense.value[l] = 0;
        }
    }

    sparse_remove_small(temp,ga->precision,&tsize);
    sparse = sparse_dense_to_sparse_sparse(temp,tsize);
    sparse_free_(dense);
    sparse_free_(temp);
    return sparse;
}
/*

   
static int unary_sparse_gradeproject__(void *out, void *data0, int *grades, Py_ssize_t size){
    SparseMultivector *sparse_out = (SparseMultivector*)out;
    SparseMultivector *sparse0 = (SparseMultivector*)data0;

    *sparse_out = unary_sparse_gradeproject_(*sparse0,data0->GA, grades, size);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        return 0;
    }

    return 1;
}

  
static int unary_sparse_reverse__(void *out, void *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)out;
    SparseMultivector *sparse0 = (SparseMultivector*)data0;

    *sparse_out = unary_sparse_reverse_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        return 0;
    }

    return 1;
}

  
static int unary_sparse_dual__(void *out, void *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)out;
    SparseMultivector *sparse0 = (SparseMultivector*)data0;

    *sparse_out = unary_sparse_dual_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        return 0;
    }

    return 1;
}

  
static int unary_sparse_undual__(void *out, void *data0){
    SparseMultivector *sparse_out = (SparseMultivector*)out;
    SparseMultivector *sparse0 = (SparseMultivector*)data0;

    *sparse_out = unary_sparse_undual_(*sparse0,data0->GA);

    if(sparse_out->size == -1){
        PyMem_RawFree(sparse_out);
        return 0;
    }

    return 1;
}

   
static int unary_dense_gradeproject__(void *out, void *data0, int *grades, Py_ssize_t size){
    DenseMultivector *dense_out = (DenseMultivector*)out;
    DenseMultivector *dense0 = (DenseMultivector*)data0;

    *dense_out = unary_dense_gradeproject_(*dense0,data0->GA, grades, size);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        return 0;
    }

    return 1;
}

  
static int unary_dense_reverse__(void *out, void *data0){
    DenseMultivector *dense_out = (DenseMultivector*)out;
    DenseMultivector *dense0 = (DenseMultivector*)data0;

    *dense_out = unary_dense_reverse_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        return 0;
    }

    return 1;
}

  
static int unary_dense_dual__(void *out, void *data0){
    DenseMultivector *dense_out = (DenseMultivector*)out;
    DenseMultivector *dense0 = (DenseMultivector*)data0;

    *dense_out = unary_dense_dual_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        return 0;
    }

    return 1;
}

  
static int unary_dense_undual__(void *out, void *data0){
    DenseMultivector *dense_out = (DenseMultivector*)out;
    DenseMultivector *dense0 = (DenseMultivector*)data0;

    *dense_out = unary_dense_undual_(*dense0,data0->GA);

    if(dense_out->size == -1){
        PyMem_RawFree(dense_out);
        return 0;
    }

    return 1;
}

     
static int unary_blades_dual__(void *out, void *data0){
    BladesMultivector *blades_out = (BladesMultivector*)out;
    BladesMultivector *blades0 = (BladesMultivector*)data0;

    *blades_out = unary_blades_dual_(*blades0,data0->GA);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        return 0;
    }

    return 1;
}

  
static int unary_blades_undual__(void *out, void *data0){
    BladesMultivector *blades_out = (BladesMultivector*)out;
    BladesMultivector *blades0 = (BladesMultivector*)data0;

    *blades_out = unary_blades_undual_(*blades0,data0->GA);

    if(blades_out->size == -1){
        PyMem_RawFree(blades_out);
        return 0;
    }

    return 1;
}


*/

  static int binary_sparse_product(void *out,void *data0, void *data1,  PyAlgebraObject *GA, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)data0;
    SparseMultivector *psparse1 = (SparseMultivector*)data1;
    SparseMultivector *psparse  = (SparseMultivector*)out;
    if(!psparse0 ||!psparse1 ||!psparse || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = binary_sparse_geometricproduct_(*psparse0,*psparse1,GA);
            break;
        case ProductType_inner:
            *psparse = binary_sparse_innerproduct_(*psparse0,*psparse1,GA);
            break;
        case ProductType_outer:
            *psparse = binary_sparse_outerproduct_(*psparse0,*psparse1,GA);
            break;
        case ProductType_regressive:
            *psparse = binary_sparse_regressiveproduct_(*psparse0,*psparse1,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
 static int ternary_sparse_product(void *out,void *data0, void *data1, void *data2,  PyAlgebraObject *GA, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)data0;
    SparseMultivector *psparse1 = (SparseMultivector*)data1;
    SparseMultivector *psparse2 = (SparseMultivector*)data2;
    SparseMultivector *psparse  = (SparseMultivector*)out;
    if(!psparse0 ||!psparse1 ||!psparse2 ||!psparse || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = ternary_sparse_geometricproduct_(*psparse0,*psparse1,*psparse2,GA);
            break;
        case ProductType_inner:
            *psparse = ternary_sparse_innerproduct_(*psparse0,*psparse1,*psparse2,GA);
            break;
        case ProductType_outer:
            *psparse = ternary_sparse_outerproduct_(*psparse0,*psparse1,*psparse2,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
  static int binary_dense_product(void *out,void *data0, void *data1,  PyAlgebraObject *GA, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)data0;
    DenseMultivector *pdense1 = (DenseMultivector*)data1;
    DenseMultivector *pdense  = (DenseMultivector*)out;
    if(!pdense0 ||!pdense1 ||!pdense || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = binary_dense_geometricproduct_(*pdense0,*pdense1,GA);
            break;
        case ProductType_inner:
            *pdense = binary_dense_innerproduct_(*pdense0,*pdense1,GA);
            break;
        case ProductType_outer:
            *pdense = binary_dense_outerproduct_(*pdense0,*pdense1,GA);
            break;
        case ProductType_regressive:
            *pdense = binary_dense_regressiveproduct_(*pdense0,*pdense1,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
 static int ternary_dense_product(void *out,void *data0, void *data1, void *data2,  PyAlgebraObject *GA, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)data0;
    DenseMultivector *pdense1 = (DenseMultivector*)data1;
    DenseMultivector *pdense2 = (DenseMultivector*)data2;
    DenseMultivector *pdense  = (DenseMultivector*)out;
    if(!pdense0 ||!pdense1 ||!pdense2 ||!pdense || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = ternary_dense_geometricproduct_(*pdense0,*pdense1,*pdense2,GA);
            break;
        case ProductType_inner:
            *pdense = ternary_dense_innerproduct_(*pdense0,*pdense1,*pdense2,GA);
            break;
        case ProductType_outer:
            *pdense = ternary_dense_outerproduct_(*pdense0,*pdense1,*pdense2,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
  static int binary_blades_product(void *out,void *data0, void *data1,  PyAlgebraObject *GA, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)data0;
    BladesMultivector *pblades1 = (BladesMultivector*)data1;
    BladesMultivector *pblades  = (BladesMultivector*)out;
    if(!pblades0 ||!pblades1 ||!pblades || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = binary_blades_geometricproduct_(*pblades0,*pblades1,GA);
            break;
        case ProductType_inner:
            *pblades = binary_blades_innerproduct_(*pblades0,*pblades1,GA);
            break;
        case ProductType_outer:
            *pblades = binary_blades_outerproduct_(*pblades0,*pblades1,GA);
            break;
        case ProductType_regressive:
            *pblades = binary_blades_regressiveproduct_(*pblades0,*pblades1,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
 static int ternary_blades_product(void *out,void *data0, void *data1, void *data2,  PyAlgebraObject *GA, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)data0;
    BladesMultivector *pblades1 = (BladesMultivector*)data1;
    BladesMultivector *pblades2 = (BladesMultivector*)data2;
    BladesMultivector *pblades  = (BladesMultivector*)out;
    if(!pblades0 ||!pblades1 ||!pblades2 ||!pblades || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = ternary_blades_geometricproduct_(*pblades0,*pblades1,*pblades2,GA);
            break;
        case ProductType_inner:
            *pblades = ternary_blades_innerproduct_(*pblades0,*pblades1,*pblades2,GA);
            break;
        case ProductType_outer:
            *pblades = ternary_blades_outerproduct_(*pblades0,*pblades1,*pblades2,GA);
            break;
        default:
            return 0;
    }

    return 1;
}


static int atomic_sparse_product(void *out, void *data0, PyAlgebraObject *GA, Py_ssize_t size, ProductType ptype){
    SparseMultivector *psparse0 = (SparseMultivector*)data0;
    SparseMultivector *psparse  = (SparseMultivector*)out;
    if(!psparse0 || !psparse){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = atomic_sparse_geometricproduct_(psparse0,size,GA);
            break;
        case ProductType_inner:
            *psparse = atomic_sparse_innerproduct_(psparse0,size,GA);
            break;
        case ProductType_outer:
            *psparse = atomic_sparse_outerproduct_(psparse0,size,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
static int atomic_dense_product(void *out, void *data0, PyAlgebraObject *GA, Py_ssize_t size, ProductType ptype){
    DenseMultivector *pdense0 = (DenseMultivector*)data0;
    DenseMultivector *pdense  = (DenseMultivector*)out;
    if(!pdense0 || !pdense){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pdense = atomic_dense_geometricproduct_(pdense0,size,GA);
            break;
        case ProductType_inner:
            *pdense = atomic_dense_innerproduct_(pdense0,size,GA);
            break;
        case ProductType_outer:
            *pdense = atomic_dense_outerproduct_(pdense0,size,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
static int atomic_blades_product(void *out, void *data0, PyAlgebraObject *GA, Py_ssize_t size, ProductType ptype){
    BladesMultivector *pblades0 = (BladesMultivector*)data0;
    BladesMultivector *pblades  = (BladesMultivector*)out;
    if(!pblades0 || !pblades){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *pblades = atomic_blades_geometricproduct_(pblades0,size,GA);
            break;
        case ProductType_inner:
            *pblades = atomic_blades_innerproduct_(pblades0,size,GA);
            break;
        case ProductType_outer:
            *pblades = atomic_blades_outerproduct_(pblades0,size,GA);
            break;
        default:
            return 0;
    }

    return 1;
}
/*
  static int atomic_blades_add(void *out, void *data, PyAlgebraObject *GA, Py_ssize_t size){
    BladesMultivector *blades = (BladesMultivector*)out;
    BladesMultivector *blades_array = (BladesMultivector*)data;

    *blades = atomic_blades_add_(blades_array,size,GA);

    if(blades->size == -1){
        return 0;
    }

    return out;
}
*/

static int atomic_mixed_product(void *out, PyMultivectorIter *iter, PyAlgebraObject *GA, Py_ssize_t size, ProductType ptype){    
    SparseMultivector *psparse = (SparseMultivector*)out;

    if(!iter || !psparse || !out){
        return 0; // raise error
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = atomic_mixed_geometricproduct_(iter,size,GA);
            break;
        case ProductType_inner:
            *psparse = atomic_mixed_innerproduct_(iter,size,GA);
            break;
        case ProductType_outer:
            *psparse = atomic_mixed_outerproduct_(iter,size,GA);
            break;
        default:
            return 0;
    }

    return 1;
}

static int binary_mixed_product(void *out, PyMultivectorIter *iter0, PyMultivectorIter *iter1, PyAlgebraObject *GA, ProductType ptype){
    SparseMultivector *psparse = (SparseMultivector*)out;
    if(!iter1 || !iter0 || !psparse || !out ){
        return 0;
    }

    switch(ptype){
        case ProductType_geometric:
            *psparse = binary_mixed_geometricproduct_(iter0,iter1,GA);
            break;
        case ProductType_inner:
            *psparse = binary_mixed_innerproduct_(iter0,iter1,GA);
            break;
        case ProductType_outer:
            *psparse = binary_mixed_outerproduct_(iter0,iter1,GA);
            break;
        case ProductType_regressive:
            *psparse = binary_mixed_regressiveproduct_(iter0,iter1,GA);
            break;
        default:
            return 0;
    }

    if(psparse->size == -1){
        return 0;
    }
    
    return 1;
}


static int blades_init(void *data, PyAlgebraObject *ga, int *bitmap, ga_float *value, Py_ssize_t size){
    BladesMultivector *blades = data;
    *blades = blades_init_(bitmap,value,size,ga);
    return 1;
}


PyMultivectorMixedMath_Funcs largemultivector_mixed_fn = {
  .add = NULL,
  .product = binary_mixed_product,
  .atomic_add = NULL,
  .atomic_product = atomic_mixed_product,
  .type_names = {"sparselarge","denselarge","bladeslarge",NULL},
};


PyMultivectorMath_Funcs largemultivector_dense_math_fn = {
    .product = binary_dense_product,
    .atomic_product = atomic_dense_product,
    .ternary_product = ternary_dense_product,
    .grade_project = unary_dense_gradeproject,
    .reverse = unary_dense_reverse,
    .add = NULL,
    .atomic_add = NULL,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_dense_dual,
    .undual = unary_dense_undual,
};

PyMultivectorMath_Funcs largemultivector_sparse_math_fn = {
    .product = binary_sparse_product,
    .atomic_product = atomic_sparse_product,
    .ternary_product = ternary_sparse_product,
    .grade_project = unary_sparse_gradeproject,
    .reverse = unary_sparse_reverse,
    .add = binary_sparse_add,
    .atomic_add = atomic_sparse_add,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_sparse_dual,
    .undual = unary_sparse_undual,
};

PyMultivectorMath_Funcs largemultivector_blades_math_fn = {
    .product = binary_blades_product,
    .atomic_product = atomic_blades_product,
    .ternary_product = ternary_blades_product,
    .grade_project = NULL,
    .reverse = NULL,
    .add = binary_blades_add,
    .atomic_add = atomic_blades_add,
    .scalar_product = NULL,
    .scalar_add = NULL,
    .dual = unary_blades_dual,
    .undual = unary_blades_undual,
};


 
PyMultivectorData_Funcs largemultivector_sparse_data_fn = {
    .free = NULL,
    .init = NULL,
    .iter_next = NULL,
    .iter_init = NULL,
};
 
PyMultivectorData_Funcs largemultivector_dense_data_fn = {
    .free = NULL,
    .init = NULL,
    .iter_next = NULL,
    .iter_init = NULL,
};

PyMultivectorData_Funcs largemultivector_blades_data_fn = {
    .free = NULL,
    .init = blades_init,
    .cast = cast_to_blades,
    .iter_next = NULL,
    .iter_init = NULL,
};


 const PyMultivectorSubType largesparse_subtype = {
    .math_funcs = &largemultivector_sparse_math_fn,
    .data_funcs = &largemultivector_sparse_data_fn,
    .name = "",
    .type_name = "sparselarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_sparse,
    .basic_size = sizeof(SparseMultivector),
};
 const PyMultivectorSubType largedense_subtype = {
    .math_funcs = &largemultivector_dense_math_fn,
    .data_funcs = &largemultivector_dense_data_fn,
    .name = "",
    .type_name = "denselarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_dense,
    .basic_size = sizeof(DenseMultivector),
};
 const PyMultivectorSubType largeblades_subtype = {
    .math_funcs = &largemultivector_blades_math_fn,
    .data_funcs = &largemultivector_blades_data_fn,
    .name = "",
    .type_name = "bladeslarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_blades,
    .basic_size = sizeof(BladesMultivector),
};


PyMultivectorSubType largemultivector_subtypes_array[3] = {
  {
    .math_funcs = &largemultivector_sparse_math_fn,
    .data_funcs = &largemultivector_sparse_data_fn,
    .name = "",
    .type_name = "sparselarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_sparse,
    .basic_size = sizeof(SparseMultivector),
},
  {
    .math_funcs = &largemultivector_dense_math_fn,
    .data_funcs = &largemultivector_dense_data_fn,
    .name = "",
    .type_name = "denselarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_dense,
    .basic_size = sizeof(DenseMultivector),
},
  {
    .math_funcs = &largemultivector_blades_math_fn,
    .data_funcs = &largemultivector_blades_data_fn,
    .name = "",
    .type_name = "bladeslarge", 
    .generated = 0,
    .metric = {-2},
    .msize = -1,
    .ntype = MultivectorType_blades,
    .basic_size = sizeof(BladesMultivector),
},
};

// PyMultivectorSubType largemultivector_subtypes_array[3] = {largesparse_subtype,largedense_subtype,largeblades_subtype};


void fill_missing_funcs(void){
     for(Py_ssize_t i = 0; i < 3; i++){
         if(largemultivector_subtypes_array[i].data_funcs->free == NULL)
            largemultivector_subtypes_array[i].data_funcs->free = multivector_subtypes_array[i].data_funcs->free;
         if(largemultivector_subtypes_array[i].data_funcs->init == NULL)
            largemultivector_subtypes_array[i].data_funcs->init = multivector_subtypes_array[i].data_funcs->init;
         if(largemultivector_subtypes_array[i].data_funcs->iter_next == NULL)
            largemultivector_subtypes_array[i].data_funcs->iter_next = multivector_subtypes_array[i].data_funcs->iter_next;
         if(largemultivector_subtypes_array[i].data_funcs->iter_init == NULL)
            largemultivector_subtypes_array[i].data_funcs->iter_init = multivector_subtypes_array[i].data_funcs->iter_init;
         if(largemultivector_subtypes_array[i].data_funcs->cast == NULL)
            largemultivector_subtypes_array[i].data_funcs->cast = multivector_subtypes_array[i].data_funcs->cast;
         if(largemultivector_subtypes_array[i].math_funcs->product == NULL)
            largemultivector_subtypes_array[i].math_funcs->product = multivector_subtypes_array[i].math_funcs->product;
         if(largemultivector_subtypes_array[i].math_funcs->atomic_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->atomic_product = multivector_subtypes_array[i].math_funcs->atomic_product;
         if(largemultivector_subtypes_array[i].math_funcs->ternary_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->ternary_product = multivector_subtypes_array[i].math_funcs->ternary_product;
         if(largemultivector_subtypes_array[i].math_funcs->grade_project == NULL)
            largemultivector_subtypes_array[i].math_funcs->grade_project = multivector_subtypes_array[i].math_funcs->grade_project;
         if(largemultivector_subtypes_array[i].math_funcs->reverse == NULL)
            largemultivector_subtypes_array[i].math_funcs->reverse = multivector_subtypes_array[i].math_funcs->reverse;
         if(largemultivector_subtypes_array[i].math_funcs->add == NULL)
            largemultivector_subtypes_array[i].math_funcs->add = multivector_subtypes_array[i].math_funcs->add;
         if(largemultivector_subtypes_array[i].math_funcs->atomic_add == NULL)
            largemultivector_subtypes_array[i].math_funcs->atomic_add = multivector_subtypes_array[i].math_funcs->atomic_add;
         if(largemultivector_subtypes_array[i].math_funcs->scalar_product == NULL)
            largemultivector_subtypes_array[i].math_funcs->scalar_product = multivector_subtypes_array[i].math_funcs->scalar_product;
         if(largemultivector_subtypes_array[i].math_funcs->scalar_add == NULL)
            largemultivector_subtypes_array[i].math_funcs->scalar_add = multivector_subtypes_array[i].math_funcs->scalar_add;
         if(largemultivector_subtypes_array[i].math_funcs->dual == NULL)
            largemultivector_subtypes_array[i].math_funcs->dual = multivector_subtypes_array[i].math_funcs->dual;
         if(largemultivector_subtypes_array[i].math_funcs->undual == NULL)
            largemultivector_subtypes_array[i].math_funcs->undual = multivector_subtypes_array[i].math_funcs->undual;
     }

     largemultivector_mixed_fn.add = multivector_mixed_fn.add;
     largemultivector_mixed_fn.atomic_add = multivector_mixed_fn.atomic_add;
}
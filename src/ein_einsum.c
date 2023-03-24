#include "ein_einsum.h"

iterator_t init_iterator_t(tensor_strides ts, void **data, size_t sizeof_data){
    iterator_t iter;
    iter.ts = ts;
    iter.data = data;
    iter.sizeof_data = sizeof_data;
    iter.index = (size_t*)malloc(ts.n_subscripts*sizeof(size_t));
    iter.depth = (int*)malloc(ts.n_subscripts*sizeof(int));
    for(size_t k = 0; k < ts.n_subscripts; k++){ // loop over each symbol
        iter.index[k] = 0;
        if(ts.strides[ts.n_tensors-1][k] == 0) // inner iterator sum of products
            iter.depth[k] = 1;
        else // outer iterator
            iter.depth[k] = 0;
    }

    return iter;
}

int general_iterator_t(iterator_t iter, int depth){
    tensor_strides ts = iter.ts;
    void **data = iter.data;
    size_t *index = iter.index;
    int k = -1;
    while(k < (int)(ts.n_subscripts - 1) && iter.depth[++k] != depth); // find first symbol
    if(iter.depth[k] == depth)
        index[k]++;
    else // no symbols at this depth
        return 0; // don't iterate


    if(index[k] < ts.n_strides[k]){ // first symbol not overflown
        for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
            data[j] += ts.strides[j][k]*iter.sizeof_data;
        return 1;  // continue incrementing
    }

    size_t i = k;
    while(i < ts.n_subscripts){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            // go to the beggining
            for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                data[j] -= (ts.n_strides[i]-1)*ts.strides[j][i]*iter.sizeof_data;
            index[i] = 0; // reset symbol i
            if(i == ts.n_subscripts-1) // this is the last symbol
                return 0;
            i++;
            while(i < ts.n_subscripts && iter.depth[i++] != depth); // find next symbol
            if(i >= ts.n_subscripts)// no more symbols, found the last one
                if(iter.depth[i-1] != depth)
                    return 0; // stop incrementing
            index[--i]++; // increment next symbol
            if(index[i] < ts.n_strides[i]){ // symbol i not overflown
                for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                    data[j] += ts.strides[j][i]*iter.sizeof_data;
                return 1; // continue incrementing
            }
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }

    return 1; // continue incrementing
}

int main_einsum_t(
    tensor_multivectors tmvs,
    void* extra,
    operator_functions opfs,
    subscript_struct s,
    tensor* out){

    subscript_shape sp = get_all_subscripts(&s, tmvs.shapes, tmvs.shape_size, tmvs.size); // check symbol-shape consistency
    if(sp.size == 0)
        return 0;

    tensor_multivectors new_tmvs =
        append_out_tensor_t(sp,s.subscripts[s.size_-1],s.size[s.size_-1],tmvs);
    opfs.init(new_tmvs.data[new_tmvs.size-1],new_tmvs.data_size[new_tmvs.size-1]);

    tensor_strides ts = compute_tensor_strides(new_tmvs.shapes,s,sp);

    einsum_sum_prods_t(ts,new_tmvs,opfs,extra);

    out->data = new_tmvs.data[new_tmvs.size-1];
    out->shapes = new_tmvs.shapes[new_tmvs.size-1];
    out->shape_size = new_tmvs.shape_size[new_tmvs.size-1];
    out->data_size = new_tmvs.data_size[new_tmvs.size-1];

    free(new_tmvs.data);
    free(new_tmvs.shapes);
    free(new_tmvs.shape_size);
    free(new_tmvs.data_size);
    free_subscript_shape(sp);
    free_tensor_strides(ts);

    return 1;
}


void einsum_sum_prods_t(
    tensor_strides ts,
    tensor_multivectors tmvs,
    operator_functions opfs,
    void* extra){

    iterator_t iter = init_iterator_t(ts,tmvs.data,tmvs.type_size);
    do{
        sum_of_products_t(tmvs,opfs,extra,iter);
    }while(general_iterator_t(iter,0));

    free(iter.index);
    free(iter.depth);
}

size_t get_nbr_iters_t(iterator_t iter, int depth){
    size_t n_iter = 1;
    for(size_t k = 0; k < iter.ts.n_subscripts; k++)
        if(iter.depth[k] == depth)
            n_iter *= iter.ts.n_strides[k];
    return n_iter;
}


void sum_of_products_t(
    tensor_multivectors tmvs,
    operator_functions opfs,
    void* extra,
    iterator_t iter){

    size_t n_iter = get_nbr_iters_t(iter,1);
    void **sum_mvs = (void*)malloc(n_iter*sizeof(void*));
    size_t j = 0;
    do{
        void *temp = tmvs.data[0];
        void *temp_;

        for(size_t i = 1; i < tmvs.size-1; i++){

            temp_ = temp;
                temp = opfs.product(temp,tmvs.data[i],extra); // compute the geometric product
            if(i > 1){ // free old alocated memory by the product
                opfs.free(temp_,1);
                free(temp_);
            }
        }
        sum_mvs[j] = temp;
        j++;
    }while(general_iterator_t(iter,1));

    void *added = opfs.atomic_add(sum_mvs,j,extra); // adds j multivectors
    void *data_ = tmvs.data[tmvs.size-1];
    void *temp_add = opfs.add(data_,added,extra);
    opfs.assign(data_,temp_add); // copies the contents of temp to the output tensor

    for(size_t i = 0; i < j; i++){
        opfs.free(sum_mvs[i],1);
        free(sum_mvs[i]);
    }

    opfs.free(added,1);
    free(added);
    free(temp_add);
    free(sum_mvs);
}


void free_tensor_strides(tensor_strides ts){
    for(size_t i = 0; i < ts.n_tensors; i++)
        free(ts.strides[i]);
    free(ts.strides);

}

tensor_multivectors append_out_tensor_t(subscript_shape sp, char *symbols, size_t n_symbols, tensor_multivectors tmvs){
   size_t *shape = (size_t*)malloc(n_symbols*sizeof(size_t)); // shape of the output tensor per symbol

    for(size_t i = 0; i < sp.size; i++){
        for(size_t j = 0; j < n_symbols; j++){
            if(symbols[j] == sp.subscripts[i]){
                shape[j] = sp.shape[i];
            }
        }
    }

    size_t size = 1;
    for(size_t i = 0; i < n_symbols; i++)
        size *= shape[i];

    void *out_data = (void*)malloc(size*tmvs.type_size);
    void **data = (void*)malloc((tmvs.size+1)*sizeof(void*));
    size_t **shapes = (size_t**)malloc((tmvs.size+1)*sizeof(size_t*));
    size_t *shape_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t));
    size_t *data_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t));
    for (size_t i = 0; i < tmvs.size; i++) {
        data[i] = tmvs.data[i];
        shapes[i] = tmvs.shapes[i];
        shape_size[i] = tmvs.shape_size[i];
        data_size[i] = tmvs.data_size[i];
    }
    // initialize output tensor to empty
    data[tmvs.size] = out_data;
    shapes[tmvs.size] = shape;
    shape_size[tmvs.size] = n_symbols;
    data_size[tmvs.size] = size;

    tensor_multivectors new_tmvs =
        {data,shapes,shape_size,data_size,tmvs.size+1,tmvs.type_size};
    return new_tmvs;
}

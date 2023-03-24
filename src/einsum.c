#include "einsum.h"

labels parse_subscripts(char *subscripts, size_t size, size_t ndim){
    char *op_labels = (char*)malloc(ndim*sizeof(char));
    size_t dim = 0;
    int n_symbols = 0;
    int ellipses_dim = ndim;
    labels l = {0,op_labels};

    for(size_t i = 0; i < size; i++)
        if(subscripts[i] != '.')
            ellipses_dim--;
    if(ellipses_dim < 0)
        return l;

    int ellipse = -1;
    for(size_t i = 0; i < size; i++){
        if(subscripts[i] == '.'){
            if(i+2 < size){
                if(subscripts[i+1] == '.' && subscripts[i+2] == '.'){
                    for(int i = 0; i < ellipses_dim; i++)
                        op_labels[(ndim-1)-(dim+i)] = 0;
                    dim += ellipses_dim;
                    ellipse = i;
                }
            }
            if(ellipse == -1)
                return l;
        }else{
            if(subscripts[i] > 'z' || subscripts[i] < 'a')
                return l;
            char symbol = subscripts[i];
            int flag = 0;
            for(size_t j = 0; j < i; j++){
                if(symbol == subscripts[j]){
                    if(ellipse != -1){
                        if(j < (size_t)ellipse){
                            op_labels[(ndim-1)-dim] = -(i-j-3 + ellipses_dim);
                            flag = 1;
                            break;
                        }
                    }
                    else{
                        op_labels[(ndim-1)-dim] = (int)j-(int)i;
                        flag = 1;
                        break;
                    }
                }
            }
            if(!flag){ // first occurrence
                op_labels[(ndim-1)-dim] = symbol;
                n_symbols++;
            }
            dim++;
        }
    }
    l.size = ndim;
    return l;
}


symbols parse_args(char *args, size_t size){
    size_t nsubs = 1;
    symbols sym = {NULL,NULL,0};
    for(size_t i = 0; i < size; i++){
        if(args[i]== ',')
            nsubs++;

        if(args[i] == '-'){
            if(i+1 < size){
                if(args[i+1] == '>')
                    nsubs++;
                else{
                    nsubs = 0;
                    break;
                }
            }
        }
    }

    if(nsubs == 0)
        return sym;

    sym.size_ = nsubs;
    sym.size = (size_t*)malloc(nsubs*sizeof(size_t));
    sym.subscripts = (char**)malloc(nsubs*sizeof(char*));

    sym.size[0] = 0;
    size_t sub_j = 0;
    for(size_t i = 0; i < size; i++){
        if(args[i]== ',')
            sym.size[++sub_j] = 0;
        else if(args[i] == '-')
            sym.size[++sub_j] = 0;
        else if(args[i] != '>')
            sym.size[sub_j]++;
    }
    for(size_t i = 0; i < nsubs; i++)
        sym.subscripts[i] = (char*)malloc(sym.size[i]*sizeof(char));

    sub_j = 0;
    size_t k = 0;
    for(size_t i = 0; i < size; i++){
        if(args[i]== ',')
            sub_j++, k = 0;
        else if(args[i] == '-')
            sub_j++, k = 0;
        else if(args[i] != '>')
            sym.subscripts[sub_j][k++] = args[i];
    }
    return sym;
}

symbols parse_all(char *args, size_t size, size_t *ndims, size_t n){
    symbols s = parse_args(args,size);
    symbols z = {NULL,NULL,0};
    if(s.size_ != n) return z;
    z.size_ = s.size_;
    z.size = (size_t*)malloc(z.size_*sizeof(size_t));
    z.subscripts = (char**)malloc(z.size_*sizeof(char*));
    for(size_t i = 0; i < z.size_; i++){
        labels l = parse_subscripts(s.subscripts[i],s.size[i],ndims[i]);
        z.subscripts[i] = l.op_labels;
        z.size[i] = l.size;
    }
    free_symbols(s);
    return z;
}


void free_symbols(symbols s){
    for(size_t i = 0; i < s.size_; i++)
        free(s.subscripts[i]);
    free(s.size);
    free(s.subscripts);
}

int main_einsum_t(
    tensor_multivectors tmvs,
    void* extra,
    operator_functions opfs,
    symbols s,
    tensor* out){

    symbol_shape sp = get_all_symbols(s, tmvs.shapes, tmvs.shape_size, tmvs.size); // check symbol-shape consistency
    if(sp.size == 0)
        return 0;

    tensor_multivectors new_tmvs =
        append_out_tensor_t(sp,s.subscripts[s.size_-1],s.size[s.size_-1],tmvs);
    opfs.init(new_tmvs.data[new_tmvs.size-1],new_tmvs.data_size[new_tmvs.size-1]);

    tensor_strides ts = compute_strides(new_tmvs.shapes,s,sp);

    einsum_sum_prods_t(ts,new_tmvs,opfs,extra);

    out->data = new_tmvs.data[new_tmvs.size-1];
    out->shapes = new_tmvs.shapes[new_tmvs.size-1];
    out->shape_size = new_tmvs.shape_size[new_tmvs.size-1];
    out->data_size = new_tmvs.data_size[new_tmvs.size-1];

    free(new_tmvs.data);
    free(new_tmvs.shapes);
    free(new_tmvs.shape_size);
    free(new_tmvs.data_size);
    free_symbol_shape(sp);
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

void sum_of_products_t(
    tensor_multivectors tmvs,
    operator_functions opfs,
    void* extra,
    iterator_t iter){

    size_t n_iter = get_nbr_inner_iters(iter);
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

tensor_multivectors append_out_tensor_t(symbol_shape sp, char *symbols, size_t n_symbols, tensor_multivectors tmvs){
   size_t *shape = (size_t*)malloc(n_symbols*sizeof(size_t)); // shape of the output tensor per symbol

    for(size_t i = 0; i < sp.size; i++){
        for(size_t j = 0; j < n_symbols; j++){
            if(symbols[j] == sp.symbols[i]){
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


tensor_strides compute_strides(size_t **shapes, symbols sym, symbol_shape sp){
    size_t **strides;
    char *symbols;
    size_t symbols_size = 0;
    tensor_strides ts = {NULL,NULL,0,0};

    symbols_size = sp.size;
    symbols = sp.symbols;
    // strides is of shape [n_tensors,n_symbols]
    strides = (size_t**)malloc(sym.size_*sizeof(size_t*));
    for(size_t i = 0; i < symbols_size; i++){
        strides[i] = (size_t*)malloc(symbols_size*sizeof(size_t));
        for(size_t j = 0; j < symbols_size; j++){
            strides[i][j] = 0;
        }
    }

    for(size_t k = 0; k < sym.size_; k++){
        char *subs = sym.subscripts[k];
        size_t size = sym.size[k];
        size_t stride = 1;

        for(size_t i = 0; i < size; i++){ // loop over each symbol
            // look for symbol in the symbol list
            for(size_t j = 0; j < symbols_size; j++){
                if(subs[i] == symbols[j]){
                    // if symbol found increment stride for that symbol
                    strides[k][j] += stride;
                    break;
                }
            }
            stride *= shapes[k][i];
        }
    }

    ts.n_tensors = sym.size_;
    ts.n_symbols = symbols_size;
    ts.strides = strides;
    ts.n_strides = sp.shape;
    return ts;
}

iterator_t init_iterator_t(tensor_strides ts, void **data, size_t sizeof_data){
    iterator_t iter;
    iter.ts = ts;
    iter.data = data;
    iter.sizeof_data = sizeof_data;
    iter.index = (size_t*)malloc(ts.n_symbols*sizeof(size_t));
    iter.depth = (int*)malloc(ts.n_symbols*sizeof(int));
    for(size_t k = 0; k < ts.n_symbols; k++){ // loop over each symbol
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
    while(k < (int)(ts.n_symbols - 1) && iter.depth[++k] != depth); // find first symbol
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
    while(i < ts.n_symbols){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            // go to the beggining
            for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                data[j] -= (ts.n_strides[i]-1)*ts.strides[j][i]*iter.sizeof_data;
            index[i] = 0; // reset symbol i
            if(i == ts.n_symbols-1) // this is the last symbol
                return 0;
            i++;
            while(i < ts.n_symbols && iter.depth[i++] != depth); // find next symbol
            if(i >= ts.n_symbols)// no more symbols, found the last one
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

size_t get_nbr_inner_iters(iterator_t iter){
    size_t n_iter = 1;
    for(size_t k = 0; k < iter.ts.n_symbols; k++)
        if(iter.ts.strides[iter.ts.n_tensors-1][k] == 0)
            n_iter *= iter.ts.n_strides[k];
    return n_iter;
}

/* This function determines all the different symbols and the shape of each one */
symbol_shape get_all_symbols(symbols sym, size_t **shapes, size_t *shape_size, size_t size){
    size_t symbols_size = 0;
    char *symbols;
    symbol_shape sp = {0,NULL,NULL};
    // check if each tensor symbols are consistent with each corresponding tensor
    // if the number of tensors is different from the symbols
    if(sym.size_-1 != size)
        return sp;

    // number of symbols for each tensor must be equal to the number of dimensions
    for(size_t j = 0; j < size; j++)
        if(sym.size[j] != shape_size[j])
            return sp;

    // initialize list of different symbols
    for(size_t j = 0; j < size; j++)
        symbols_size += sym.size[j];
    symbols = (char*)malloc(symbols_size*sizeof(char));
    size_t *s_shape = (size_t*)malloc(symbols_size*sizeof(size_t)); // shape for each symbol
    symbols_size = 0;

    // determine all the different symbols;
    for(size_t i = 0; i < size; i++){ // loop over each tensor
        for(size_t j = 0; j < sym.size[i]; j++){ // loop over each symbol
            int found = 0;
            for(size_t k = 0; k < symbols_size; k++){ // loop over found symbols
                if(sym.subscripts[i][j] == symbols[k]){ // found symbol
                    found = 1;
                    if(s_shape[k] != shapes[i][j]) // check if shape is consistent
                        return sp;
                }
            }
            if(!found){ // if symbol not in the list
                symbols[symbols_size] = sym.subscripts[i][j]; // add symbol to the list
                s_shape[symbols_size] = shapes[i][j]; // add symbol shape to the list
                symbols_size++;
            }
        }
    }
    sp.size = symbols_size;
    sp.symbols = symbols;
    sp.shape = s_shape;
    return sp;
}

void free_symbol_shape(symbol_shape sp){
    free(sp.symbols);
    free(sp.shape);
}

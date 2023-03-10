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

int main_einsum(graded_tensor_multivectors tmvs, symbols s){
    symbol_shape sp = get_all_symbols(s, tmvs.shapes, tmvs.shape_size, tmvs.size); // check symbol-shape consistency
    if(sp.size == 0)
        return 0;
    // alloc and append output tensor to the end of the tensor list
    graded_tensor_multivectors new_tmvs = append_out_tensor(sp,s.subscripts[s.size_-1],s.size[s.size_-1],tmvs);
    /* free_tensors_holder(tmvs); */
    tensor_strides ts = compute_strides(new_tmvs.shapes,s,sp);

    blades *out_tensor = new_tmvs.data[new_tmvs.size-1];
    einsum_sum_prods(ts,new_tmvs);

    return 1;
}

void free_tensors_holder(graded_tensor_multivectors tmvs){
    free(tmvs.data);
    free(tmvs.shapes);
    free(tmvs.shape_size);
    free(tmvs.data_size);
}

graded_tensor_multivectors append_out_tensor(symbol_shape sp, char *symbols, size_t n_symbols, graded_tensor_multivectors tmvs){
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

    blades *out_data = (blades*)malloc(size*sizeof(blades));
    blades **data = (blades**)malloc((tmvs.size+1)*sizeof(blades*));
    size_t **shapes = (size_t**)malloc((tmvs.size+1)*sizeof(size_t*));
    size_t *shape_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t*));
    size_t *data_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t*));
    for (size_t i = 0; i < tmvs.size; i++) {
        data[i] = tmvs.data[i];
        shapes[i] = tmvs.shapes[i];
        shape_size[i] = tmvs.shape_size[i];
        data_size[i] = tmvs.data_size[i];
    }
    data[tmvs.size] = out_data;
    shapes[tmvs.size] = shape;
    shape_size[tmvs.size] = n_symbols;
    data_size[tmvs.size] = size;
    for(size_t i = 0; i < size; i++)
        out_data[i].grade = NULL;

    graded_tensor_multivectors new_tmvs =
        {data,shapes,shape_size,data_size,tmvs.size+1,tmvs.m,tmvs.gm,tmvs.precision};
    return new_tmvs;
}

tensor_strides compute_strides(size_t **shapes, symbols sym, symbol_shape sp){
    size_t **strides; //= (size_t**)malloc(sym.size_*sizeof(size_t*));
    char *symbols;
    size_t symbols_size = 0;
    size_t *shape;
    tensor_strides ts = {NULL,NULL,0,0};

    symbols_size = sp.size;
    symbols = sp.symbols;
    shape = sp.shape;
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

void einsum_sum_prods(tensor_strides ts, graded_tensor_multivectors tmvs){
    iterator iter = init_iterator(ts,(void**)tmvs.data,sizeof(tmvs.data[0][0]));
    do{
        sum_of_products(tmvs,iter);
    }while(outter_iterator(iter));
}

void einsum_no_sum_prods(tensor_strides ts, graded_tensor_multivectors tmvs){
    iterator iter = init_iterator(ts,(void**)tmvs.data,sizeof(tmvs.data[0][0]));
    do{
        blades temp = *tmvs.data[0];
        blades temp_;
        for(size_t i = 1; i < tmvs.size-1; i++){ // left multiply all the multivectors
            temp_ = temp;
            temp = graded_product_(temp,*tmvs.data[i],tmvs.m,tmvs.gm,tmvs.precision);
            if(i > 1) free_blades(temp_); // free old alocated memory by the product;
        }
        *tmvs.data[tmvs.size-1] = temp;
    }while(outter_iterator(iter));
}

iterator init_iterator(tensor_strides ts, void **data, size_t sizeof_data){
    iterator iter;
    iter.ts = ts;
    iter.data = data;
    iter.sizeof_data = sizeof_data;
    iter.index = (size_t*)malloc(ts.n_symbols*sizeof(size_t));
    for(size_t k = 0; k < ts.n_symbols; k++) // loop over each symbol
        iter.index[k] = 0;
    return iter;
}

// iterates over the output tensor symbols
int outter_iterator(iterator iter){
    tensor_strides ts = iter.ts;
    void **data = iter.data;
    size_t *index = iter.index;
    int k = -1;
    while(k < (int)(ts.n_symbols - 1) && ts.strides[ts.n_tensors-1][++k] == 0); // find first not hidden symbol
    index[k]++;

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
            while(i < ts.n_symbols && ts.strides[ts.n_tensors-1][i++] == 0); // find next not hidden symbol
            if(i >= ts.n_symbols)// no more hidden symbols, found the last one
                if(ts.strides[ts.n_tensors-1][i-1] == 0) // the last symbol is hidden
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

// iterates over the hidden output tensor symbols
int inner_iterator(iterator iter){
    tensor_strides ts = iter.ts;
    void **data = iter.data;
    size_t *index = iter.index;
    int k = -1;
    while(ts.strides[ts.n_tensors-1][++k] != 0); // find first hidden symbol
    index[k]++;

    if(index[k] < ts.n_strides[k]){ // first symbol not overflown
        for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
            data[j] += ts.strides[j][k]*iter.sizeof_data;
        return 1;
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
            while(i < ts.n_symbols && ts.strides[ts.n_tensors-1][i++] != 0); // find next hidden symbol
            if(i >= ts.n_symbols) // no more hidden symbols, found the last one
                if(ts.strides[ts.n_tensors-1][i-1] != 0) // the last symbol is not hidden
                    return 0; // stop incrementing

            index[--i]++; // increment next symbol
            if(index[i] < ts.n_strides[i]){ // symbol i not overflown
                for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                    data[j] += ts.strides[j][i]*iter.sizeof_data;
                return 1;
            }
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }
    return 1; // continue incrementing
}

size_t get_nbr_inner_iters(iterator iter){
    size_t n_iter = 1;

    for(size_t k = 0; k < iter.ts.n_symbols; k++)
        if(iter.ts.strides[iter.ts.n_tensors-1][k] == 0)
            n_iter *= iter.ts.n_strides[k];

    return n_iter;

}

void sum_of_products(graded_tensor_multivectors tmvs, iterator iter){
    size_t n_iter = get_nbr_inner_iters(iter);
    blades **sum_mvs = (blades**)malloc(n_iter*sizeof(blades*));
    size_t j = 0;
    do{
        blades temp = *tmvs.data[0];
        blades temp_;
        for(size_t i = 1; i < tmvs.size-1; i++){ // left multiply the multivectors
            temp_ = temp;
            temp = graded_product_(temp,*tmvs.data[i],tmvs.m,tmvs.gm,tmvs.precision);
            if(i > 1) free_blades(temp_); // free old alocated memory by the product;
        }
        sum_mvs[j] = (blades*)malloc(sizeof(blades));
        *sum_mvs[j] = temp;
        j++;
    }while(inner_iterator(iter));

    blades added = graded_atomic_add_add_(sum_mvs,j,tmvs.gm,tmvs.precision);
    *tmvs.data[tmvs.size-1] = graded_add_add_(*tmvs.data[tmvs.size-1],added,tmvs.gm,tmvs.precision);

    for(size_t i = 0; i < j; i++){
        free_blades(*sum_mvs[i]);
        free(sum_mvs[i]);
    }
    free_blades(added);
    free(sum_mvs);
}

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

graded_tensor vector_matrix_mult(graded_tensor_multivectors tmvs){
    blades *matrix = tmvs.data[0];
    blades *vector = tmvs.data[1];
    blades *vec_out = (blades*)malloc(tmvs.shapes[0][0]*sizeof(blades));
    blades **vec_temp = (blades**)malloc(tmvs.shapes[0][1]*sizeof(blades*));
    graded_tensor out;

    for(size_t i = 0; i < tmvs.shapes[0][0]; i++){
        for(size_t j = 0; j < tmvs.shapes[0][1]; j++){
            blades *temp = (blades*)malloc(sizeof(blades));
            *temp = graded_product_(matrix[i*tmvs.shapes[0][0]+j],vector[j],tmvs.m,tmvs.gm,tmvs.precision);
            vec_temp[j] = temp;
        }
        vec_out[i] = graded_atomic_add_add_(vec_temp,tmvs.shapes[0][1],tmvs.gm,tmvs.precision);
        for(size_t j = 0; j < tmvs.shapes[0][1]; j++){
            free_blades(*vec_temp[j]);
            free(vec_temp[j]);
        }
    }
    free(vec_temp);

    out.data = vec_out;
    out.shapes = (size_t*)malloc(sizeof(size_t));
    out.shapes[0] = tmvs.shapes[0][0];
    out.shape_size = 1;

    return out;
}

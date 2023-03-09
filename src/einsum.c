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
                        op_labels[dim+i] = 0;
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
                            op_labels[dim] = -(i-j-3 + ellipses_dim);
                            flag = 1;
                            break;
                        }
                    }
                    else{
                        op_labels[dim] = (int)j-(int)i;
                        flag = 1;
                        break;
                    }
                }
            }
            if(!flag){ // first occurrence
                op_labels[dim] = symbol;
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


tensor_strides compute_strides(size_t **shapes, size_t *shape_size, size_t size, symbols sym){
    size_t **strides; //= (size_t**)malloc(sym.size_*sizeof(size_t*));
    char *symbols;
    size_t symbols_size = 0;
    tensor_strides ts = {NULL,NULL,0,0};
    // if the number of tensors is different from the symbols
    if(sym.size_ != size)
        return ts;

    // number of symbols for each tensor must be equal to the number of dimensions
    for(size_t j = 0; j < sym.size_; j++)
        if(sym.size[j] != shape_size[j])
            return ts;


    for(size_t j = 0; j < sym.size_; j++)
        symbols_size += sym.size[j];
    symbols = (char*)malloc(symbols_size*sizeof(char));
    symbols_size = 0;


    // determine all the different symbols;
    for(size_t j = 0; j < sym.size_; j++){
        for(size_t i = 0; i < sym.size[j]; j++){
            int found;
            for(size_t k = 0; k < symbols_size; k++){
                if(sym.subscripts[j][i] == symbols[k]){
                    found = 1;
                }
            }
            if(!found){
                symbols[symbols_size] = sym.subscripts[i][j];
                symbols_size++;
            }
        }
    }
    strides = (size_t**)malloc(symbols_size*sizeof(size_t*));
    for(size_t j = 0; j < symbols_size; j++){
        strides[j] = (size_t*)malloc(sym.size_*sizeof(size_t));
        for(size_t i = 0; i < sym.size_; i++){
            strides[j][i] = 0;
        }
    }

    for(size_t j = 0; j < sym.size_; j++){
        size_t stride = 1;
        // run through all symbols
        for(size_t i = 0; i < sym.size[j]; i++){
            for(size_t k = 0; k < symbols_size; k++){
                if(sym.subscripts[j][i] == symbols[k]){
                    // if symbols are repeated also increment stride by corrensponding shape
                    strides[i][j] += stride;
                    break;
                }
            }
            stride *= shapes[j][i];
        }
    }
    ts.n_tensors = sym.size_;
    ts.n_symbols = symbols_size;
    ts.strides = strides;
    return ts;
}

void einsum_sum_prods(tensor_strides ts, graded_tensor_multivectors tmvs){

    size_t size = ts.n_tensors - 1;

    size_t *index = (size_t*)malloc(ts.n_symbols*sizeof(size_t));

    for(size_t k = 0; k < ts.n_symbols; k++) // loop over each symbol
        index[k] = 0;


    while(iterator((void**)tmvs.data,ts.strides,ts.n_tensors,ts.n_strides,index,ts.n_symbols)){
        sum_of_products(tmvs,ts.strides[],size,ts.n_strides[]);
    }


    for(size_t k = 0; k < ts.n_symbols; k++){ // loop over each symbol
        if(ts.strides[k][size] != 0){ // symbol is in the output
            for(size_t l = 0; l < ts.n_strides[k]; l++){
                // loop over the symbols that are not in the output
                for(size_t j = 0; j < ts.n_symbols; j++){ // loop over each symbol
                    if(ts.strides[j][size] == 0)  // symbol not in the output
                        sum_of_products(tmvs,ts.strides[j],size,ts.n_strides[j]);
                }

                for(size_t i = 0; i < ts.n_tensors; i++) // loop over each tensor
                    tmvs.data[i] += ts.strides[k][i];
            }
        }
    }
}


// iterates over the output tensor symbols
int outter_iterator(tensor_strides ts, void **data, size_t *index){
    int k = -1;
    while(ts.strides[++k][ts.n_tensors-1] == 0); // find first not hidden symbol
    index[k]++;

    if(index[k] < ts.n_strides[k]) // first symbol not overflown
        for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
            data[j] += ts.strides[k][j];

    size_t i = k;
    while(i < ts.n_symbols){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            index[i++] = 0; // reset symbol i
            while(i < ts.n_symbols - 1 && ts.strides[i++][ts.n_tensors-1] == 0); // find next not hidden symbol
            if(i >= ts.n_symbols - 1) // no more hidden symbols, found the last one
                return 0; // stop incrementing
            index[i]++; // increment next symbol
            if(index[i] < ts.n_strides[i]) // symbol i not overflown
                for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                    data[j] += ts.strides[i][j];
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }

    return 1; // continue incrementing
}

// iterates over the hidden output tensor symbols
int inner_iterator(tensor_strides ts, void **data, size_t *index){
    int k = -1;
    while(ts.strides[++k][ts.n_tensors-1] != 0); // find first hidden symbol
    index[k]++;

    if(index[k] < ts.n_strides[k]) // first symbol not overflown
        for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
            data[j] += ts.strides[k][j];

    size_t i = k;
    while(i < ts.n_symbols){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            index[i++] = 0; // reset symbol i
            while(i < ts.n_symbols - 1 && ts.strides[i++][ts.n_tensors-1] != 0); // find next hidden symbol
            if(i >= ts.n_symbols - 1) // no more hidden symbols, found the last one
                return 0; // stop incrementing
            index[i]++; // increment next symbol
            if(index[i] < ts.n_strides[i]) // symbol i not overflown
                for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                    data[j] += ts.strides[i][j];
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }

    return 1; // continue incrementing
}


void einsum_no_sum_prods(tensor_strides ts, graded_tensor_multivectors tmvs){

    size_t size = ts.n_tensors - 1;
    for(size_t j = 0; j < ts.n_symbols; j++){ // loop over each symbol
        for(size_t k = 0; k < ts.n_strides[j]; k++){
            blades temp = *tmvs.data[0];
            blades temp_;
            for(size_t i = 1; i < size; i++){ // left multiply all the multivectors
                temp_ = temp;
                temp = graded_product_(temp,*tmvs.data[i],tmvs.m,tmvs.gm,tmvs.precision);
                if(i > 1) free_blades(temp_); // free old alocated memory by the product;
            }
            *tmvs.data[size] = temp;

            for(size_t i = 0; i < ts.n_tensors; i++) // loop over each tensor
                tmvs.data[i] += ts.strides[j][i];
        }
    }
}

void sum_of_products(graded_tensor_multivectors tmvs, size_t *strides, size_t size, size_t n_strides){

    blades **sum_mvs = (blades**)malloc(n_strides*sizeof(blades*));
    for(size_t j = 0; j < n_strides; j++){
        blades temp = *tmvs.data[0];
        blades temp_;
        for(size_t i = 1; i < size; i++){ // left multiply the multivectors
            temp_ = temp;
            temp = graded_product_(temp,*tmvs.data[i],tmvs.m,tmvs.gm,tmvs.precision);
            if(i > 1) free_blades(temp_); // free old alocated memory by the product;
        }

        sum_mvs[j] = (blades*)malloc(sizeof(blades));
        *sum_mvs[j] = temp;
        for(size_t i = 0; i < size; i++)
            tmvs.data[i] += strides[i];
    }

    blades added = graded_atomic_add_add_(sum_mvs,n_strides,tmvs.gm,tmvs.precision);
    *tmvs.data[size] = graded_add_add_(*tmvs.data[size],added,tmvs.gm,tmvs.precision);
    for(size_t j = 0; j < n_strides; j++){
        free_blades(*sum_mvs[j]);
        free(sum_mvs[j]);
    }
    free_blades(added);
    free(sum_mvs);

}

graded_tensor initialize_graded_out_tensor(symbols sym, size_t **shapes, size_t *shape_size, size_t size){
    size_t out_sub_size = sym.size[sym.size_-1];
    char *out_sub = sym.subscripts[sym.size_-1];
    size_t *out_shapes = (size_t*)malloc(out_sub_size*sizeof(size_t));
    graded_tensor out_tensor = {NULL,NULL,0};

    for(size_t i = 0; i < out_sub_size; i++)
        out_shapes[i] = 0;

    for(size_t i = 0; i < out_sub_size; i++){
        for(size_t j = 0; j < sym.size_-1; j++){
            for(size_t k = 0; k < sym.size[j]; k++){
                if(out_sub[i] == sym.subscripts[j][k]){
                    if(out_shapes[i] != 0){
                        out_shapes[i] = shapes[j][k];
                    }else{
                        // check if shapes are equal
                        if(out_shapes[i] != shapes[j][k]){
                            return out_tensor;
                        }
                    }
                    break;
                }
            }
        }
    }

    size_t out_size = 1;

    for(size_t i = 0; i < out_sub_size; i++)
        out_size *= out_shapes[i];

    out_tensor.data = (blades*)malloc(out_size*sizeof(blades));
    out_tensor.shapes = out_shapes;
    out_tensor.shape_size = out_sub_size;

    return out_tensor;

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

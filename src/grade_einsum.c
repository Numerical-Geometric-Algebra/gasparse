#include "grade_einsum.h"



int general_iterator(iterator iter, int depth){
    tensor_strides ts = iter.ts;
    grades_strides gs = iter.gs;
    void **data = iter.data;
    int **grades = iter.grades;
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
        for(size_t j = 0; j < gs.n_operations; j++) // loop over operations
            grades[j][k] += gs.strides[j][k]; // update grades
        return 1;  // continue incrementing
    }

    size_t i = k;
    while(i < ts.n_symbols){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            // go to the beggining
            for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                data[j] -= (ts.n_strides[i]-1)*ts.strides[j][i]*iter.sizeof_data;
            for(size_t j = 0; j < gs.n_operations; j++) // loop over operations
                grades[j][i] -= (gs.n_strides[i]-1)*gs.strides[j][i]; // reset grades
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
                for(size_t j = 0; j < gs.n_operations; j++) // loop over operations
                    grades[j][i] += gs.strides[j][i]; // update grades
                return 1; // continue incrementing
            }
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }

    return 1; // continue incrementing
}

strides_info compute_strides_info(expression_struct es){
    expression_struct *temp_es = &es;
    // find first expression to be computed
    while(!is_symbol(temp_es->left_sub.symbol) || !is_symbol(temp_es->right_sub.symbol)){
        if(temp_es->left != NULL && !temp_es->left->visited)
            temp_es = temp_es->left;
        else if(temp_es->right != NULL && !temp_es->right->visited)
             temp_es = temp_es->right;
        else if(temp_es->up != NULL)
            temp_es->visited = 1, temp_es = temp_es->up;
        else break;
    }

    // data -> pointer to tensor
    // operator -> string with the operation information
    // operation order -> left or right multiplication
    // subscripts -> array of char
    // grades -> array of char
    /*
    ** Need to know
    ** - all symbols
    ** - all subscripts
    ** - all grades
     */



    strides_info si;
    return si;
}


grades_strides compute_strides(size_t **shapes, symbols sym, symbol_shape sp){
    size_t **strides;
    char *symbols;
    size_t symbols_size = 0;
    grades_strides gs = {NULL,NULL,0,0};

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

    gs.n_operations = sym.size_;
    gs.n_symbols = symbols_size;
    gs.strides = strides;
    gs.n_strides = sp.shape;
    return gs;
}

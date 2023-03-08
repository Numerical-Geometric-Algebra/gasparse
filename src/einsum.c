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

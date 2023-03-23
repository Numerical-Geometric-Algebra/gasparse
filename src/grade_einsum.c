#include "grade_einsum.h"


int general_iterator(iterator iter, int depth){
    tensor_strides ts = iter.ts;
    grades_strides gs = iter.gs;
    void ***data = iter.data;
    size_t **grades = iter.grades;
    size_t *index = iter.index;
    int k = -1;
    while(k < (int)(ts.n_subscripts - 1) && iter.depth[++k] != depth); // find first symbol
    if(iter.depth[k] == depth)
        index[k]++;
    else // no symbols at this depth
        return 0; // don't iterate


    if(index[k] < ts.n_strides[k]){ // first symbol not overflown
        for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
            *(data[j]) += ts.strides[j][k]*iter.sizeof_data;
        for(size_t j = 0; j < gs.n_grades; j++) // loop over grades
            *(grades[j]) += gs.strides[j][k]; // update grades
        return 1;  // continue incrementing
    }

    size_t i = k;
    while(i < ts.n_subscripts){ // loop over all symbols
        if(index[i] >= ts.n_strides[i]){ // symbol i overflown
            // go to the beggining
            for(size_t j = 0; j < ts.n_tensors; j++) // loop over tensors
                *(data[j]) -= (ts.n_strides[i]-1)*ts.strides[j][i]*iter.sizeof_data;
            for(size_t j = 0; j < gs.n_grades; j++) // loop over operations
                *(grades[j]) -= gs.max_grade*gs.strides[j][i]; // reset grades
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
                    *(data[j]) += ts.strides[j][i]*iter.sizeof_data;
                for(size_t j = 0; j < gs.n_grades; j++) // loop over grades
                    *(grades[j]) += gs.strides[j][i]; // update grades
                return 1; // continue incrementing
            }
        }else{ // no more symbols to increment
            return 1; // continue incrementing
        }
    }

    return 1; // continue incrementing
}

int set_symbols(char *symbols, size_t *size, char sub_symbol, int *pos){
    *pos = -1;
    if(*size > MAX_SYMBOLS)
        return -1;
    if(sub_symbol != '\0'){
        int flag = 1;
        for(size_t i = 0; i < *size; i++){
            if(symbols[i] == sub_symbol){
                flag = 0;
                *pos = i;
                break;
            }
        }

        if(flag){
            symbols[*size] = sub_symbol;
            (*size)++;
        }
        return flag;
    }else
        return -1;
}

int next_symbol(symbol_iterator iter, char symbol){
   if(symbol == '\0')
        return -1;
   size_t k = 0;
   for(size_t i = 0; i < iter.size; i++){
        if(symbol == iter.symbols[k]){
            int index =  k+iter.pos[i];
            iter.pos[i]++;
            return index;
        }
        k += iter.repeated[i];
   }
   return -1;
}

void reset_symbol_iterator(symbol_iterator iter){
    for(size_t i = 0; i < iter.size; i++)
        iter.pos[i] = 0;
}


// visited variable to not having to reset all the visited values
// as long as all values of visited in the tree are set in a loop
int iterate_expression_symbols(expression_struct **es, int visited){
    if(*es == NULL)
        return 0;

    // check if node was visited
    if((*es)->visited == visited)
        return 0;

    // go down the tree until a symbol is found
    while((*es)->left_sub.symbol == '\0' && (*es)->right_sub.symbol == '\0'){
        if((*es)->left != NULL && (*es)->left->visited != visited)
            *es = (*es)->left;
        else if((*es)->right != NULL && (*es)->right->visited != visited)
            (*es) = (*es)->right;
        else if((*es)->up != NULL)
            (*es)->visited = visited, (*es) = (*es)->up;
        else return 0;
    }

    // set this node to visited
    if((*es)->left_sub.symbol != '\0' || (*es)->right_sub.symbol != '\0')
        (*es)->visited = visited;

    return 1;
}

int iterate_expression_operations(operation_tree **ot, expression_struct **es, int visited){
    if(*es == NULL)
        return 0;

    // if reached bottom go up and then one down
    if((*es)->left == NULL && (*es)->right == NULL){
        // try to go up and then right until found not visited
        do{
            (*es)->visited = visited; // set node to visited
            (*ot)->visited = visited;
            (*es) = (*es)->up; // go up
            (*ot) = (*ot)->up_op;
            if((*es) == NULL)
                return 0;
            if((*ot) == NULL)
                return 0;
            if((*es)->right != NULL && (*es)->right->visited != visited){
                (*es) = (*es)->right; // then try go right

                // add right node to the operattion tree
                if((*ot)->right_op == NULL){
                    (*ot)->right_op = (operation_tree*)malloc(sizeof(operation_tree));
                    initialize_operation_tree((*ot)->right_op);
                }
                (*ot)->right_op->up_op = (*ot);
                (*ot) = (*ot)->right_op;
            }
            else{
                (*es)->visited = visited; // if can't go right go up again
                (*ot)->visited = visited;
            }
        }while((*es)->visited == visited);

        return 1;
    }

    if((*es)->left != NULL && (*es)->left->visited != visited){ // try to go left
        *es = (*es)->left;

        // add left node to the operation tree
        if((*ot)->left_op == NULL){
            (*ot)->left_op = (operation_tree*)malloc(sizeof(operation_tree));
            initialize_operation_tree((*ot)->left_op);
        }
        (*ot)->left_op->up_op = (*ot);
        (*ot) = (*ot)->left_op;
    }
    else if((*es)->right != NULL && (*es)->right->visited != visited){
        (*es) = (*es)->right;
        if((*ot)->right_op == NULL){
            (*ot)->right_op = (operation_tree*)malloc(sizeof(operation_tree));
            initialize_operation_tree((*ot)->right_op);
        }
        (*ot)->right_op->up_op = (*ot);
        (*ot) = (*ot)->right_op;
    }
    else if((*es)->up != NULL){
        (*es)->visited = visited, (*es) = (*es)->up;
        (*ot)->visited = visited, (*ot) = (*ot)->up_op;
    }
    else return 0;

    return 1;
}


int iterate_expression(expression_struct **es, int visited){
    if(*es == NULL)
        return 0;

    // if reached bottom go up and then one down
    if((*es)->left == NULL && (*es)->right == NULL){
        // try to go up and then right until found not visited
        do{
            (*es)->visited = visited; // set node to visited
            (*es) = (*es)->up; // go up
            if((*es) == NULL)
                return 0;
            if((*es)->right != NULL && (*es)->right->visited != visited){
                (*es) = (*es)->right; // then try go right

            }
            else{
                (*es)->visited = visited; // if can't go right go up again
            }
        }while((*es)->visited == visited);

        return 1;
    }

    if((*es)->left != NULL && (*es)->left->visited != visited){ // try to go left
        *es = (*es)->left;
    }
    else if((*es)->right != NULL && (*es)->right->visited != visited){
        (*es) = (*es)->right;
    }
    else if((*es)->up != NULL){
        (*es)->visited = visited, (*es) = (*es)->up;
    }
    else return 0;

    return 1;
}


void set_operator_index(char *operators, size_t size, expression_struct *es){
    int visited = !es->visited;
    do{
        int flag = 0;
        for(size_t i = 0; i < size; i++){
            if(operators[i] == es->operator){
                es->operator = i;
                flag = 1;
                break;
            }
        }
        if(!flag) es->operator = 0; // if the operator is not found set the index to 0
    }while(iterate_expression(&es,visited));
}

tensor_multivectors repeat_tensors(char *symbols, size_t *repeated, tensor_multivectors tmvs){
    tensor_multivectors new_tmvs;
    size_t rep_size = 0;
    for(size_t i = 0; i < tmvs.size; i++)
        rep_size += repeated[i];

    new_tmvs.data = (void**)malloc(rep_size*sizeof(void*));
    new_tmvs.shapes = (size_t**)malloc(rep_size*sizeof(size_t*));
    new_tmvs.shape_size = (size_t*)malloc(rep_size*sizeof(size_t));
    new_tmvs.data_size = (size_t*)malloc(rep_size*sizeof(size_t));
    // associate data by alphabetical order
    char smallest = 'a';
    size_t k = 0;
    size_t *index = (size_t*)malloc(tmvs.size*sizeof(size_t));
    while(is_symbol(smallest)){
        for(size_t i = 0; i < tmvs.size; i++){
            if(symbols[i] == smallest){
                index[i] = k;
                k++;
                break;
            }
        }
        smallest++;
    }

    k = 0;
    for(size_t i = 0; i < tmvs.size; i++){
        for(size_t j = 0; j < repeated[i]; j++){
            new_tmvs.data[k] = tmvs.data[index[i]];
            new_tmvs.shapes[k] = tmvs.shapes[index[i]];
            new_tmvs.shape_size[k] = tmvs.shape_size[index[i]];
            new_tmvs.data_size[k] = tmvs.data_size[index[i]];
            k++;
        }
    }

    new_tmvs.size = rep_size;
    new_tmvs.type_size = tmvs.type_size;
    free(index);
    return new_tmvs;

}

void free_tensor_multivectors(tensor_multivectors tmvs){
    for(size_t i = 0; i < tmvs.size; i++){
        free(tmvs.data[i]);
        free(tmvs.shapes[i]);
    }free(tmvs.shapes); free(tmvs.data);
    free(tmvs.shape_size);
    free(tmvs.data_size);
}

int initialize_einsum(expression_struct *es,
                      char *output_subscripts,
                      tensor_multivectors tmvs,
                      size_t max_grade,
                      iterator *iter,
                      operation_tree **op_tree__,
                      tensor *out){

    char symbols[MAX_SYMBOLS];
    size_t subscripts_length[MAX_SYMBOLS];
    size_t sym_size = 0;
    size_t repeated[MAX_SYMBOLS];
    expression_struct *temp_es = es;
    int pos = -1;

    char **subscripts;
    char *repeated_symbols;
    size_t rep_size = 0;
    size_t *symbol_pos;
    /* order *op_order; */
    int index;

    for(size_t i = 0; i < MAX_SYMBOLS; i++){
        subscripts_length[i] = 0;
        repeated[i] = 0;
    }

    // get list of all symbols
    int visited = !temp_es->visited;
    while(iterate_expression_symbols(&temp_es,visited)){
        // add left symbol to list and check subscripts
        int flag = set_symbols(symbols,&sym_size, temp_es->left_sub.symbol,&pos);
        if(flag == 0){
            // not a new symbol
            if(subscripts_length[pos] != strlen(temp_es->left_sub.subscripts)){
                // same symbol but different shape (number of subscripts)
                return -1;
            }else{
                repeated[pos]++; // repeated symbol
            }
        }
        else if(flag == 1){ // add subscripts length to table
            subscripts_length[sym_size-1] = strlen(temp_es->left_sub.subscripts);
            repeated[sym_size-1] = 1;
        }


        flag = set_symbols(symbols,&sym_size, temp_es->right_sub.symbol,&pos);
        // add right symbol to list and check subscripts
        if(flag == 0){
            // not a new symbol
            if(subscripts_length[pos] != strlen(temp_es->right_sub.subscripts)){
                // same symbol but different shape (number of subscripts)
                return -1;
            }else{
                repeated[pos]++; // repeated symbol
            }
        }else if(flag == 1){ // add subscripts length to table
            subscripts_length[sym_size-1] = strlen(temp_es->right_sub.subscripts);
            repeated[sym_size-1] = 1;
        }
    }
    // the number of input tensors should be the same as the number of different symbols
    if(tmvs.size != sym_size)
        return -1;

    // if there are repeated symbols copy corresponding repeated tensors to a new tensor array
    tensor_multivectors new_tmvs = repeat_tensors(symbols,repeated,tmvs);

    for(size_t i = 0; i < sym_size; i++)
        rep_size += repeated[i];

    /* op_order = (order*)malloc(rep_size*sizeof(order)); */
    subscripts = (char**)malloc((rep_size+1)*sizeof(char*)); // +1 to include output subscripts
    repeated_symbols = (char*)malloc(rep_size*sizeof(char));
    symbol_pos = (size_t*)malloc(sym_size*sizeof(size_t));
    size_t k = 0;
    for(size_t i = 0; i < sym_size; i++){
        for(size_t j = 0; j < repeated[i]; j++){
            subscripts[k] = (char*)malloc((subscripts_length[i] + 1)*sizeof(char));
            repeated_symbols[k] = symbols[i];
            k++;
        }
    }

    symbol_iterator sym_iter = {repeated,symbol_pos,repeated_symbols,sym_size};
    size_t ***grades = (size_t***)malloc(rep_size*3*sizeof(size_t**));
    size_t **grade_strides = (size_t**)malloc(rep_size*3*sizeof(size_t*));
    char **grade_subscripts = (char**)malloc(rep_size*3*sizeof(char*));
    size_t size = 0;

    void ***data = (void***)malloc((rep_size+1)*sizeof(void**)); // include output tensor
    operation_tree *op_tree = (operation_tree*)malloc(sizeof(operation_tree));
    initialize_operation_tree(op_tree);
    operation_tree *parent_node = op_tree;

    for(size_t i = 0; i < rep_size*3; i++){
        grades[i] = NULL;
        grade_strides[i] = NULL;
        grade_subscripts[i] = NULL;
    }

    reset_symbol_iterator(sym_iter);
    temp_es = es;

    // determine the output grades
    // the loop only determines sub childs output grades
    index = set_grades(temp_es->grades,
                       grades,
                       grade_strides,
                       grade_subscripts,
                       &size);
    if(index != -1){
        op_tree->grades.out = grades[index];
        op_tree->grades.out_size = strlen(grade_subscripts[index]);
    }

    visited = !temp_es->visited;
    do{
        // add data to the tree and to the list
        index = next_symbol(sym_iter,temp_es->left_sub.symbol);
        if(index != -1){
            data[index] = (void**)malloc(sizeof(void*));
            *(data[index]) = new_tmvs.data[index];
            op_tree->left = data[index];
            strcpy(subscripts[index], temp_es->left_sub.subscripts);
        }

        index = next_symbol(sym_iter,temp_es->right_sub.symbol);
        if(index != -1){
            data[index] = (void**)malloc(sizeof(void*));
            *(data[index]) = new_tmvs.data[index];
            op_tree->right = data[index];
            strcpy(subscripts[index], temp_es->right_sub.subscripts);
        }

        if(temp_es->left != NULL){
            index = set_grades(temp_es->left->grades,
                               grades,
                               grade_strides,
                               grade_subscripts,
                               &size);
            if(index != -1){
                op_tree->left_op = (operation_tree*)malloc(sizeof(operation_tree));
                initialize_operation_tree(op_tree->left_op);
                op_tree->left_op->grades.out = op_tree->grades.left = grades[index];
                op_tree->left_op->grades.out_size = op_tree->grades.left_size = strlen(grade_subscripts[index]);
            }
        }else if(temp_es->left_sub.symbol != '\0'){
            index = set_grades(temp_es->left_sub.grades,
                               grades,
                               grade_strides,
                               grade_subscripts,
                               &size);
            if(index != -1){
                op_tree->grades.left = grades[index];
                op_tree->grades.left_size = strlen(grade_subscripts[index]);
            }
        }


        if(temp_es->right != NULL){
            index = set_grades(temp_es->right->grades,
                               grades,
                               grade_strides,
                               grade_subscripts,
                               &size);
            if(index != -1){
                op_tree->right_op = (operation_tree*)malloc(sizeof(operation_tree));
                initialize_operation_tree(op_tree->right_op);
                op_tree->right_op->grades.out = op_tree->grades.right = grades[index];
                op_tree->right_op->grades.out_size = op_tree->grades.right_size = strlen(grade_subscripts[index]);
            }
        }else if(temp_es->right_sub.symbol != '\0'){
            index = set_grades(temp_es->right_sub.grades,
                               grades,
                               grade_strides,
                               grade_subscripts,
                               &size);
            if(index != -1){
                op_tree->grades.right = grades[index];
                op_tree->grades.right_size = strlen(grade_subscripts[index]);
            }
        }

        op_tree->grades.operator = temp_es->operator;
    }while(iterate_expression_operations(&op_tree,&temp_es,visited));

    // add output subscripts to the end of the list
    size_t len = strlen(output_subscripts) + 1;
    subscripts[rep_size] = (char*)malloc(len*sizeof(char));
    strcpy(subscripts[rep_size],output_subscripts);
    subscript_struct sub = {NULL,subscripts,rep_size + 1};
    subscript_shape sps = get_all_subscripts(&sub,new_tmvs.shapes,new_tmvs.shape_size,new_tmvs.size);
    char *diff_grade_subscripts = get_grade_subscripts(grade_subscripts,size);
    tensor_multivectors tmvs_out = append_out_tensor(new_tmvs,
                                                     &sps,output_subscripts,
                                                     diff_grade_subscripts,max_grade,out);
    free(new_tmvs.data);
    free(new_tmvs.shapes);
    free(new_tmvs.shape_size);
    free(new_tmvs.data_size);

    tensor_strides ts = compute_tensor_strides(tmvs_out.shapes, sub, sps);

    size_t **grades_list = NULL;
    size_t **grade_strides_list = NULL;
    size_t total_size = 0; // the number of all the different grade subscripts

    // compute grades strides
    if(size > 0){
        for(size_t i = 0; i < size; i++)
            total_size += strlen(grade_subscripts[i]);

        grades_list = (size_t**)malloc(total_size*sizeof(size_t*));
        grade_strides_list = (size_t**)malloc(total_size*sizeof(size_t*));


        k = 0;
        for(size_t i = 0; i < size; i++){
            size_t len  = strlen(grade_subscripts[i]);
            for(size_t j = 0; j < len; j++){
                grades_list[k] = grades[i][j];
                grade_strides_list[k] = (size_t*)malloc(sps.size*sizeof(size_t));
                // check if subscript is in the list of different subscripts and determine the strides
                for(size_t l = 0; l < sps.size; l++){
                    if(grade_subscripts[i][j] == sps.subscripts[l])
                        grade_strides_list[k][l] = 1;
                    else
                        grade_strides_list[k][l] = 0;
                }
                k++;
            }
        }
    }

    data[rep_size] = (void**)malloc(sizeof(void*)); // alloc memory for the output tensor
    *(data[rep_size]) = out->data;
    iterator iter__ =
        {ts,{grade_strides_list,max_grade,total_size,sps.size},
         data,grades_list,NULL,NULL,
         ts.n_tensors-1,tmvs.type_size};

    initialize_iterator(&iter__);
    *iter = iter__;
    *op_tree__ = parent_node;

    free_subscript_shape(sps);
    free_subscript_struct(sub);
    free(sym_iter.pos);
    free(sym_iter.symbols);
    free(tmvs_out.data);
    free(tmvs_out.shapes);
    free(tmvs_out.shape_size);
    free(tmvs_out.data_size);

    if(diff_grade_subscripts != NULL)
        free(diff_grade_subscripts);

    for(size_t i = 0; i < rep_size*3; i++){
        if(grade_strides[i] != NULL)
            free(grade_strides[i]);
        if(grade_subscripts[i] != NULL)
            free(grade_subscripts[i]);
    }
    free(grades);
    free(grade_strides);
    free(grade_subscripts);

    return 1;
}

void free_subscript_struct(subscript_struct sub){
    for(size_t i = 0; i < sub.size_; i++){
        free(sub.subscripts[i]);
    }free(sub.subscripts);
    free(sub.size);
}

void free_subscript_shape(subscript_shape sp){
    free(sp.subscripts);
    free(sp.shape);
}

void free_symbol_iterator(symbol_iterator iter){
    free(iter.repeated);
    free(iter.pos);
    free(iter.symbols);
}

void free_operation_tree_recursive(operation_tree *op){
    if(op->right_op != NULL)
        free_operation_tree_recursive(op->right_op);
    if(op->left_op != NULL)
        free_operation_tree_recursive(op->left_op);

    if(op->grades.left != NULL) free(op->grades.left);
    if(op->grades.right != NULL) free(op->grades.right);
    if(op->grades.out != NULL) free(op->grades.out);

    free(op);
}


int main_einsum(tensor_multivectors tmvs,
                void *extra,
                operator_functions opfs,
                expression_struct *es,
                char *output_subscripts,
                size_t max_grade,
                tensor *out){

    iterator iter;
    operation_tree *op_tree;

    initialize_einsum(es,output_subscripts,tmvs,max_grade,&iter,&op_tree,out);
    opfs.init(out->data,out->data_size); // initialize output tensor to 0
    einsum_sum_prods(opfs,extra,iter,op_tree);
    free_iterator(iter);
    free_operation_tree_recursive(op_tree);

    return 1;
}

void free_iterator(iterator iter){
    free_tensor_strides(iter.ts);
    free_grades_strides(iter.gs);
    for(size_t i = 0; i < iter.ts.n_tensors; i++){
        free(iter.data[i]);
    }free(iter.data);

    if(iter.grades != NULL){
        for(size_t i = 0; i < iter.gs.n_grades; i++){
            if(iter.grades[i] != NULL)
                free(iter.grades[i]);
        }free(iter.grades);
    }

    free(iter.index);
    free(iter.depth);
}


void initialize_operation_tree(operation_tree *tree){
    tree->left = NULL;
    tree->right = NULL;
    tree->left_op = NULL;
    tree->right_op = NULL;
    tree->up_op = NULL;
    tree->grades.left = NULL;
    tree->grades.right = NULL;
    tree->grades.out = NULL;
    tree->grades.operator = '\0';
    tree->grades.left_size = 0;
    tree->grades.right_size = 0;
    tree->grades.out_size = 0;
    tree->result = NULL;
    tree->visited = 0;
}

size_t get_nbr_iters(iterator iter, int depth){
    size_t n_iter = 1;
    for(size_t k = 0; k < iter.ts.n_subscripts; k++)
        if(iter.depth[k] == depth)
            n_iter *= iter.ts.n_strides[k];
    return n_iter;
}

void sum_of_products(
    operator_functions opfs,
    void* extra, iterator iter, operation_tree *op_tree){

    int n_iter = get_nbr_iters(iter,1);

    void **sum_mvs = (void*)malloc(n_iter*sizeof(void*));
    int *free_result = (int*)malloc(n_iter*sizeof(int));
    size_t j = 0;
    do{
        sum_mvs[j] = recursive_products(opfs,extra,op_tree,free_result+j);
        j++;
    }while(general_iterator(iter,1));

    void *added = opfs.atomic_add(sum_mvs,n_iter,extra); // adds n_iter multivectors
    void *data = *iter.data[iter.out_index];
    void *temp_add = opfs.add(data,added,extra);
    opfs.assign(data,temp_add);
    for(int i = 0; i < n_iter; i++){
        if(free_result[i]){
            opfs.free(sum_mvs[i],1);
            free(sum_mvs[i]);
        }
    }

    opfs.free(added,1);
    free(added);
    free(temp_add);
    free(free_result);
    free(sum_mvs);

}



void einsum_sum_prods(
    operator_functions opfs,
    void* extra, iterator iter, operation_tree *op_tree){

    do{
        sum_of_products(opfs,extra,iter,op_tree);
    }while(general_iterator(iter,0));

}

void *recursive_products(operator_functions opfs, void *extra, operation_tree *op_tree, int *free_result){
    void *right,*left,*result;
    int right_free = 0, left_free = 0;
    int free_result_left = 0, free_result_right = 0;
    *free_result = 0;

    if(op_tree->left_op == NULL)
        if(op_tree->left != NULL) left = *op_tree->left;
        else left = NULL;
    else left = recursive_products(opfs,extra,op_tree->left_op,&free_result_left), left_free = 1;

    if(op_tree->right_op == NULL)
        if(op_tree->right != NULL) right = *op_tree->right;
        else right = NULL;
    else right = recursive_products(opfs,extra,op_tree->right_op,&free_result_right), right_free = 1;

    if(right == NULL && left == NULL) return NULL;

    if(right != NULL && left != NULL)
        result = opfs.product(left,right,extra,op_tree->grades), *free_result = 1; // only free the result when the product is computed
    else if(right != NULL)
        result = right, *free_result = right_free, right_free = 0; // I don't want to free the result here
    else
        result = left, *free_result = left_free, left_free = 0; // I don't want to free the result here

    if(free_result_left && left_free && left != NULL) opfs.free(left,1), free(left);
    if(free_result_right && right_free && right != NULL) opfs.free(right,1), free(right);

    return result;
}

void initialize_iterator(iterator *iter){
    iter->index = (size_t*)malloc(iter->ts.n_subscripts*sizeof(size_t));
    iter->depth = (int*)malloc(iter->ts.n_subscripts*sizeof(int));
    for(size_t k = 0; k < iter->ts.n_subscripts; k++){ // loop over each symbol
        iter->index[k] = 0;
        if(iter->ts.strides[iter->ts.n_tensors-1][k] == 0) // inner iterator sum of products
            iter->depth[k] = 1;
        else // outer iterator
            iter->depth[k] = 0;
    }
}

tensor_multivectors append_out_tensor(tensor_multivectors tmvs,
                                      subscript_shape *sps, char *output_subscripts, // subscripts information
                                      char *grade_subscript, size_t max_grade, /* grade subscripts information */
                                      tensor *out){

    tensor_multivectors new_tmvs = {NULL,NULL,NULL,NULL,0,0};
    size_t grade_len = 0;
    if(grade_subscript != NULL)
        grade_len = strlen(grade_subscript);
    size_t size = sps->size + grade_len;
    char *subscripts = (char*)malloc(size*sizeof(char));
    size_t *shape = (size_t*)malloc(size*sizeof(size_t));

    for(size_t i = 0; i < sps->size; i++){
        subscripts[i] = sps->subscripts[i];
        shape[i] = sps->shape[i];

    }size = sps->size;

    // join grade subscripts and subscripts checking if shapes are consistent
    for(size_t j = 0; j < grade_len; j++){
        int found = 0;
        for(size_t i = 0; i < sps->size; i++){
            if(subscripts[i] == grade_subscript[j]){
                if(sps->shape[i] != max_grade+1){ // if the same subscript -> shape must be equal to max_grade
                    return new_tmvs;
                }
                found = 1;
            }
        }
        if(!found){
            subscripts[size] = grade_subscript[j];
            shape[size] = max_grade+1;
            size++;
        }
    }

    // update all different subscripts list
    free(sps->subscripts);
    free(sps->shape);
    sps->size = size;
    sps->subscripts = subscripts;
    sps->shape = shape;

    size_t out_sub_size = strlen(output_subscripts);
    size_t *output_shape = (size_t*)malloc(out_sub_size*sizeof(size_t));
    // check if output has some of the subscripts of the input
    // also compute the shape of the output
    for(size_t i = 0; i < out_sub_size; i++){
        int found = 0;
        for(size_t j = 0; j < size; j++){
            if(subscripts[j] == output_subscripts[i]){
                found = 1;
                output_shape[i] = shape[j];
                break;
            }
        }
        if(!found)
            return new_tmvs;
    }

    size_t output_size = 1;
    for(size_t i = 0; i < out_sub_size; i++)
        output_size *= output_shape[i];

    void *out_data = (void*)malloc(output_size*tmvs.type_size);
    new_tmvs.data = (void*)malloc((tmvs.size+1)*sizeof(void*));
    new_tmvs.shapes = (size_t**)malloc((tmvs.size+1)*sizeof(size_t*));
    new_tmvs.shape_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t));
    new_tmvs.data_size = (size_t*)malloc((tmvs.size+1)*sizeof(size_t));
    for (size_t i = 0; i < tmvs.size; i++) {
        new_tmvs.data[i] = tmvs.data[i];
        new_tmvs.shapes[i] = tmvs.shapes[i];
        new_tmvs.shape_size[i] = tmvs.shape_size[i];
        new_tmvs.data_size[i] = tmvs.data_size[i];
    }

    out->data = new_tmvs.data[tmvs.size] = out_data;
    out->shapes = new_tmvs.shapes[tmvs.size] = output_shape;
    out->shape_size = new_tmvs.shape_size[tmvs.size] = out_sub_size;
    out->data_size = new_tmvs.data_size[tmvs.size] = output_size;
    new_tmvs.size = tmvs.size+1;
    new_tmvs.type_size = tmvs.type_size;

    return new_tmvs;
}


char *get_grade_subscripts(char **subscripts, size_t size){
    size_t diff_size = 0;
    char *diff_subscripts;
    if(!size) return NULL;
    for(size_t i = 0; i < size; i++)
        diff_size += strlen(subscripts[i]);
    diff_subscripts = (char*)malloc((diff_size+1)*sizeof(char));
    diff_size = 0;

    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < strlen(subscripts[i]); j++){
            int found = 0;
            for(size_t k = 0; k < diff_size; k++){
                if(subscripts[i][j] == diff_subscripts[k]){
                    found = 1;
                    break;
                }
            }
            if(!found){
                if(!is_constant_grade(subscripts[i][j])){ // if it is not a number
                    diff_subscripts[diff_size] = subscripts[i][j];
                    diff_size++;
                }
            }
        }
    }
    diff_subscripts[diff_size] = '\0';
    return diff_subscripts;
}

int set_grades(char *sub,
               size_t ***grades,
               size_t **grade_strides,
               char **grade_subscripts,
               size_t *size){

    if(sub == NULL)
        return -1;

    size_t len = strlen(sub);
    if(len != 0){
        grade_strides[*size] = (size_t*)malloc(len*sizeof(size_t));
        grade_subscripts[*size] = (char*)malloc((len+1)*sizeof(char));
        strcpy(grade_subscripts[*size],sub);
        grades[*size] = (size_t**)malloc(len*sizeof(size_t*));
        for(size_t i = 0; i < len; i++){
            grades[*size][i] = (size_t*)malloc(sizeof(size_t));
            grade_strides[*size][i] = 1;
            if(is_constant_grade(sub[i])){
                *(grades[*size][i]) = (size_t)(sub[i] - '0');
                grade_strides[*size][i] = 0;
            }else
                *(grades[*size][i]) = 0; // initialize grades to zero
        }
        return (*size)++;
    }
    return -1;
}


void *get_data_from_symbol(char symbol, char *subscript, char **subscripts, char *symbols, size_t size, tensor_multivectors tmvs, size_t *index){
    for(size_t i = 0; i < size; i++){
        if(!strcmp(subscripts[i],subscript)){
            if(symbol == symbols[i]){
                *index = i;
                return tmvs.data[i];

            }
        }
    }
    return NULL;
}

// it is assumed that **shapes is of size sub.size
tensor_strides compute_tensor_strides(size_t **shapes, subscript_struct sub, subscript_shape sp){
    size_t **strides;
    char *subscripts;
    size_t subscripts_size = 0;
    tensor_strides ts = {NULL,NULL,0,0};

    subscripts_size = sp.size;
    subscripts = sp.subscripts;
    // strides is of shape [n_tensors,n_subscripts]
    strides = (size_t**)malloc(sub.size_*sizeof(size_t*)); // allocate for n_tensors
    for(size_t i = 0; i < sub.size_; i++){
        strides[i] = (size_t*)malloc(subscripts_size*sizeof(size_t)); // allocate for n_subscripts
        for(size_t j = 0; j < subscripts_size; j++){
            strides[i][j] = 0;
        }
    }

    for(size_t k = 0; k < sub.size_; k++){
        char *subs = sub.subscripts[k];
        size_t size = sub.size[k];
        size_t stride = 1;

        for(size_t i = 0; i < size; i++){ // loop over each subscripts
            // look for subscripts in the subscript list
            for(size_t j = 0; j < subscripts_size; j++){
                if(subs[i] == subscripts[j]){
                    // if subscript found increment stride for that subscripts
                    strides[k][j] += stride;
                    break;
                }
            }
            stride *= shapes[k][i];
        }
    }

    // copy shape to n_strides
    ts.n_strides = (size_t*)malloc(sp.size*sizeof(size_t));
    for(size_t i = 0; i < sp.size; i++)
        ts.n_strides[i] = sp.shape[i];


    ts.n_tensors = sub.size_;
    ts.n_subscripts = subscripts_size;
    ts.strides = strides;
    return ts;
}

void free_tensor_strides(tensor_strides ts){
    for(size_t i = 0; i < ts.n_tensors; i++)
        free(ts.strides[i]);
    free(ts.strides);
    free(ts.n_strides);
}

void free_grades_strides(grades_strides gs){
    if(gs.strides != NULL)
        for(size_t i = 0; i < gs.n_grades; i++){
            if(gs.strides[i] != NULL)
                free(gs.strides[i]);
        }free(gs.strides);
}

subscript_shape get_all_subscripts(subscript_struct *sub, size_t **shapes, size_t *shape_size, size_t size){
    size_t subscripts_size = 0;
    char *subscripts;
    subscript_shape sp = {0,NULL,NULL};
    // check if each tensor subscripts are consistent with each corresponding tensor
    // if the number of tensors is different from the subscripts
    if(sub->size_-1 != size)
        return sp;

    sub->size = (size_t*)malloc(sub->size_*sizeof(size_t));

    // number of subscripts for each tensor must be equal to the number of dimensions
    for(size_t j = 0; j < sub->size_; j++){
        sub->size[j] = strlen(sub->subscripts[j]);
        if(j < size && sub->size[j] != shape_size[j]) // don't compare output
            return sp;
    }

    // initialize list of different subscripts
    for(size_t j = 0; j < size; j++)
        subscripts_size += sub->size[j];
    subscripts = (char*)malloc(subscripts_size*sizeof(char));
    size_t *s_shape = (size_t*)malloc(subscripts_size*sizeof(size_t)); // shape for each subscripts
    subscripts_size = 0;

    // determine all the different subscripts;
    for(size_t i = 0; i < size; i++){ // loop over each tensor
        for(size_t j = 0; j < sub->size[i]; j++){ // loop over each symbol
            int found = 0;
            for(size_t k = 0; k < subscripts_size; k++){ // loop over found subscripts
                if(sub->subscripts[i][j] == subscripts[k]){ // found symbol
                    found = 1;
                    if(s_shape[k] != shapes[i][j]) // check if shape is consistent
                        return sp;
                }
            }
            if(!found){ // if symbol not in the list
                subscripts[subscripts_size] = sub->subscripts[i][j]; // add subscript to the list
                s_shape[subscripts_size] = shapes[i][j]; // add subscripts shape to the list
                subscripts_size++;
            }
        }
    }
    sp.size = subscripts_size;
    sp.subscripts = subscripts;
    sp.shape = s_shape;
    return sp;
}

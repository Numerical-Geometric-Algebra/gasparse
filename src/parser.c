#include "parser.h"
#define OPERATORS  "|.^*&"

// a,b,c,d,e,f,g,h
int is_symbol(char c){
    return c <= 'h' && c >= 'a';
}

// i,j,k,l,m,n,o,p,q,r,s,t
int is_subscript(char c){
    return c <= 't' && c >= 'i';
}

int is_constant_grade(char c){
    return c <= '9' && c >= '0';
}


int is_grade_subscripts(char c){
    if(!is_subscript(c))
        return c <= '9' && c >= '0';
    else
        return 1;
}

int is_operator(char c){
    char *operators = OPERATORS;
    for (size_t i = 0; i < strlen(operators); i++){
        if(c == operators[i])
            return 1;
    }
    return 0;
}

expression_struct *parse_expression(char *expression, size_t size, char **output_subscripts){
    expression_struct *es = (expression_struct*)malloc(sizeof(expression_struct));
    init_expression_struct(es);
    char *output_expression = NULL, *input_expression = NULL;
    if(separate_expression(expression,size,&output_expression,&input_expression) == -1)
        return NULL;

    size_t len = strlen(output_expression);
    *output_subscripts = get_subscripts(output_expression,0,len);
    free(output_expression);
    if(*output_subscripts == NULL)
        return NULL;

    recursive_parser(input_expression,strlen(input_expression),es);
    free(input_expression);
    return es;
}

void free_expression_struct_recursive(expression_struct *es){
    if(es->left_sub.grades != NULL)
        free(es->left_sub.grades);
    if(es->left_sub.subscripts != NULL)
        free(es->left_sub.subscripts);

    if(es->right_sub.grades != NULL)
        free(es->right_sub.grades);
    if(es->right_sub.subscripts != NULL)
        free(es->right_sub.subscripts);

    if(es->left != NULL)
        free_expression_struct_recursive(es->left);
    if(es->right != NULL)
        free_expression_struct_recursive(es->right);
    if(es->grades != NULL)
        free(es->grades);
    free(es);
}

int separate_expression(char *expression, size_t size, char **output, char **input){
    int index = -1;
    for(size_t i = 0; i < size; i++){
        if(expression[i] == '='){
            index = i;
            break;
        }
    }

    if(index == -1)
        return -1;

    size_t out_len = index;
    size_t in_len = size - index - 1;

    *output = (char*)malloc((out_len + 1)*sizeof(char));
    *input = (char*)malloc((in_len + 1)*sizeof(char));

    strncpy(*output,expression,out_len);
    (*output)[out_len] = '\0';

    strcpy(*input,expression + out_len + 1);
    return 1;
}

void init_subexpression(sub_expression *sub){
    sub->grades = NULL;
    sub->subscripts = NULL;
    sub->symbol = '\0';
}

void init_expression_struct(expression_struct *es){
    es->left = NULL;
    es->right = NULL;
    es->grades = NULL;
    es->up = NULL;
    es->visited = 0;
    es->operator = '\0';
    init_subexpression(&(es->left_sub));
    init_subexpression(&(es->right_sub));
}


int recursive_parser(char *expression, size_t size, expression_struct *es){
    sub_expression left_sub;
    sub_expression right_sub;
    expression_struct es_sub;

    int beg = 0, end;
    char operator;
    int flag = 0;
    int beg_right;
    es->operator = '\0';

    init_expression_struct(&es_sub);
    end = parse_expression_struct(expression,size,&left_sub,&right_sub,&operator);

    if(end == -1) // No subexpression found
        flag = 1;
    else if(right_sub.symbol == '\0'){ // found left subexpression
        expression += end, size -= end;
        es->left_sub = left_sub;
    }
    else{
        if(size-end > 0){ // there is still expression to be parsed
            expression_struct *left_es = (expression_struct*)malloc(sizeof(expression_struct));
            expression_struct *right_es = (expression_struct*)malloc(sizeof(expression_struct));
            init_expression_struct(left_es);
            init_expression_struct(right_es);

            left_es->left_sub = left_sub;
            left_es->right_sub = right_sub;
            if(is_operator(expression[end])){
                es->operator = expression[end];
                end++;
            } else es->operator = '\0';

            end = recursive_parser(expression + end, size - end, right_es);
            if(end == -1) return -1;

            es->left = left_es;
            es->right = right_es;

            es->right->up = es;
            es->left->up = es;
        }else{
            es->left_sub = left_sub;
            es->right_sub = right_sub;
        }
        return 1;
    }

    if(size <= 0) // nothing more to parse
        return 1;

    end = find_matching_angle_brackets(expression,size,&beg);
    if(end == -1) return -1; // No matching angle brackets found

    if(end + 1 < (int)size){
        int end_grades = end + 1;
        if(expression[end_grades] == '['){
            es->grades = get_subscripts(expression,end+1,size);
            end_grades += 2 + strlen(es->grades);
        }
        if(flag) // only look for operator if no expression was found
            if(end_grades < (int)size && is_operator(expression[end_grades]))
                es->operator = expression[end_grades];
    }
    if(!flag)
        if(is_operator(*expression))
            es->operator = *expression;

    beg_right = end;
    end = recursive_parser(expression + beg + 1,end - beg - 1,&es_sub);
    if(end == -1) return -1;

    if(flag){
        es->left = (expression_struct*)malloc(sizeof(expression_struct));
        *(es->left) = es_sub;
        es->left->up = es;
        init_expression_struct(&es_sub); // reset

        if(es->operator != '\0')
            expression++, size--;
        if(size - beg_right - 1 > 0){ // no more expression to parse
            end = recursive_parser(expression + beg_right + 1, size - beg_right - 1, &es_sub);
            if(end == -1) return -1;
            es->right = (expression_struct*)malloc(sizeof(expression_struct));
            *(es->right) = es_sub;
            es->right->up = es;
        }
    }else{
         es->right = (expression_struct*)malloc(sizeof(expression_struct));
        *(es->right) = es_sub;
        es->right->up = es;
    }

    return 1;
}

int parse_expression_struct(char *expression, size_t size, sub_expression *left_sub,  sub_expression *right_sub, char *operator){
    sub_expression m;
    int end = -1;
    end = sub_expression_parser(expression,size,&m);
    if(end != -1){ // if first argument is a subexpression
        *left_sub = m;
        expression += end; size -= end;
        if(is_operator(*expression)){ *operator = *expression; expression++; size--; }
        else *operator = '\0';
        int end_right = sub_expression_parser(expression,size,&m);
        if(end_right != -1){ // if second argument is a subexpression
            *right_sub = m;
            end += end_right;
        }else{
            right_sub->symbol = '\0'; // no right subexpression found
        }
    }else{
        left_sub->symbol = '\0'; // no left subexpression found
    }
    return end;
}



int sub_expression_parser(char *expression, size_t size, sub_expression *m){
    int beg = 0;
    int end = 0;
    char *subscripts = NULL , *grades = NULL;
    char symbol = -1;
    if(size){
        if(is_symbol(*expression)){
            beg = 1;
            symbol = *expression;
            if(expression[beg] == '['){ // if symbol has subscripts
                subscripts = get_subscripts(expression,beg,size);
                if(subscripts == NULL)
                    return -1;
                end = beg + strlen(subscripts) + 1;
            }
        }else{
            end = find_matching_angle_brackets(expression,size,&beg);
            if(end == -1) return -1;
            if(beg != 0)  return -1;
            if(is_symbol(expression[1])){
                beg = 2;
                symbol = expression[1];
            }else return -1;

            // only want to parse one symbol
            for(int i = beg; i < end; i++){
                if(is_symbol(expression[i]))
                    return -1;
            }
            if(end + 1 < (int)size && expression[end+1] == '[')
                grades = get_subscripts(expression,end,size);
            subscripts = get_subscripts(expression,beg,end);
            if(grades != NULL)
                end += (int)strlen(grades)+2;
        }
    }else return -1;

    m->symbol = symbol;
    m->grades = grades;
    m->subscripts = subscripts;

    return end + 1;
}


char *get_subscripts(char *expression, int beg, int end){
    int sub_count = -1;
    // find square bracket and matching
    for(int i = beg; i < end; i++){
        if(sub_count != -1 && expression[i] == '[')
            return NULL;
        if(expression[i] == '[')
            sub_count = i;
        if(expression[i] == ']'){
            sub_count = i-sub_count-1;
            break;
        }
        if(sub_count != -1){
            if(is_symbol(expression[i]))
                return NULL;
        }
    }

    if(sub_count == -1) return NULL;

    char *subscripts = (char*)malloc((sub_count+1)*sizeof(char));
    for(int i = beg; i < end; i++){
        if(expression[i] == '['){
            int j = 0;
            while(++i < end && expression[i] != ']'){
                subscripts[j++] = expression[i];
            }break;
        }
    }
    subscripts[sub_count] = '\0';
    return subscripts;
}


int find_matching_angle_brackets(char *expression, size_t size, int *beg){
    size_t i;
    *beg = -1;
    for(i = 0; i < size; i++)
        if(expression[i] == '<')
            break;
    if(i == size)
        return -1;
    *beg = i;
    expression = expression + i;
    size_t count = 1;
    for(size_t j = 1; j < size; j++){ // find matching angle bracket
        if(expression[j] == '<')
            count++;
        if(expression[j] == '>'){
            count--;
            if(!count){
                return j + i;
            }
        }
    }
    return -1;
}

#include "cayley.h"


void free_map(map m){
    for(size_t i = 0; i < m.size; i++){
        free(m.sign[i]);
        free(m.bitmap[i]);
    }
    free(m.sign);
    free(m.bitmap);
}

// Determines the number of one bits in an integer
// Which is equivalent to determining the grade of a bitset
unsigned int grade(unsigned int v){
    v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}

void free_grade_map(grade_map m){
    free(m.grade);
    free(m.position);
    free(m.grade_size);
}

// Determines the position of the bitmap on the corresponding grade
grade_map bitmap_grade_map(size_t size){
    grade_map m;
    unsigned int max_grade = grade(size-1);
    unsigned int *g_pos = (unsigned int*)malloc((max_grade + 1)*sizeof(unsigned int));
    m.grade = (unsigned int*)malloc(size*sizeof(unsigned int));
    m.position = (unsigned int*)malloc(size*sizeof(unsigned int));
    for(size_t i = 0; i <= max_grade; i++)
        g_pos[i] = 0;

    for(size_t i = 0; i < size; i++){
        m.grade[i] = grade(i);
        m.position[i] = g_pos[m.grade[i]]++; // assign value then increment
    }
    m.size = size;
    m.grade_size = g_pos;
    m.max_grade = max_grade;
    return m;
}

map init_map(size_t n){
    map m;
    int **sign = (int**)malloc(n*sizeof(int*));
    unsigned int **bitmap = (unsigned int**)malloc(n*sizeof(unsigned int*));

    for(size_t i = 0; i < n; i++){
        sign[i] = (int*)malloc(n*sizeof(int));
        bitmap[i] = (unsigned int*)malloc(n*sizeof(unsigned int));
    }
    m.bitmap = bitmap;
    m.sign = sign;
    m.size = n;
    return m;
}

map cayley_table(size_t p, size_t q, size_t r){
    map m;
    size_t n = 1 << (p+q+r);
    // allocate memory for the whole table
    m = init_map(n);

    m.sign[0][0] = 1;// initialize algebra of scalars
    for(size_t i = 0; i < p; i++)//loop through the identity-square basis vectors
        sub_algebra(i,m.sign,1);
    for(size_t i = p; i < p+q; i++)//loop through the negative-square basis vectors
        sub_algebra(i,m.sign,-1);
    for(size_t i = p+q; i < p+q+r; i++)//loop through the zero-square basis vectors
        sub_algebra(i,m.sign,0);

    // determine each basis blade and its grade
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < n; j++){
            unsigned int bitmap_ij = i ^ j;
            m.bitmap[i][j] = bitmap_ij;
        }
    }

    return m;
}

void sub_algebra(unsigned int k, int **s, int metric){
    size_t m = 1 << k; // same as 2^k
    size_t n = m << 1; // same as 2^(k+1)
    int sign;
    // This could be improved by checking if the element in the array is zero
    for(size_t i = m; i < n; i++){// loop through the new elements
        for(size_t j = 0; j < m; j++){// loop through old elements
            // j is indepedent of the new basis vector
            sign = ((grade(j) & 1) == 0) ? 1 : -1;// return minus one if grade is odd
            s[i][j] = sign*s[i-m][j];
            s[j][i] = s[j][i-m];
        }
        if(metric != 0 ){
            for(size_t j = m; j < n; j++){// loop through new elements
                // These elements have the new basis vector in common
                sign = metric;
                // remove the new basis vector then determine sign
                sign *= ((grade(j-m) & 1) == 0) ? 1 : -1;
                sign *= s[i-m][j-m];// remove the new vector part
                s[i][j] = sign;
            }
        }else{//if null metric -> set all elements to zero
            for(size_t j = m; j < n; j++)
                s[i][j] = 0;
        }
    }
}

unsigned int* get_grade_bool(unsigned int *grades, size_t size, size_t max_grade){
    unsigned int *g = (unsigned int*)malloc(max_grade*sizeof(unsigned int));
    if(size == 0){ // if size is 0 project to all grades
        for(size_t i = 0; i < max_grade; i++)
            g[i] = 1;
    }else{
        for(size_t i = 0; i < max_grade; i++)
            g[i] = 0;

        for(size_t i = 0; i < size; i++)
            g[grades[i]] = 1;
    }
    return g;
}


map invert_map(map m){
    map m_inv = init_map(m.size);
    for(size_t i = 0; i < m.size; i++){
        for(size_t j = 0; j < m.size; j++){
            m_inv.bitmap[i][m.bitmap[i][j]] = j;
            m_inv.sign[i][m.bitmap[i][j]] = m.sign[i][j];
        }
    }
    return m_inv;

}

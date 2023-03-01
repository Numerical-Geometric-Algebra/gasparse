#include "cayley.h"


void free_map(map m){
    for(size_t i = 0; i < m.size; i++){
        free(m.sign[i]);
        free(m.bitmap[i]);
        free(m.grade[i]);
    }
    free(m.sign);
    free(m.bitmap);
    free(m.grade);
}


void print_map(map m){
     for(size_t i = 0; i < m.size; i++ ){
        for(size_t j = 0; j < m.size; j++ ){
            printf(" (%zu,%zu) | ", i, j);
            printf("(%d,%d,%d)\n", m.sign[i][j], m.grade[i][j], m.bitmap[i][j]);
        }
        printf("\n");
    }
}

// Determines the number of one bits in an integer
// Which is equivalent to determining the grade of a bitset
unsigned int grade(unsigned int v){
    v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}

map cayley_table(size_t p, size_t q, size_t r){
    map algebra_map;
    size_t n = 1 << (p+q+r);
    int **sign = (int**)malloc(n*sizeof(int*));
    unsigned int **bitmap = (unsigned int**)malloc(n*sizeof(unsigned int*));
    int **g = (int**)malloc(n*sizeof(int*));
    // allocate memory for the whole table
    for(size_t i = 0; i < n; i++){
        sign[i] = (int*)malloc(n*sizeof(int));
        bitmap[i] = (unsigned int*)malloc(n*sizeof(unsigned int));
        g[i] = (int*)malloc(n*sizeof(int));
    }

    sign[0][0] = 1;// initialize algebra of scalars
    for(size_t i = 0; i < p; i++)//loop through the identity-square basis vectors
        sub_algebra(i,sign,1);
    for(size_t i = p; i < p+q; i++)//loop through the negative-square basis vectors
        sub_algebra(i,sign,-1);
    for(size_t i = p+q; i < p+q+r; i++)//loop through the zero-square basis vectors
        sub_algebra(i,sign,0);

    // determine each basis blade and its grade
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < n; j++){
            unsigned int bitmap_ij = i ^ j;
            bitmap[i][j] = bitmap_ij;
            g[i][j] = grade(bitmap_ij);
        }
    }
    algebra_map.sign = sign;
    algebra_map.bitmap = bitmap;
    algebra_map.grade = g;
    algebra_map.size = n;
    return algebra_map;
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

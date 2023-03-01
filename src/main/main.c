#include "main.h"
#include <time.h>
#define MAX_METRIC_SIZE 32

int main(){

    clock_t start, end;
    double cpu_time_used1;
    double cpu_time_used2;
    size_t size;
    size_t p,q,r;

    // For metric with more than 14 basis vector the cpu runs out of memory
    for(size_t n = 1; n <= 14; n++){
        p = n; q = 0; r = 0;
        start = clock();
        map m = cayley_table_map(p,q,r);
        end = clock();
        cpu_time_used1 = ((double)(end-start))/CLOCKS_PER_SEC;
        free_map_map(m);

        start = clock();
        algebra_map **map = cayley_table(p,q,r,&size);
        end = clock();
        cpu_time_used2 = ((double)(end-start))/CLOCKS_PER_SEC;
        printf("(p,q,r)=(%zu,%zu,%zu) ",p,q,r);
        printf("(old,subalgebra)=(%f,%f)\n",cpu_time_used2,cpu_time_used1);
        free_map(map,size);
    }

    /* p = 2; q = 2; r = 2; */
    /* map m = cayley_table_map(p,q,r); */
    /* print_map_map(m); */
    /* free_map_map(m); */

    /* print_map(map,mv_size); */
    /* algebra_map **map = cayley_table(2,1,1,&mv_size); */


}

void free_map_map(map m){
    for(size_t i = 0; i < m.size; i++){
        free(m.sign[i]);
        free(m.bitmap[i]);
        free(m.grade[i]);
    }
    free(m.sign);
    free(m.bitmap);
    free(m.grade);
}


void print_map_map(map m){
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
unsigned int compute_grade(unsigned int v){
    v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}


// Determines the sign when reordering basis blades
int sign_reorder(mixed_bitmap a,mixed_bitmap b){
    unsigned int x[3] = {a.p,a.q,a.r};
    unsigned int y[3] = {b.p,b.q,b.r};

    unsigned int sum[3] = {0,0,0};
    unsigned int n_swaps = 0;
    for(int i=0; i<3; i++){
        x[i] = x[i] >> 1;
        while(x[i] != 0){
            sum[i] += compute_grade(x[i] & y[i]);
            x[i] = x[i] >> 1;
        }
    }

    // If it is a scalar commutes with everything
    // Only add the grade if it's not a scalar
    n_swaps = sum[0] + sum[1] + sum[2];
    if(a.q != 0)
        n_swaps += compute_grade(b.p);

    if(a.r != 0)
        n_swaps += compute_grade(b.q) + compute_grade(b.p);

    // returns 1 if n_swaps is even -1 otherwise
    return  ((n_swaps & 1) == 0) ? 1 : -1;

}


algebra_map geo_prod(mixed_bitmap a, mixed_bitmap b){
    mixed_bitmap y;
    unsigned int sign = 1;

    y.p = a.p ^ b.p;
    y.q = a.q ^ b.q;
    y.r = a.r ^ b.r;

    sign *= ((a.r & b.r) == 0) ? 1 : 0;
    if(sign){
        sign *= ((compute_grade(a.q & b.q) & 1) == 0) ? 1 : -1;
        sign *= sign_reorder(a,b);
    }
    unsigned int grade = compute_grade(y.p) + compute_grade(y.q) + compute_grade(y.r);
    algebra_map map = {sign,y,grade};

    return map;

}

algebra_map **cayley_table(size_t p, size_t q, size_t r, size_t *mv_size){
    size_t n = *mv_size = (size_t)pow((double)2,(double)(p+q+r));

    size_t p_max = (size_t)pow((double)2,(double)(p));
    size_t q_max = (size_t)pow((double)2,(double)(q));
    size_t r_max = (size_t)pow((double)2,(double)(r));

    /* mixed_bitmap metric = {p,q,r}; */

    algebra_map **map = (algebra_map**)malloc(n*sizeof(algebra_map*));
    mixed_bitmap a,b;
    for(a.p = 0; a.p < p_max; a.p++){
        for(a.q = 0; a.q < q_max; a.q++){
            for(a.r = 0; a.r < r_max; a.r++){
                size_t i = (a.p*q_max + a.q)*r_max + a.r;
                map[i] = malloc(n*sizeof(algebra_map));
                for(b.p = 0; b.p < p_max; b.p++){
                    for(b.q = 0; b.q < q_max; b.q++){
                        for(b.r = 0; b.r < r_max; b.r++){
                            size_t j = (b.p*q_max + b.q)*r_max + b.r;
                            map[i][j] = geo_prod(a,b);
                            /* printf("{(%d,%d,%d),(%d,%d,%d)} | ",a.p,a.q,a.r,b.p,b.q,b.r); */
                            /* printf("(%d,%d) | ", map[i][j].sign, map[i][j].grade); */
                            /* printf("(%d,%d,%d)\n",map[i][j].index.p,map[i][j].index.q,map[i][j].index.r); */
                        }
                    }
                }
            }
        }
    }
    return map;
}

map cayley_table_map(size_t p, size_t q, size_t r){
    map algebra_map;
    size_t n = 1 << (p+q+r);
    int **sign = (int**)malloc(n*sizeof(int*));
    int **bitmap = (int**)malloc(n*sizeof(int*));
    int **grade = (int**)malloc(n*sizeof(int*));
    // allocate memory for the whole table
    for(size_t i = 0; i < n; i++){
        sign[i] = (int*)malloc(n*sizeof(int));
        bitmap[i] = (int*)malloc(n*sizeof(int));
        grade[i] = (int*)malloc(n*sizeof(int));
    }

    sign[0][0] = 1;// initialize algebra of scalars
    for(size_t i = 0; i < p; i++)//loop through the identity-square basis vectors
        sub_algebra_cayley(i,sign,1);
    for(size_t i = p; i < p+q; i++)//loop through the negative-square basis vectors
        sub_algebra_cayley(i,sign,-1);
    for(size_t i = p+q; i < p+q+r; i++)//loop through the zero-square basis vectors
        sub_algebra_cayley(i,sign,0);

    // determine each basis blade and its grade
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < n; j++){
            int bitmap_ij = i ^ j;
            bitmap[i][j] = bitmap_ij;
            grade[i][j] = compute_grade(bitmap_ij);
        }
    }
    algebra_map.sign = sign;
    algebra_map.bitmap = bitmap;
    algebra_map.grade = grade;
    algebra_map.size = n;
    return algebra_map;
}

void sub_algebra_cayley(unsigned int k, int **s, int metric){
    size_t m = 1 << k; // same as 2^k
    size_t n = m << 1; // same as 2^(k+1)
    int sign;
    // This could be improved by checking if the element in the array is zero
    for(size_t i = m; i < n; i++){// loop through the new elements
        for(size_t j = 0; j < m; j++){// loop through old elements
            // j is indepedent of the new basis vector
            sign = ((compute_grade(j) & 1) == 0) ? 1 : -1;// return minus one if grade is odd
            s[i][j] = sign*s[i-m][j];
            s[j][i] = s[j][i-m];
        }
        if(metric != 0 ){
            for(size_t j = m; j < n; j++){// loop through new elements
                // These elements have the new basis vector in common
                sign = metric;
                // remove the new basis vector then determine sign
                sign *= ((compute_grade(j-m) & 1) == 0) ? 1 : -1;
                sign *= s[i-m][j-m];// remove the new vector part
                s[i][j] = sign;
            }
        }else{//if null metric -> set all elements to zero
            for(size_t j = m; j < n; j++)
                s[i][j] = 0;
        }
    }
}

void free_map(algebra_map **map, size_t mv_size){
    for(size_t i = 0; i < mv_size; i++ )
        free(map[i]);
    free(map);
}

void print_map(algebra_map **map, size_t mv_size){
     for(size_t i = 0; i < mv_size; i++ ){
        for(size_t j = 0; j < mv_size; j++ ){
            printf("(%zu,%zu) | ", i, j);
            printf("(%d,%d) | ", map[i][j].sign, map[i][j].grade);
            printf("(%d,%d,%d)\n",map[i][j].index.p,map[i][j].index.q,map[i][j].index.r);
        }
        printf("\n");
    }
}

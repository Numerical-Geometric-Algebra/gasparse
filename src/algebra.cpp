// #include <iostream>
// using namespace std;
#include "algebra.hpp"
#include <Python.h>
// #include <iostream> // header in standard library


#define GRADE(value) (__builtin_popcountll(value))

using namespace CliffordAlgebra;


// class Algebra {
    // public:

    
    
// void LinkedArrayList<ItemType>::insert (int index, const ItemType& item)
    template<typename T,typename M>
    void Algebra<T,M>::clifford_sub_algebra(Py_ssize_t k, char **s, int metric) {
        Py_ssize_t m = 1 << k; // same as 2^k
        Py_ssize_t n = m << 1; // same as 2^(k+1)
        int sign;
        // This could be improved by checking if the element in the array is zero
        for (Py_ssize_t i = m; i < n; i++) {   // loop through the new elements
            for (Py_ssize_t j = 0; j < m; j++) { // loop through old elements
                // j is indepedent of the new basis vector
                sign = ((GRADE(j) & 1) == 0) ? 1 : -1; // return minus one if grade is odd
                s[i][j] = sign * s[i - m][j];
                s[j][i] = s[j][i - m];
            }
            if (metric != 0) {
                for (Py_ssize_t j = m; j < n; j++) { // loop through new elements
                    // These elements have the new basis vector in common
                    sign = metric;
                    // remove the new basis vector then determine sign
                    sign *= ((GRADE(j - m) & 1) == 0) ? 1 : -1;
                    sign *= s[i - m][j - m]; // remove the new vector part
                    s[i][j] = sign;
                }
            } else { // if null metric -> set all elements to zero
                for (Py_ssize_t j = m; j < n; j++)
                    s[i][j] = 0;
            }
        }
    }

    template<typename T,typename M>
    void Algebra<T,M>::map_dealloc(Algebra::CliffordMap *self) {
        if (self->sign) {
            for (Py_ssize_t i = 0; i < self->size; i++)
                PyMem_RawFree(self->sign[i]), self->sign[i] = NULL;
            PyMem_RawFree(self->sign);
            self->sign = NULL;
        }

        if (self->bitmap) {
            for (Py_ssize_t i = 0; i < self->size; i++)
                PyMem_RawFree(self->bitmap[i]), self->bitmap[i] = NULL;
            PyMem_RawFree(self->bitmap);
            self->bitmap = NULL;
        }
    }

    template<typename T, typename W>
    int Multivector<T, W>::alloc_mvarray_data(){
        // Alloc memory for the data
        if(this->strides != NULL){
            this->MultivectorData = (void*)PyMem_RawMalloc(this->strides[0]*sizeof(W));
            if(!this->data) return 0;

            for(Py_ssize_t i = 0; i < this->strides[0]; i++){
                if(!Multivector::init(i,NULL,NULL,0))
                    return 0;
            }
        }
    return 1;
    }


    template <typename T>
    inline int MultivectorSparse<T>::init(Py_ssize_t i, int *bitmap, T *value, Py_ssize_t size){
        this->MultivectorData->values = (T*)PyMem_RawMalloc(size*sizeof(T));
        this->MultivectorData->bitmap = (int*)PyMem_RawMalloc(size*sizeof(int));
        this->MultivectorData->size = size;
        for(Py_ssize_t i = 0; i < size; i++){
            this->MultivectorData->value[i] = value[i];
            this->MultivectorData->bitmap[i] = bitmap[i];
        }
        return 1;
    }

// };
#include "common.h"


Py_ssize_t parse_list_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades){
    Py_ssize_t size = -1;
    if(PyLong_Check(grades_obj)){ // check if object is an integer
        int grade = (int)PyLong_AsLong(grades_obj);
        if(grade > MAX_GRADE(ga) || grade < 0) // not a valid grade
            return -1; // raise error
        *(grades) = (int*)PyMem_RawMalloc(sizeof(int));
        **grades = grade;
        size = 1;
    }else if(PyList_Check(grades_obj)){ // check if object is a list type
        size = PyList_Size(grades_obj);
        if(!size) return -1;
        *grades = (int*)PyMem_RawMalloc(size*sizeof(int));
        for(Py_ssize_t i = 0; i < size; i++){
            PyObject *grade_obj = PyList_GetItem(grades_obj,i);
            if(!PyLong_Check(grade_obj))
                return -1; // raise error
            (*grades)[i] = (int)PyLong_AsLong(grade_obj);
            if((*grades)[i] > MAX_GRADE(ga)){
                PyMem_RawFree(grades);
                return -1; // raise error
            }
        }
    }
    return size;
}

Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades){
    Py_ssize_t *g = (Py_ssize_t*)PyMem_RawMalloc(n_grades*sizeof(Py_ssize_t));
    if(!g){
        PyErr_SetString(PyExc_MemoryError,"Error allocating memory");
        return NULL;
    }
    if(size == 0){ // if size is 0 project to all grades
        for(Py_ssize_t i = 0; i < n_grades; i++)
            g[i] = 1;
    }else{
        for(Py_ssize_t i = 0; i < n_grades; i++)
            g[i] = 0;
        for(Py_ssize_t i = 0; i < size; i++)
            g[grades[i]] = 1;
    }
    return g;
}
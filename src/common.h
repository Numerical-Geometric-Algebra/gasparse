#ifndef COMMON_H_
#define COMMON_H_
#include "types.h"

#define GRADE(value) (__builtin_popcountll(value))
typedef double ga_float;

Py_ssize_t parse_list_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades);


#endif // COMMON_H_

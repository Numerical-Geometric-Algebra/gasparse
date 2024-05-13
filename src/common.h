#ifndef COMMON_H_
#define COMMON_H_
#include "types.h"

#define GRADE(value) (__builtin_popcountll(value))


Py_ssize_t parse_list_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
Py_ssize_t parse_arguments_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
Py_ssize_t* get_grade_bool(int *grades, Py_ssize_t size, Py_ssize_t n_grades);
Py_ssize_t parse_list_tuple_as_grades(PyAlgebraObject *ga, PyObject *grades_obj, int **grades);
ProductType string_to_product_type(char *type_str);


#endif // COMMON_H_

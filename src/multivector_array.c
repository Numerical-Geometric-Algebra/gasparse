#include "multivector_array.h"
#include "common.h"
#include "gasparse.c"
#include "gasparse.h"
#include "listobject.h"
#include "pyport.h"
#include "pytypedefs.h"




static Py_ssize_t get_ndims(PyObject *data){
    Py_ssize_t ndims = 0;
    PyObject *sublist = data;
    while(PyList_Check(sublist)){
        if((sublist = PyList_GetItem(sublist,0)) == NULL) break; // Break on empty list
        ndims++;
    }
    return ndims;
}

static Py_ssize_t *get_shapes(PyObject *data, Py_ssize_t ndims){
    Py_ssize_t *shapes = (Py_ssize_t*)PyMem_RawMalloc(ndims*sizeof(Py_ssize_t));
    if(!shapes) return NULL;
    PyObject *sublist = data;
    Py_ssize_t i = 0;
    while(PyList_Check(sublist)){
        shapes[i] = PyList_Size(sublist);
        if(!shapes[i]) break;
        sublist = PyList_GetItem(sublist,0);
        i++;
        if(i >= ndims) break;
    }
    return shapes;
}

static int iterate_nested_lists(PyObject *list,
                                ga_float **array, 
                                Py_ssize_t *strides,
                                Py_ssize_t *shape,
                                Py_ssize_t index, 
                                Py_ssize_t dim,
                                Py_ssize_t nbasis){
    
    if(!PyList_Check(list))
        return 1; // Deepest nested element 
    
    PyObject *sublist = list;
    if(PyList_Size(sublist) != shape[dim]) return -1; // Not the right shape
    for(Py_ssize_t i = 0; i < PyList_Size(sublist); i++){
        sublist = PyList_GetItem(sublist, i);
        int flag = iterate_nested_lists(sublist, array, strides, shape, index + i*strides[dim], dim+1, nbasis);
        if(flag == 1){
            // most inner list
            Py_ssize_t size = parse_list_as_values(list, &array[index + i*strides[dim]]);
            if(size <= 0 && size != nbasis) return -1;
        }else if(flag == -1){
            return -1;
        }
    }

    return 0;
}

static Py_ssize_t *get_strides(Py_ssize_t ndims, Py_ssize_t *shapes){
    Py_ssize_t *strides = (Py_ssize_t*)PyMem_RawMalloc(ndims*sizeof(Py_ssize_t));
    if(!strides) return NULL;
    strides[ndims-1] = 1;
    for(Py_ssize_t i = ndims - 2; i >= 1; i--)
        strides[i] = strides[i+1]*shapes[i];
    return strides;
}





static PyObject *algebra_multivector_array(PyAlgebraObject *self, PyObject *args,
																		 PyObject *kwds) {
	static char *kwlist[] = {"values", "basis", "grades","dtype", NULL};
	PyObject *values = NULL, *basis = NULL, *grades = NULL;
	int *bitmaps_int = NULL;
	int *grades_int = NULL;
	ga_float *values_float = NULL;
	ga_float *values_conv = NULL;
	Py_ssize_t size, bsize;
	PyMultivectorObject *multivector;
	char *type_name = NULL;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OOs", kwlist, &values, &basis, &grades, &type_name))
		return NULL;

	if (!values)
		PyErr_SetString(PyExc_ValueError, "Values must be non empty");

	if (grades && basis) {
		PyErr_SetString(PyExc_ValueError, "Can only define multivectors through "
																			"basis blades or grades, not both!");
		return NULL;
	}

	// Defining multivectors through some basis
	if (basis) {
		if (PyList_Check(basis)) {
			PyErr_SetString(PyExc_ValueError,
											"Basis must be a list and of the same size as values");
			return NULL;
		}
		
        bsize = parse_list_as_bitmaps(basis, &bitmaps_int);
        if (bsize <= 0) {
            PyErr_SetString(PyExc_TypeError, "Error parsing basis list as bitmaps");
            return NULL;
        }
		
	} else if (grades) {
		
		Py_ssize_t gsize = parse_list_as_grades(self, grades, &grades_int);
		if (gsize <= 0) {
			PyMem_RawFree(values_float);
			PyErr_SetString(PyExc_ValueError,
											"Error parsing grades, invalid value or empty");
			return NULL;
		}
		Py_ssize_t mv_size =
				parse_list_as_basis_grades(*self, grades_int, &bitmaps_int, gsize);

		if (mv_size != size) {
			PyMem_RawFree(grades_int);
			PyMem_RawFree(bitmaps_int);
			PyMem_RawFree(values_float);
			PyErr_SetString(PyExc_ValueError,
											"Basis grades must be of the same size as values");
			return NULL;
		}
        bsize = mv_size;

	} else {
		// Defining a multivector through the whole algebra basis
		size = self->asize;
		bitmaps_int = (int *)PyMem_RawMalloc(size * sizeof(int));
		for (int i = 0; i < size; i++)
			bitmaps_int[i] = i;
        bsize = size;
	}
	if(!type_name)
		type_name = self->mdefault.type_name;

    size = parse_list_as_values(values, &values_float);
	if (size){ 
        // It's a one dimentional list

        multivector = new_multivector(self,type_name);
        if (!multivector) {
            PyMem_RawFree(values_float);
            PyMem_RawFree(bitmaps_int);
            PyErr_SetString(PyExc_ValueError,
                                                "dtype name invalid, select from...");
            return NULL;
        }

        gainitfunc init = multivector->type->data_funcs->init;
        if (init)
            multivector->data = init(bitmaps_int, values_float, size, self);
        else {
            PyMem_RawFree(values_float);
            PyMem_RawFree(bitmaps_int);
            return NULL; // raise not implemented error
        }

        PyMem_RawFree(values_float);
        PyMem_RawFree(bitmaps_int);
        PyMem_RawFree(grades_int);

        return (PyObject *)multivector;
    }else{
        Py_ssize_t ndims = get_ndims(values);
        if(ndims <= 0){
            PyErr_SetString(PyExc_ValueError,"Error geting dims");
            return NULL;
        }
        Py_ssize_t *shapes = get_shapes(values,ndims);
        if(!shapes){
            PyErr_SetString(PyExc_ValueError,"Error geting shapes");
            return NULL;
        }
        Py_ssize_t *strides = get_strides(ndims,shapes);
        if(!strides){
            PyErr_SetString(PyExc_ValueError,"Error geting strides");
            return NULL;
        }
         
        
        ga_float **values_float = (ga_float**)PyMem_RawMalloc(strides[0]*sizeof(ga_float*));
        if(iterate_nested_lists(values,values_float,strides,shapes,0,ndims,bsize) == -1) {
            PyErr_SetString(PyExc_ValueError,"Error iterating nested lists");
            return NULL;
        }
    }
    Py_RETURN_NONE;
}



/*
static PyObject* read_numpy_array(PyObject* self, PyObject* args) {

    PyObject* numpy_array;
    if (!PyArg_ParseTuple(args, "O", &numpy_array)) {
        return NULL;
    }

    int ndim = PyArray_NDIM(numpy_array);
    npy_intp *dims = PyArray_DIMS(numpy_array);


     if (!PyArray_Check(numpy_array)) {
        switch(PyArray_TYPE((PyTypeObject*)numpy_array)){
            case NPY_INT:
            break;
        }
     }
   

}

static PyObject *get_nested_list_item(PyObject *list, Py_ssize_t *index,  Py_ssize_t ndims){
    PyObject *sublist = list;
    Py_ssize_t i = 0;
    while(PyList_Check(sublist)){
        if(index[i] >= PyList_Size(sublist)) return NULL;
        sublist = PyList_GetItem(sublist,index[i]);
        i++;
    }
    return sublist;
}

static Py_ssize_t *get_indices(Py_ssize_t i, Py_ssize_t ndims, Py_ssize_t *shapes){
    Py_ssize_t *sizes = (Py_ssize_t*)PyMem_RawMalloc(ndims*sizeof(Py_ssize_t));
    sizes[0] = 1;
    for(Py_ssize_t i = 1; i < ndims; i++)
        sizes[i] = sizes[i-1]*shapes[i];
}


static PyMultivectorObject* multivector_array_from_lists(PyObject *list, Py_ssize_t ndims,Py_ssize_t *shapes){
    Py_ssize_t size = 1; // size of the multivector array

    for(Py_ssize_t i = 0; i < ndims; i++)
        size *= shapes[i];
    
    PyMultivectorObject* multivector_array_data = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    Py_ssize_t j = 0;
    PyObject *sublist = list;
    PyObject *rootlist = NULL;
    while(PyList_Check(sublist)){
        sublist = PyList_GetItem(sublist,0);
        if(!sublist) break;
        rootlist = PyList_GetItem(sublist,0);
        if(!PyList_Check(rootlist)){
            rootlist = sublist;
            break;
        }
    }
    while(1){
        for(Py_ssize_t i = 0; i < shapes[j]; i++){

        }
    }

}
*/

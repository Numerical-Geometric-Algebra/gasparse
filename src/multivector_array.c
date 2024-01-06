#include "multivector_array.h"
#include "common.h"
#include "gasparse.h"
#include "listobject.h"
#include "modsupport.h"
#include "pyerrors.h"
#include "pyport.h"
#include "pytypedefs.h"


/* Multivector array functions */
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

static int iterate_nested_lists_1(PyObject *list,
                                ga_float **array, 
                                Py_ssize_t *strides,
                                Py_ssize_t *shape,
                                Py_ssize_t index, 
                                Py_ssize_t dim,
								Py_ssize_t ndims,
                                Py_ssize_t nbasis){
    
    if(PyList_Size(list) != shape[dim]) return -1; // Not the right shape
    if(!PyList_Check(list)) return -1; // Deepest nested element 
    
    if(dim == ndims){
        Py_ssize_t size = parse_list_as_values(list, &array[index]);
        if(size <= 0 && size != nbasis) return -1;
    }else{
        for(Py_ssize_t i = 0; i < (Py_ssize_t)PyList_Size(list); i++){
			PyObject *sublist = PyList_GetItem(list, i);
			int flag = iterate_nested_lists_1(sublist, array, strides, shape, index + i*strides[dim+1], dim+1,ndims,nbasis);
			if(flag == -1) return -1;
		}
    }

    return 0;
}


static int iterate_nested_lists(PyObject *list,
                                ga_float **array, 
                                Py_ssize_t *strides,
                                Py_ssize_t *shape,
                                Py_ssize_t index, 
                                Py_ssize_t dim,
								Py_ssize_t ndims,
                                Py_ssize_t nbasis){
    
    if(!PyList_Check(list))
        return -1; // Deepest nested element 
    
    if(PyList_Size(list) != shape[dim]) return -1; // Not the right shape
	if(dim == ndims-1){ // The innermost dimension
		for(Py_ssize_t i = 0; i < (Py_ssize_t)PyList_Size(list); i++){
        	PyObject *sublist = PyList_GetItem(list, i);
			Py_ssize_t size = parse_list_as_values(sublist, &array[index + i*strides[dim+1]]);
            if(size <= 0 && size != nbasis) return -1;
		}
	}else if(ndims-1 > 0){
		for(Py_ssize_t i = 0; i < (Py_ssize_t)PyList_Size(list); i++){
			PyObject *sublist = PyList_GetItem(list, i);
			int flag = iterate_nested_lists(sublist, array, strides, shape, index + i*strides[dim+1], dim+1,ndims,nbasis);
			if(flag == -1) return -1;
		}
	}else{
		Py_ssize_t size = parse_list_as_values(list,array);
        if(size <= 0 && size != nbasis) return -1;
    }

    return 0;
}

static Py_ssize_t *get_strides(Py_ssize_t ndims, Py_ssize_t *shapes){
    Py_ssize_t *strides = (Py_ssize_t*)PyMem_RawMalloc((ndims+1)*sizeof(Py_ssize_t));
    if(!strides) return NULL;
    strides[ndims] = 1;

    for(Py_ssize_t i = ndims; i >= 1; i--)
        strides[i-1] = strides[i]*shapes[i-1];
    return strides;
}

static void multivector_array_dealloc(PyMvArrayObj *self){
	void **data = self->data;
    gafreefunc free_type = self->type->data_funcs->free;
    if(free_type)
    for(Py_ssize_t i = 0; i < self->strides[0]; i++){
        if(self->data[i]){
            free_type(self->data[i]);
            PyMem_RawFree(self->data[i]);
        }
    }

	Py_XDECREF((PyObject*)self->GA);
	PyMem_RawFree(self->strides);
	PyMem_RawFree(self->shapes);
	PyMem_RawFree(data);
	PyMem_RawFree(self);
}

static PyMvArrayObj *init_multivector_array(PyAlgebraObject *ga, Py_ssize_t ndims, Py_ssize_t *strides,Py_ssize_t *shapes){
	PyMvArrayObj *array_obj = (PyMvArrayObj*)PyMem_RawMalloc(sizeof(PyMvArrayObj));
	array_obj->strides = strides;
	array_obj->shapes = shapes;
	array_obj->data = NULL;
	array_obj->ns = ndims;
	// set type and increase reference count
    array_obj->GA = ga;
    // Set the mixed type operations table
    array_obj->mixed = ga->mixed;
    Py_XINCREF((PyObject*)array_obj->GA);
	Py_SET_TYPE(array_obj, &PyMultivectorArrayType);
	Py_SET_REFCNT((PyObject*)array_obj,1);

	return array_obj;
}

static int multivector_array_iter_next(PyMvArrayIter *iter, Py_ssize_t dim){
    if(dim < 0) return 1; // Ignore iteration
    if(dim >= iter->ns) return 0;
    
    iter->data += iter->strides[dim+1];
    iter->index[dim]++;
    if(iter->index[dim] >= iter->shapes[dim]){
        iter->data -= iter->shapes[dim]*iter->strides[dim+1];
        iter->index[dim] = 0;
        return 0; // last iteration
    }
    return 1;
}

static int mtp_arrays_iterate(PyMultipleArrayIter iter, Py_ssize_t dim){
    Py_ssize_t **dims = iter.dims[dim];
    Py_ssize_t *repeat = iter.repeat[dim];
    for(Py_ssize_t j = 0; j < iter.nm; j++)
        for(Py_ssize_t i = 0; i < repeat[j]; i++)
            iter.array_iter[j].next(&iter.array_iter[j],dims[j][i]);
    return 0;
}

int multiple_arrays_iter_next(PyMultipleArrayIter *iter){
    iter->index[0]++;
    mtp_arrays_iterate(*iter,0);
    for(Py_ssize_t i = 0; i < iter->ns && iter->index[i] >= iter->shapes[i]; i++){
        if(i == iter->ns - 1)
            return 0; // last iteration
        iter->index[i] = 0;
        iter->index[i+1]++;
        iter->dim = i+1;
        mtp_arrays_iterate(*iter,i+1);
    }
    
    return 1;
}

static PyMultipleArrayIter init_artificial_mtp_arrays_iter(PyMvArrayObject *self){
    PyMultipleArrayIter iter;
    iter.array_iter = (PyMvArrayIter*)PyMem_RawMalloc(sizeof(PyMvArrayIter));
    iter.array_iter->index = (Py_ssize_t*)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t));
    iter.array_iter->data = self->data;
    iter.array_iter->strides = self->strides;
    iter.array_iter->shapes = self->shapes;
    iter.array_iter->next = multivector_array_iter_next;
    iter.array_iter->ns = self->ns;
    iter.nm = 1;
    iter.shapes = self->shapes;
    iter.dims = (Py_ssize_t***)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t**));
    iter.repeat = (Py_ssize_t**)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t*));
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t));
    for(Py_ssize_t i = 0; i < self->ns; i++){
        iter.dims[i] = (Py_ssize_t**)PyMem_RawMalloc(iter.nm*sizeof(Py_ssize_t*));
        iter.repeat[i] = (Py_ssize_t*)PyMem_RawMalloc(iter.nm*sizeof(Py_ssize_t));
        for(Py_ssize_t j = 0; j < iter.nm; j++){
            iter.dims[i][j] = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
            *(iter.dims[i][j]) = self->ns - 1 - i; // Change order of dimensions
            iter.repeat[i][j] = 1; // No repeated symbols
        }
        iter.array_iter->index[i] = 0;
        iter.index[i] = 0;
    }
    iter.dim = -1;
    iter.ns = self->ns;
    
    return iter;
}

static void free_artificial_mtp_arrays_iter(PyMultipleArrayIter iter){
    PyMem_RawFree(iter.array_iter->index);
    for(Py_ssize_t i = 0; i < iter.ns; i++){
        PyMem_RawFree(iter.repeat[i]);
        for(Py_ssize_t j = 0; j < iter.nm; j++){
            PyMem_RawFree(iter.dims[i][j]);
        }
        PyMem_RawFree(iter.dims[i]);
    }
    PyMem_RawFree(iter.index);
    PyMem_RawFree(iter.dims);
    PyMem_RawFree(iter.repeat);
    PyMem_RawFree(iter.array_iter);
}

PyMultivectorIter *init_multivector_array_iter(PyMvArrayObj *self, Py_ssize_t size){
    PyMultivectorIter *iter = (PyMultivectorIter*)PyMem_RawMalloc(size*sizeof(PyMultivectorIter));
    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;
    for(Py_ssize_t i = 0; i < size; i++)
        iter[i] = iter_init(self->data[i],self->type);
    return iter;
}

PyObject *multivector_array_repr(PyMvArrayObj *self){
    PrintTypeMV ptype = self->GA->print_type_mv;
    PyObject *out = Py_BuildValue("s","");
    char *mv_str;
    char out_str[10000];
    out_str[0] = '\0';

    const char *lbracket = "[";
    const char *rbracket_nl = "],\n";
    const char *rbracket = "]";

    Py_ssize_t last_dim = -9;
    
    //strcat(out_str,"multivector_array(");

    for(Py_ssize_t i = 0; i < self->ns; i++)
        strcat(out_str,lbracket);

    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;
    if(self->ns > 0){
        PyMultipleArrayIter arr_iter = init_artificial_mtp_arrays_iter(self);
        do{
            PyMultivectorIter iter = iter_init(*arr_iter.array_iter->data,self->type);
            mv_str = type_iter_repr_1(&iter,ptype,iter.niters);
            
            if(last_dim != arr_iter.dim && last_dim != -9){
                out_str[strlen(out_str) - 1] = ']';
                for(Py_ssize_t i = 0; i < arr_iter.dim - 1; i++)
                    strcat(out_str,rbracket);
                strcat(out_str,",");
                strcat(out_str,"\n");
                strcat(out_str,lbracket);
                for(Py_ssize_t i = 0; i < arr_iter.dim - 1; i++)
                    strcat(out_str,lbracket);
            }
            strcat(out_str,lbracket);
            strcat(out_str,mv_str);
            strcat(out_str,rbracket);
            strcat(out_str,",");
            PyMem_RawFree(mv_str);
            PyMem_RawFree(iter.index);
            
            last_dim = arr_iter.dim;
        }while(multiple_arrays_iter_next(&arr_iter));
        free_artificial_mtp_arrays_iter(arr_iter);
        
        out_str[strlen(out_str) - 1] = ']';
        for(Py_ssize_t i = 0; i < self->ns - 1; i++)
           strcat(out_str,rbracket);
        //strcat(out_str,")");
    }else{
        PyMultivectorIter iter = iter_init(*self->data,self->type);
        mv_str = type_iter_repr_1(&iter,ptype,iter.niters);
        strcpy(out_str,mv_str);
        PyMem_RawFree(mv_str);
    }

    out = Py_BuildValue("s",out_str);

    if(ptype == PrintTypeMV_normal){
        PyErr_SetString(PyExc_NotImplementedError, "Not implemented for multivector arrays");
        return NULL;
    }else if(ptype == PrintTypeMV_reduced){
        return out;
    }else {
        PyErr_SetString(PyExc_ValueError,"The selected print type is not valid");
        return NULL;
    }
}

PyTypeObject PyMultivectorArrayType = {
		PyVarObject_HEAD_INIT(NULL, 0).tp_name = "gasparse.multivector_array",
		.tp_doc = PyDoc_STR(
				"Builds an array of multivectors (PyMultivectorObject)"),
		.tp_basicsize = sizeof(PyMvArrayObj),
		.tp_itemsize = 0,
		.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
		.tp_dealloc = (destructor)multivector_array_dealloc,
		.tp_repr = (reprfunc)multivector_array_repr,
        .tp_str = (reprfunc)multivector_array_repr,
		//.tp_str = (reprfunc)multivector_repr,
		//.tp_call = (ternaryfunc)multivector_grade_project,
		.tp_new = NULL,
		//.tp_as_number = &PyMultivectorNumberMethods,
		//.tp_methods = multivector_methods
};

PyMvArrayObj *new_multivector_array(PyAlgebraObject *ga, char *type, Py_ssize_t ndims, Py_ssize_t *strides,Py_ssize_t *shapes){
    PyMvArrayObj *self = init_multivector_array(ga,ndims,strides,shapes);
    if(!self) return NULL;
    if(!get_multivector_type_table(ga, type, &self->type)) return NULL;
    return self;
}

PyObject *algebra_multivector_array(PyAlgebraObject *self, PyObject *args,PyObject *kwds) {
	static char *kwlist[] = {"values", "basis", "grades","dtype", NULL};
	PyObject *values = NULL, *basis = NULL, *grades = NULL;
	int *bitmaps_int = NULL;
	int *grades_int = NULL;
	ga_float *values_float = NULL;
	ga_float *values_basis = NULL;
	Py_ssize_t size, bsize;
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

    int basis_is_mv = 0;
	
	if (basis) {
        // Defining multivectors through some basis
		if (!PyList_Check(basis)) {
			PyErr_SetString(PyExc_ValueError,"Basis must be a list and of the same size as values");
			return NULL;
		}
        bsize = parse_list_as_multivectors_1(basis, &values_basis, &bitmaps_int);
		if (bsize <= 0){
            bsize = parse_list_as_bitmaps(basis, &bitmaps_int);
                if(bsize <= 0) {
                    PyErr_SetString(PyExc_TypeError, "Error parsing basis list as bitmaps");
                    return NULL;
                }
        }else
            basis_is_mv = 1;
        
	} else if (grades) {
		
		Py_ssize_t gsize = parse_list_as_grades(self, grades, &grades_int);
		if (gsize <= 0) {
			PyMem_RawFree(values_float);
			PyErr_SetString(PyExc_ValueError,"Error parsing grades, invalid value or empty");
			return NULL;
		}
		Py_ssize_t mv_size =
				parse_list_as_basis_grades(*self, grades_int, &bitmaps_int, gsize);

        bsize = mv_size;
		PyMem_RawFree(grades_int);

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
    
    // Determine strides and shapes from the nested lists
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
    Py_ssize_t *strides = get_strides(ndims-1,shapes);
    if(!strides){
        PyErr_SetString(PyExc_ValueError,"Error geting strides");
        return NULL;
    }
    
    ga_float **values_float_array = (ga_float**)PyMem_RawMalloc(strides[0]*sizeof(ga_float*));
    
    // Get the values from the nested list
    if(iterate_nested_lists_1(values,values_float_array,strides,shapes,0,0,ndims-1,bsize) == -1) {
        PyErr_SetString(PyExc_ValueError,"Error iterating nested lists");
        return NULL;
    }

    Py_ssize_t *shapes_ = NULL;
    if(ndims > 1){
        shapes_ =  PyMem_RawMalloc((ndims-1)*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < ndims-1; i++) // Discard the innermost shape
            shapes_[i] = shapes[i];
    }
    PyMvArrayObj *mv_array = new_multivector_array(self,type_name,ndims-1,strides,shapes_);
    if(!mv_array){
        PyErr_SetString(PyExc_ValueError,"Error creating new multivector array");
        return NULL;
    }
    
    gainitfunc init = mv_array->type->data_funcs->init;
    if (!init){
        PyMem_RawFree(values_float_array);
        PyMem_RawFree(bitmaps_int);
        return NULL; // raise not implemented error
    }

    // If the basis is a multivector list multiply each basis element by the values array. 
    if(basis_is_mv){
        for(Py_ssize_t i = 0; i < strides[0]; i++){
            for(Py_ssize_t j = 0;j < bsize; j++){
                values_float_array[i][j] *= values_basis[j];
            }
        }
        PyMem_RawFree(values_basis);
    }
    mv_array->data = (void**)PyMem_RawMalloc(strides[0]*sizeof(void*));

    for(Py_ssize_t i = 0; i < strides[0]; i++)
        mv_array->data[i] = NULL;
    
    for(Py_ssize_t i = 0; i < strides[0]; i++){
        mv_array->data[i] = init(bitmaps_int,values_float_array[i],bsize,self);
        if(!mv_array->data[i]){
            // dealloc the memory from all the initialized objects
            multivector_array_dealloc(mv_array);
            PyErr_SetString(PyExc_ValueError,"Error initializing a single multivector!");
            mv_array = NULL;
            break;
        }
    }
    
    for(Py_ssize_t i = 0; i < strides[0]; i++)
        PyMem_RawFree(values_float_array[i]);
    PyMem_RawFree(values_float_array);
    PyMem_RawFree(bitmaps_int);
    PyMem_RawFree(shapes);

    return (PyObject*)mv_array;
}


/*

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
*/


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

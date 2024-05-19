#include "multivector_object.h"
#include "multivector_types.h"
#include "common.h"
#include "types.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define RECTIFIER(x) ((x<0) ? 0 : x)

// Multivector Array Initializers
static PyMvObject *init_multivector_array(PyAlgebraObject *GA, Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes){
	if(!GA) return NULL; 
    PyMvObject *array_obj = (PyMvObject*)PyMem_RawMalloc(sizeof(PyMvObject));
    if(!array_obj) return NULL;

    // Copy strides and copy shapes if not null
    if(shapes != NULL){
        array_obj->shapes = (Py_ssize_t*)PyMem_RawMalloc(ndims*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < ndims; i++)
            array_obj->shapes[i] = shapes[i];
    }else
        array_obj->shapes = NULL;
    if(strides != NULL){
        array_obj->strides = (Py_ssize_t*)PyMem_RawMalloc((ndims+1)*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < ndims + 1; i++)
            array_obj->strides[i] = strides[i];
    }else{
        array_obj->strides = NULL;
    }

	array_obj->data = NULL;
	array_obj->ns = ndims;
	// set type and increase reference count
    array_obj->GA = GA;
    // Set the mixed type operations table
    array_obj->mixed = GA->mixed;

    Py_XINCREF((PyObject*)array_obj->GA);
	Py_SET_TYPE(array_obj, &PyMultivectorType);
	Py_SET_REFCNT((PyObject*)array_obj,1);

	return array_obj;
}

static int alloc_mvarray_data(PyMultivectorObject *obj){
    // Alloc memory for the data
    if(obj->strides != NULL){
        obj->data = (void*)PyMem_RawMalloc(obj->strides[0]*obj->type->basic_size);
        if(!obj->data) return 0;

        //Initializing each multivector to NULL otherwise freeing can cause trouble
        gainitfunc init = obj->type->data_funcs->init;
        if(!init) return 0;
        for(Py_ssize_t i = 0; i < obj->strides[0]; i++){
            if(!init(obj->data + i*obj->type->basic_size,obj->GA,NULL,NULL,0))
                return 0;
        }
    }
    return 1;
}

PyMultivectorObject *new_multivector_array(PyAlgebraObject *GA, char *type,  Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes){
    PyMultivectorObject *self = init_multivector_array(GA, ndims, strides, shapes);
    if(!self) return NULL;
    if(!get_multivector_type_table(GA, type, &self->type)) return NULL;
    if(!alloc_mvarray_data(self)) return NULL;
    return self;
}

PyMultivectorObject *new_mvarray_inherit_type(PyAlgebraObject *GA,  Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes, PyMultivectorSubType *type){
    PyMultivectorObject *self = init_multivector_array(GA, ndims, strides, shapes);
    if(!self) return NULL;
    self->type = type;
    if(!alloc_mvarray_data(self)) return NULL;
    return self;
}

// Single multivector initializers and destructors
PyMultivectorObject *init_multivector(PyAlgebraObject *GA){
    // Initializes a single multivector, also allocs memory for the data
    Py_ssize_t strides = 1;
    PyMultivectorObject *self = init_multivector_array(GA,0,&strides,NULL);
    if(!self) return NULL;
    return self;
}

PyMultivectorObject *new_multivector(PyAlgebraObject *GA, char *type){
    PyMultivectorObject *self = init_multivector(GA);
    if(!self) return NULL;
    if(!get_multivector_type_table(GA, type, &self->type)) return NULL;
    if(!alloc_mvarray_data(self)) return NULL;
    return self;
}

PyMultivectorObject *new_multivector_inherit_type(PyAlgebraObject *GA, PyMultivectorSubType *type){
    PyMultivectorObject *self = init_multivector(GA);
    if(!self || !type) return NULL;
    self->type = type;
    // Can only allocate memory after type is set
    if(!alloc_mvarray_data(self)) return NULL;
    return self;
}

// Creates a multivector array of the same shape as the mvarray
PyMultivectorObject *new_mvarray_from_mvarray(PyMvObject *mvarray){
    PyMultivectorObject *self = init_multivector_array(mvarray->GA, mvarray->ns, mvarray->strides, mvarray->shapes);
    if(!self) return NULL;
    self->type = mvarray->type;
    if(!alloc_mvarray_data(self)) return NULL;
    return self;
}

static void multivector_array_dealloc(PyMvObject *self){
	void *data = self->data;
    gafreefunc free_type = self->type->data_funcs->free;
    if(free_type){
        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            free_type(self->data + i*self->type->basic_size);
        }
    }

	Py_XDECREF((PyObject*)self->GA);
	PyMem_RawFree(self->strides);
	PyMem_RawFree(self->shapes);
	PyMem_RawFree(data);
	PyMem_RawFree(self);
}

PyMultivectorIter *init_multivector_array_iters(PyMvObject *self){
    PyMultivectorIter *iters = (PyMultivectorIter*)PyMem_RawMalloc(self->strides[0]*sizeof(PyMultivectorIter));
    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;
    for(Py_ssize_t i = 0; i < self->strides[0]; i++){
        iters[i] = iter_init(self->data + i*self->type->basic_size, self->type);
    }
    return iters;
}

void free_multivector_array_iter(PyMultivectorIter *iters, Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        PyMem_RawFree(iters[i].index);
    }
    PyMem_RawFree(iters);
}

int cast_mvarray(PyMvObject *from, PyMvObject *to){
    PyMultivectorIter *iters = init_multivector_array_iters(from);
    gacastfunc cast = to->type->data_funcs->cast;
    for(Py_ssize_t i = 0; i < from->strides[0]; i++){
        if(!cast(iters+i,to->data + i*to->type->basic_size,to->GA)){
            free_multivector_array_iter(iters,from->strides[0]);
            return 0;
        }
    }
    free_multivector_array_iter(iters,from->strides[0]);
    return 1;
}

PyMvObject *cast_mvarray_inherit_type(PyMvObject *from, PyMultivectorSubType *type){
    PyMvObject *to = new_mvarray_inherit_type(from->GA,from->ns,from->strides,from->shapes,type);
    if(!cast_mvarray(from,to)){
        multivector_array_dealloc(to);
        return NULL;
    }
    return to;
}

PyMvObject *cast_mvarray_to_type(PyMvObject *from, char *type){
    PyMvObject *to = new_multivector_array(from->GA,type,from->ns,from->strides,from->shapes);
    if(!cast_mvarray(from,to)){
        multivector_array_dealloc(to);
        return NULL;
    }
    return to;
}


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

static int iterate_nested_lists(PyObject *list,
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
        // return error if the innermost list is not of the same size of the bitmaps/grades 
        if(size <= 0 || size != nbasis) return -1;
    }else{
        for(Py_ssize_t i = 0; i < (Py_ssize_t)PyList_Size(list); i++){
			PyObject *sublist = PyList_GetItem(list, i);
			int flag = iterate_nested_lists(sublist, array, strides, shape, index + i*strides[dim+1], dim+1,ndims,nbasis);
			if(flag == -1) return -1;
		}
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


int multiple_arrays_iter_next(PyMultipleArrayIter *iter){
    // Get out of the iterator imediatly. If index is null only need one iteration.
    if(!iter->index) return 0;
    iter->index[iter->ns-1]++;
    iter->dflag = 0; // reset the dim flag
    for(Py_ssize_t i = iter->ns-1; i >= 0 && iter->index[i] >= iter->shapes[i]; i--){
        if(i == 0) return 0; // last iteration
        iter->index[i] = 0;
        iter->index[i-1]++;
        iter->dim = i-1;
        iter->dflag = 1; // set the dim flag
    }

    for(Py_ssize_t j = 0; j < iter->nm; j++){
        void *data = iter->arrays[j].data0;
        Py_ssize_t index = 0;
        for(Py_ssize_t i = 0; i < iter->ns; i++)
            index += iter->arrays[j].strides[i+1]*iter->index[i];
        iter->arrays[j].data = data + index*iter->arrays[j].basic_size;
    }
    return 1;
}

static int set_listsoflists_element(PyObject *element, PyObject *list, Py_ssize_t *index, Py_ssize_t *shape, Py_ssize_t size){
    PyObject *sublist = list;
    PyObject *subsublist;
    for(Py_ssize_t i = 0; i < size-1; i++){ 
        subsublist = PyList_GetItem(sublist, index[i]);
        if(subsublist == NULL){ // when empty list 
            // Create list at index[i]
            subsublist = PyList_New(shape[i+1]); 
            PyList_SetItem(sublist,index[i],subsublist);
        }
        sublist = subsublist;
    }
    PyList_SetItem(sublist,index[size-1],element); // Last index
    return 1;
}

static int check_arrays(Py_ssize_t *arr, Py_ssize_t item,Py_ssize_t size){
    for(Py_ssize_t i = 0; i < size; i++){
        if(arr[i] == item) return 1;
    }
    return 0;
}

// Given the index array initialize the iterator
static PyMultipleArrayIter init_arrays_iter_index(PyMvObject *self, Py_ssize_t *pos, Py_ssize_t *index, Py_ssize_t size, PyMvObject **other){
    PyMultipleArrayIter iter = {.arrays = NULL, .index = NULL, .nm = 2, .dflag = 0};

    // check if indices are valid
    for(Py_ssize_t i = 0; i < size; i++)
        if(index[i] >= self->shapes[pos[i]])
            return iter;

    iter.arrays = (PyMvBasicArray*)PyMem_RawMalloc(2*sizeof(PyMvBasicArray));
    
    // initiallize the arrays according to the index array
    iter.arrays->data0 = self->data;
    iter.arrays->basic_size = self->type->basic_size;

    for(Py_ssize_t i = 0; i < size; i++)
        iter.arrays->data0 += index[i]*self->strides[pos[i] + 1]*iter.arrays->basic_size;
    
    iter.arrays->data = iter.arrays->data0;
    iter.arrays->ns = self->ns - size;
    iter.ns = self->ns - size;

    iter.arrays->strides = (Py_ssize_t*)PyMem_RawMalloc((iter.arrays->ns + 1)*sizeof(Py_ssize_t));
    
    iter.arrays->strides[0] = self->strides[0];
    Py_ssize_t j = 1;
    for(Py_ssize_t i = 1; i < self->ns + 1; i++){
        while(i <= self->ns + 1 && check_arrays(pos,i-1,size)) i++; // skip when i equals some pos[k] + 1
        if(i >= self->ns + 1) break; // check if end is reached
        iter.arrays->strides[j] = self->strides[i];
        j++;
    }
    
    
    iter.shapes = (Py_ssize_t*)PyMem_RawMalloc(iter.arrays->ns*sizeof(Py_ssize_t));

    j = 0;
    for(Py_ssize_t i = 0; i < self->ns; i++){
        while(i <= self->ns && check_arrays(pos,i,size)) i++; // skip when i equals some pos[k]
        if(i >= self->ns) break;
        iter.shapes[j] = self->shapes[i];
        j++;
    }

    if(iter.arrays->ns > 0){ // Only allocate memory if the size is greater than zero
        iter.index = (Py_ssize_t*)PyMem_RawMalloc(iter.arrays->ns*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < iter.arrays->ns; i++){
            iter.index[i] = 0;
        }
    }
    iter.dim = -1;

    // From the shapes of the iterator determine the strides for the new multivector
    Py_ssize_t *strides_other = get_strides(iter.arrays->ns,iter.shapes);
    
    // Allocate memory for the multivector array
    *other = new_multivector_array(self->GA, self->type->type_name, iter.arrays->ns, strides_other, iter.shapes);
    iter.arrays[1].data0 = (*other)->data;
    iter.arrays[1].data = (*other)->data;
    iter.arrays[1].basic_size = (*other)->type->basic_size;
    iter.arrays[1].ns = (*other)->ns;
    iter.arrays[1].strides = (*other)->strides;

    PyMem_RawFree(strides_other);

    return iter;
}

static void free_arrays_iter_index(PyMultipleArrayIter iter){
    PyMem_RawFree(iter.arrays->strides);
    PyMem_RawFree(iter.index);
    PyMem_RawFree(iter.arrays);
    PyMem_RawFree(iter.shapes);
}

static PyMultipleArrayIter init_single_array_iter(PyMvObject *self){
    PyMultipleArrayIter iter = {.dflag = 0};
    iter.arrays = (PyMvBasicArray*)PyMem_RawMalloc(sizeof(PyMvBasicArray));
    iter.arrays->data = self->data;
    iter.arrays->data0 = self->data;
    iter.arrays->strides = (Py_ssize_t*)PyMem_RawMalloc((self->ns+1)*sizeof(Py_ssize_t));
    
    for(Py_ssize_t i = 0; i < self->ns + 1; i++)
        iter.arrays->strides[i] = self->strides[i];
    
    iter.arrays->ns = self->ns;
    iter.ns = self->ns;
    iter.arrays->basic_size = self->type->basic_size;
    iter.nm = 1; // The number of arrays
    iter.shapes = self->shapes;
    iter.index = (Py_ssize_t*)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t));
    for(Py_ssize_t i = 0; i < self->ns; i++){
        iter.index[i] = 0;
    }
    iter.dim = -1;
    
    return iter;
}

static void free_multiple_arrays_iter(PyMultipleArrayIter iter){
    for(Py_ssize_t i = 0; i < iter.nm; i++)
        PyMem_RawFree(iter.arrays[i].strides);
    
    PyMem_RawFree(iter.arrays);
    PyMem_RawFree(iter.index);
}

static int get_value_bitmap_from_mv(PyMultivectorObject *data, ga_float *value, int *bitmap) {
	// Takes a multivector and writes the first non-zero element into value and bitmap

	PyMultivectorIter *iter = init_multivector_iter(data, 1);
	if (!iter)
		return -1;
	while(iter->next(iter)){
		*value = iter->value;
		*bitmap = iter->bitmap;
		if(*value != 0.0)
			break;
	}
	
	free_multivector_iter(iter, 1);
	return 1;
}

int parse_list_as_multivectors(PyObject *basis, ga_float **values, int **bitmaps){
    if (!PyList_Check(basis))
		return -1;
	Py_ssize_t size = PyList_Size(basis);
    *bitmaps = (int *)PyMem_RawMalloc(size * sizeof(int));
	*values = (ga_float *)PyMem_RawMalloc(size * sizeof(ga_float));
    for (Py_ssize_t i = 0; i < size; i++) {
		PyObject *basis_i = PyList_GetItem(basis, i);

		if (Py_IS_TYPE(basis_i, &PyMultivectorType)) {
			if (!get_value_bitmap_from_mv((PyMvObject *)basis_i,&(*values)[i], &(*bitmaps)[i])) {
				PyMem_RawFree(*bitmaps);
				PyMem_RawFree(*values);
				return -1;
			}
        }else if (PyFloat_Check(basis_i)) {
			(*values)[i] = (ga_float)PyFloat_AsDouble(basis_i);
			(*bitmaps)[i] = 0;
		} else if (PyLong_Check(basis_i)) {
			(*values)[i] = (ga_float)PyLong_AsLong(basis_i);
			(*bitmaps)[i] = 0;
		} else {
			PyMem_RawFree(*bitmaps);
			PyMem_RawFree(*values);
			return -1;
		}
	}
	return size;
}

// Parses a list of lists and converts it to an array of specified types
PyObject *algebra_multivector(PyAlgebraObject *self, PyObject *args, PyObject *kwds) {
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
        bsize = parse_list_as_multivectors(basis, &values_basis, &bitmaps_int);
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

    
    if(!strcmp(type_name,"scalar")){
        // If the size of the basis is not one and the basis is not a scalar
        if(!(bsize == 1 && !bitmaps_int[0])){
            PyErr_SetString(PyExc_ValueError, "The scalar type multivector needs a scalar basis!!");
            PyMem_RawFree(bitmaps_int);
            return NULL;
        }
    }
    
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
    if(iterate_nested_lists(values,values_float_array,strides,shapes,0,0,ndims-1,bsize) == -1) {
        PyErr_SetString(PyExc_ValueError,"Error iterating nested lists, shape of the list might not be correct!!!");
        return NULL;
    }

    // Copies shapes up to i=ndims-2. Copies strides up to i=ndims-1
    PyMvObject *mv_array = new_multivector_array(self,type_name,ndims-1,strides,shapes);
    if(!mv_array){
        PyErr_SetString(PyExc_ValueError,"Error creating new multivector array");
        return NULL;
    }
    
    gainitfunc init = mv_array->type->data_funcs->init;
    if (!init){
        PyMem_RawFree(values_float_array);
        PyMem_RawFree(bitmaps_int);
        PyErr_SetString(PyExc_NotImplementedError,"The initializer for this type is not implemented");
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
    
    // Write the data into the multivector array
    for(Py_ssize_t i = 0; i < strides[0]; i++){
        if(!init(INDEX_DATA(mv_array,i),self,bitmaps_int,values_float_array[i],bsize)){
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
    PyMem_RawFree(strides);

    return (PyObject*)mv_array;
}



char *type_iter_repr(PyMultivectorIter *iter, Py_ssize_t dsize){
    char *out_str;
    // if(ptype == PrintTypeMV_reduced){
    if(dsize){
        char **str_blade = (char**)PyMem_RawMalloc(dsize*sizeof(char*));
        Py_ssize_t len = 0;
        char sep[] = " + ";

        Py_ssize_t i = 0;
        while(iter->next(iter)){
            // only skip small values if its not sparse or blades type
            // if(iter->type != MultivectorType_sparse && iter->type != MultivectorType_blades)
                // maybe should add this as an option when creating the algebra
                // if(ABS(iter->value) < 1e-12) continue; // don't print small values

            char *value = PyOS_double_to_string((double)iter->value,'g',8,0,NULL);
            if(iter->bitmap){
                char *bitmap = bitmap_to_string(iter->bitmap);
                Py_ssize_t size = strlen(bitmap) + strlen(value)+3;
                str_blade[i] = (char*)PyMem_RawMalloc(size*sizeof(char));
                PyOS_snprintf(str_blade[i],size,"%s*%s",value,bitmap);
                PyMem_RawFree(bitmap);
                PyMem_Free(value);
            }else{
                Py_ssize_t size = strlen(value) + 1;
                str_blade[i] = (char*)PyMem_RawMalloc(size*sizeof(char));
                strcpy(str_blade[i],value);
                PyMem_Free(value);
            }
            len += strlen(str_blade[i]);
            i++;
        }
        if(iter->type != MultivectorType_sparse && iter->type != MultivectorType_blades)
            dsize = i;
        if(!dsize){
            PyMem_RawFree(str_blade);
            out_str = PyMem_RawMalloc(4*sizeof(char));
            out_str[3] = '\0';
            strcpy(out_str,"0.0");
            return out_str;
        }
        len += dsize*3 + 3;
        out_str = (char*)PyMem_RawMalloc(len*sizeof(char));

        Py_ssize_t j = 0;
        for(i = 0; i < dsize-1; i++){
            strcpy(out_str + j, str_blade[i]);
            j += strlen(str_blade[i]);
            strcpy(out_str + j, sep);
            j += strlen(sep);
        }
        strcpy(out_str + j, str_blade[i]);
        j += strlen(str_blade[i]);

        for(Py_ssize_t i = 0; i < dsize; i++)
            PyMem_RawFree(str_blade[i]);
        PyMem_RawFree(str_blade);

        return out_str;
    }else{ 
        out_str = PyMem_RawMalloc(4*sizeof(char));
        out_str[3] = '\0';
        strcpy(out_str,"0.0");
        return out_str;
    }
    // }else if(ptype == PrintTypeMV_normal){
        /*
        if(dsize){
            char **str_bitmap = (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            char **str_value =  (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            Py_ssize_t len_bitmap = 0, len_value = 0;
            Py_ssize_t i = 0;
            while(iter->next(iter)){
                str_bitmap[i] = bitmap_to_string(iter->bitmap);
                len_bitmap += strlen(str_bitmap[i]);
                str_value[i] = PyOS_double_to_string((double)iter->value,'g',6,0,NULL);
                len_value += strlen(str_value[i]);
                i++;
            }
            len_value += dsize + 3;
            len_bitmap += dsize + 3;
            char *format_value = (char*)PyMem_RawMalloc((len_value)*sizeof(char));
            char *format_bitmap = (char*)PyMem_RawMalloc((len_bitmap)*sizeof(char));

            Py_ssize_t j = 0, k = 0; i = 0;
            for(i = 0; i < dsize - 1; i++){
                strcpy(format_bitmap + j, str_bitmap[i]);
                j += strlen(str_bitmap[i]);
                format_bitmap[j++] = ',';

                strcpy(format_value + k, str_value[i]);
                k += strlen(str_value[i]);
                format_value[k++] = ',';
            }
            strcpy(format_bitmap + j, str_bitmap[i]);
            j += strlen(str_bitmap[i]);
            format_bitmap[j] = '\0';


            strcpy(format_value + k, str_value[i]);
            k += strlen(str_value[i]);
            format_value[k] = '\0';

            char *format = ".multivector([%s],blades=[%s],dtype=%s)";
            size_t size = len_value + len_bitmap + strlen(format);
            char *format_out = (char*)PyMem_RawMalloc(size*sizeof(char));

            PyOS_snprintf(format_out,size,format,format_value,format_bitmap,iter->type_name);

            for(Py_ssize_t i = 0; i < dsize; i++){
                PyMem_RawFree(str_bitmap[i]);
                PyMem_Free(str_value[i]);
            }
            PyMem_RawFree(str_bitmap);
            PyMem_RawFree(str_value);
            PyMem_RawFree(format_value);
            PyMem_RawFree(format_bitmap);
            return format_out;
        }else {
            const char* zero = ".multivector(0.0)";
            out_str = PyMem_RawMalloc((strlen(zero)+1)*sizeof(char));
            out_str[strlen(zero)] = '\0';
            strcpy(out_str,zero);
            return out_str;
        }*/
    // }
    return NULL; // raise error
}

/*
  Method calls of PyMultivectorType:
    print(a) -> multivector_repr(self): prints the multivector
    a(2)     -> multivector_grade_project(self): grade projection
  Method calls for the number protocol:
    a*b -> multivector_geometric_product: geometric product
    a^b -> multivector_outer_product: wedge product
    a|b -> multivector_inner_product: inner product
    a+b -> multivector_add: multivector addition
    a-b -> multivector_subtract: multivector subtraction
    ~a  -> multivector_invert: reverse the order of the basis vector by changing the sign
    a/b -> multivector_divide: Divide by scalar or scalar arrays
*/

static void repeat_string(const char *rep,char *str,Py_ssize_t n){
    for(Py_ssize_t i = 0; i < n; i++)
        strcat(str,rep);
}

// Prints the array of multivectors as a list of lists
PyObject *multivector_repr(PyMvObject *self){
    PrintTypeMV ptype = self->GA->print_type_mv;
    PyObject *out;

    Py_ssize_t out_size = RECTIFIER(self->ns) + 2;
    Py_ssize_t n_spaces = self->ns;

    const char *lbracket = "[";
    const char *rbracket = "]";

    char *mv_str;
    char *out_str = (char*)PyMem_RawMalloc(out_size*sizeof(char));
    char fflag = 1; // a first iteration flag
    
    out_str[0] = '\0';

    repeat_string(lbracket, out_str,self->ns);
    
    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;

    if(self->ns > 0){
        PyMultipleArrayIter arr_iter = init_single_array_iter(self);
        out_size += RECTIFIER(self->ns);
        
        do{            
                
            Py_ssize_t n =  arr_iter.ns - arr_iter.dim - 1;

            PyMultivectorIter iter = iter_init(arr_iter.arrays->data,self->type);
            mv_str = type_iter_repr(&iter,iter.niters);
            out_size += strlen(mv_str) + 2;
            out_size += 1;

            // the start of a new dimension
            if(arr_iter.dflag){
                
                out_size += 3*RECTIFIER(n) + 2 + RECTIFIER(n_spaces-n);
                out_str = (char*)PyMem_RawRealloc(out_str, out_size); // Realloc memory for the output string
                
                repeat_string(rbracket, out_str,n);
                strcat(out_str,",\n");
                repeat_string("\n", out_str,n);
                repeat_string(" ", out_str, n_spaces-n);
                repeat_string(lbracket, out_str,n);

            }else if(!fflag){ // Not the first iteration
                out_size += 2 + RECTIFIER(n_spaces);
                out_str = (char*)PyMem_RawRealloc(out_str, out_size); // Realloc memory for the output string

                strcat(out_str,",\n");
                repeat_string(" ", out_str, n_spaces); // add n spaces to out_str
            }else{// the first iteration
                out_str = (char*)PyMem_RawRealloc(out_str, out_size);
            }
            strcat(out_str,lbracket);
            strcat(out_str,mv_str);
            strcat(out_str,rbracket);

            fflag = 0;
            // char dims_str[10];
            // sprintf(dims_str,"%zd",arr_iter.dim);
            // strcat(out_str,dims_str);

            PyMem_RawFree(mv_str);
            PyMem_RawFree(iter.index);
            
            
        }while(multiple_arrays_iter_next(&arr_iter));

        free_multiple_arrays_iter(arr_iter);
        repeat_string(rbracket, out_str,self->ns);

    }else{
        if(self->strides[0] == 1){
            PyMultivectorIter iter = iter_init(self->data,self->type);
            mv_str = type_iter_repr(&iter,iter.niters);
            
            // Realloc memory for the output string
            out_size += strlen(mv_str);
            out_str = (char*)PyMem_RawRealloc(out_str, out_size);
            
            strcpy(out_str,mv_str);
            PyMem_RawFree(mv_str);
            PyMem_RawFree(iter.index);
        }
    }
    
    out = Py_BuildValue("s",out_str);
    PyMem_RawFree(out_str);
    
    if(ptype == PrintTypeMV_normal){
        PyObject *mv_arr_str;
        PyObject *mv_arr_end;

        if(self->ns > 0 ){
            mv_arr_str = Py_BuildValue("s",".mvarray(\n");
            mv_arr_end = Py_BuildValue("s","\n)");
        }else{
            mv_arr_str = Py_BuildValue("s",".mvarray(");
            mv_arr_end = Py_BuildValue("s",")");
        }
        PyObject *ga_repr = PyObject_Repr((PyObject*)self->GA);
        PyObject *out1 = PyUnicode_Concat(ga_repr, mv_arr_str);
        out1 = PyUnicode_Concat(out1, out);
        out1 = PyUnicode_Concat(out1, mv_arr_end);

        Py_XDECREF(ga_repr);
        Py_XDECREF(mv_arr_str);
        Py_XDECREF(mv_arr_end);
        return out1;

    }else if(ptype == PrintTypeMV_reduced){
        return out;
    }else {
        PyErr_SetString(PyExc_ValueError,"The selected print type is not valid");
        return NULL;
    }


}

static int get_scalar(PyObject *self, ga_float *value){
    if(PyFloat_Check(self)){
        *value = (ga_float)PyFloat_AsDouble(self);
        return 1;
    }
    if(PyLong_Check(self)){
        *value = (ga_float)PyLong_AsDouble(self);
        return 1;
    }
    return 0;
}


PyObject *list_from_mvarray(PyMvObject *dense, Py_ssize_t *grade_bool, Py_ssize_t size){
    gaiterinitfunc iter_init = dense->type->data_funcs->iter_init;
    PyMultipleArrayIter arr_iter = init_single_array_iter(dense);
    PyObject * values_list = PyList_New(dense->shapes[0]);

    do{
        PyMultivectorIter iter = iter_init(arr_iter.arrays->data,dense->type);
        Py_ssize_t j = 0;
        PyObject *element = PyList_New(size);
        while(iter.next(&iter)){
            if(grade_bool[GRADE(iter.bitmap)] && j < size){
                PyObject *value = PyFloat_FromDouble(iter.value);
                PyList_SetItem(element,j,value);
                j++;
            }
            
        } PyMem_RawFree(iter.index);
        set_listsoflists_element(element, values_list, arr_iter.index, arr_iter.shapes, arr_iter.ns);
       
    }while(multiple_arrays_iter_next(&arr_iter));

	free_multiple_arrays_iter(arr_iter);
    return values_list;
}

PyObject *grade_from_multivector(PyMultivectorIter iter){
    PyObject *element;
    int grade = -1;
    while(iter.next(&iter)){
            if(grade == -1){ //First iteration
                if(iter.value != 0.0)
                    grade = iter.grade;
            }else if(grade != iter.grade){
                if(iter.value != 0.0){
                    PyMem_RawFree(iter.index);
                    return PyLong_FromLong(-1);
                }
            }
        }if(grade == -1) grade = 0;
        element = PyLong_FromLong(grade);
        PyMem_RawFree(iter.index);
        return element;
}

PyObject *grade_list_from_mvarray(PyMvObject *self){
    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;
    PyMultipleArrayIter arr_iter = init_single_array_iter(self);
    PyObject * values_list = PyList_New(self->shapes[0]);
    PyObject *element;
    do{
        PyMultivectorIter iter = iter_init(arr_iter.arrays->data,self->type);
        element = grade_from_multivector(iter);
        set_listsoflists_element(element, values_list, arr_iter.index, arr_iter.shapes, arr_iter.ns);
    }while(multiple_arrays_iter_next(&arr_iter));

	free_multiple_arrays_iter(arr_iter);
    return values_list;
}

PyObject* multivector_list_(PyMvObject *self, PyObject *args, int as_bitmap){
    int *grades_int = NULL;
    PyObject *list = NULL;
    PyObject *values_list = NULL;
    PyObject *bitmap = NULL;
    PyMultivectorObject *dense;
    PyMultivectorIter iter;
	Py_ssize_t *grade_bool = NULL;
	Py_ssize_t size;
    int free_mv = 1;

    // Cast to dense if the multivector is not dense or is not a generated type
    if(!self->type->generated && strcmp(self->type->type_name,"dense")){

        dense = new_multivector_array(self->GA,"dense",self->ns,self->strides,self->shapes);

        if(!dense){
            PyErr_SetString(PyExc_TypeError,"Error populating types table");
            return NULL;
        }
        if(!cast_mvarray(self, dense)){
            PyErr_SetString(PyExc_TypeError,"Error casting the multivector");
            return NULL;    
        }
    }else {
        free_mv = 0;
        dense = self;
    }
    
    if((size = parse_arguments_as_grades(self->GA, args, &grades_int)) > 0){ 
        // if the grades are valid and non empty
        grade_bool = get_grade_bool(grades_int, size, MAX_GRADE(self->GA) + 1);
        PyMem_RawFree(grades_int);
        size = self->GA->asize;
        Py_ssize_t psize = 0;
        for (Py_ssize_t i = 0; i < size; i++) {
            if (grade_bool[GRADE(i)])
                psize++;
        }
        size = psize;
    }else{
        size = self->GA->asize;
        grade_bool  = (Py_ssize_t*)PyMem_RawMalloc((MAX_GRADE(self->GA) + 1)*sizeof(Py_ssize_t));
		for(Py_ssize_t i = 0; i < MAX_GRADE(self->GA) + 1; i++) 
			grade_bool[i] = 1;
    }
    if(dense->strides[0] == 1)
        values_list = PyList_New(size);

    bitmap = PyList_New(size);
    iter = dense->type->data_funcs->iter_init(dense->data,dense->type);
	
	Py_ssize_t j = 0;
	ga_float basis_value = 1;
    while(iter.next(&iter)){
        if(grade_bool[GRADE(iter.bitmap)] && j < size){
            if(dense->strides[0] == 1){
                PyObject *value = PyFloat_FromDouble(iter.value);
			    PyList_SetItem(values_list,j,value);
            }
			if(as_bitmap){
				PyObject *bitmap_obj = PyLong_FromLong(iter.bitmap);
            	PyList_SetItem(bitmap,j,bitmap_obj);
			}else{
				PyMultivectorObject *mv = new_multivector(self->GA,self->type->type_name);
				if(!mv){
					// free memory
					Py_XDECREF(list);
					Py_XDECREF(bitmap);
					PyErr_SetString(PyExc_TypeError, "Cannot populate the types");
					return NULL; 
				}
				gainitfunc init = mv->type->data_funcs->init;
				init(mv->data,self->GA,&iter.bitmap,&basis_value,1);
				PyList_SetItem(bitmap,j,(PyObject*)mv);
			}
			j++;
		}else if (j > size) {
        	break;
      	}
    }
    PyMem_RawFree(iter.index);

    
    if(dense->strides[0] > 1)
        values_list = list_from_mvarray(dense,grade_bool,size);
    

    PyObject *tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple,0,values_list);
    PyTuple_SetItem(tuple,1,bitmap);
	PyMem_RawFree(grade_bool);
    if(free_mv)
        Py_XDECREF((PyObject*)dense);

    return tuple;
}

PyObject* multivector_list(PyMvObject *self, PyObject *args){
    return multivector_list_(self,args,0);
}

PyObject* multivector_list_as_bitmap(PyMvObject *self, PyObject *args){
    return multivector_list_(self,args,1);
}


PyObject* multivector_grade(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    PyObject *element = NULL;
    if(self->strides[0] == 1){
        PyMultivectorIter iter = self->type->data_funcs->iter_init(self->data,self->type);
        element = grade_from_multivector(iter);
    } else {
        element = grade_list_from_mvarray(self);
    }
    return element;
}


static PyMvObject *multivector_mixed_product(PyMvObject *left, PyMvObject *right, ProductType ptype, int isleft_single){
    int isleft_bigger = is_bigger_metric(left->GA,right->GA);
    PyAlgebraObject *GA = NULL;
    PyMvObject *out = NULL;
    Py_ssize_t left_inc = 1,right_inc = 1;
    Py_ssize_t size = -1;
    gamixedprodfunc mixed_product = NULL;

    if(isleft_bigger == -1){
        return NULL;
    }else {
        gaiterinitfunc iter_init_left = left->type->data_funcs->iter_init;
        gaiterinitfunc iter_init_right = right->type->data_funcs->iter_init;
        Py_ssize_t *strides, *shapes, ns;
        if(isleft_single == 1){
            strides = right->strides;
            shapes = right->shapes;
            ns = right->ns;
            size = right->strides[0];
            left_inc = 0;
        }else if(isleft_single == 0){
            strides = left->strides;
            shapes = left->shapes;
            ns = left->ns;
            size = left->strides[0];
            right_inc = 0;
        }else{
            strides = left->strides;
            shapes = left->shapes;
            ns = left->ns;
            size = left->strides[0];
        }
        
        if(isleft_bigger){
            out = new_multivector_array(left->GA, "sparse",ns,strides,shapes);
            GA = left->GA;
            mixed_product = left->mixed->product;
        }
        else{
            out = new_multivector_array(right->GA,"sparse",ns,strides,shapes);
            GA = right->GA;
            mixed_product = right->mixed->product;
        }

        if(mixed_product){
            for(Py_ssize_t i = 0; i < size; i++){
                PyMultivectorIter iter_left = iter_init_left(left->data + i*left->type->basic_size*left_inc,left->type);
                PyMultivectorIter iter_right = iter_init_right(right->data + + i*right->type->basic_size*right_inc,right->type);
                
                if(!mixed_product(out->data + i*out->type->basic_size,&iter_left,&iter_right,GA,ptype)){
                    return NULL;
                }
                PyMem_RawFree(iter_left.index);
                PyMem_RawFree(iter_right.index);
            }

            return out;
        }else 
            return NULL;
        
    }
} 

static PyMvObject *multivector_casttype_product(PyMvObject *left, PyMvObject *right, ProductType ptype, int isleft_single){
    PyMvObject *left_cast = NULL;
    PyMvObject *right_cast = NULL;
    PyMvObject *out = NULL;
    PyAlgebraObject *GA = NULL;
    gaprodfunc product = NULL;
    int cast_left = -1;
    Py_ssize_t left_inc = 1,right_inc = 1;
    Py_ssize_t size = -1;

    Py_ssize_t *strides, *shapes, ns;
    if(isleft_single == 1){ // The left mvarray is only one mv
        strides = right->strides;
        shapes = right->shapes;
        ns = right->ns;
        size = right->strides[0];
        left_inc = 0; // Do not increment the left array indices
    }else if(isleft_single == 0){ // The right mvarray is only one mv
        strides = left->strides;
        shapes = left->shapes;
        ns = left->ns;
        size = left->strides[0];
        right_inc = 0; // Do not increment the right array indices
    }else{ // Same shape mvarrays
        strides = left->strides;
        shapes = left->shapes;
        ns = left->ns;
        size = left->strides[0];
    }
    
    // If the algebras are different and at least one is a generated type 
    if(left->GA != right->GA && (left->type->generated || right->type->generated)){
        int isleft_bigger = is_bigger_metric(left->GA,right->GA);
        if(isleft_bigger == -1)
           return NULL;
        if(isleft_bigger){
            cast_left = 0;
            left_cast = left;
            right_cast = cast_mvarray_inherit_type(right,left->type);
            GA = left->GA;
            product = left->type->math_funcs->product;
            out = new_mvarray_inherit_type(left->GA, ns, strides, shapes, left->type);
        }else{
            cast_left = 1;
            left_cast = cast_mvarray_inherit_type(left,right->type);
            right_cast = right;
            GA = right->GA;
            product = right->type->math_funcs->product;
            out = new_mvarray_inherit_type(right->GA, ns, strides, shapes, right->type);
        }
    }else if(left->GA == right->GA){
        left_cast = left;
        right_cast = right;
        GA = left->GA;
        product = left->type->math_funcs->product;
        out = new_mvarray_inherit_type(left->GA, ns, strides, shapes, left->type);
    }else goto failure;

    for(Py_ssize_t i = 0; i < size; i++){
        if(!product(out->data + i*out->type->basic_size,
                    left_cast->data + i*left_cast->type->basic_size*left_inc,
                    right_cast->data + i*right_cast->type->basic_size*right_inc,GA,ptype))
                    goto failure;
    }
    goto success;

failure:
    multivector_array_dealloc(out);
    out = NULL;

success:
    if(cast_left != -1){
        if(cast_left)
            multivector_array_dealloc(left_cast);
        else
            multivector_array_dealloc(right_cast);
    }

    return out;
}

// Multiplication by a scalar
static PyMvObject* multivector_scalar_product(PyMvObject *data, ga_float scalar, ProductType ptype, int scalar_left){
    PyMvObject *out = NULL;
    PyMvObject *scalar_mv = NULL;
    gaprodfunc product = NULL; 
    gascalarfunc scalar_product = NULL;

    out = new_mvarray_inherit_type(data->GA, data->ns, data->strides, data->shapes, data->type);

    if(ptype == ProductType_inner){ // return 0 if inner product with scalar
        return out;
    }else if(ptype == ProductType_regressive){
        product = data->type->math_funcs->product;
        if(!product){
            multivector_array_dealloc(out);
            return NULL;
        }
        int bitmap = 0;
        scalar_mv = new_multivector_inherit_type(data->GA, data->type);
        scalar_mv->type->data_funcs->init(scalar_mv->data,data->GA,&bitmap,&scalar,1);
        if(scalar_left){
            for(Py_ssize_t i = 0; i < data->strides[0]; i++){
                if(!product(out->data + i*out->type->basic_size,scalar_mv->data,data->data + i*data->type->basic_size,data->GA,ptype)){
                    multivector_array_dealloc(out);
                    return NULL;
                }
            }
            
        } else{
            for(Py_ssize_t i = 0; i < data->strides[0]; i++){
                if(!product(out->data + i*out->type->basic_size,data->data + i*data->type->basic_size,scalar_mv->data,data->GA,ptype)){
                    multivector_array_dealloc(out);
                    return NULL;
                }
            }
        }
        multivector_array_dealloc(scalar_mv);
    }else{
        scalar_product = data->type->math_funcs->scalar_product;
        if(!scalar_product){
            multivector_array_dealloc(out);
            return NULL;
        }
        for(Py_ssize_t i = 0; i < data->strides[0]; i++){
            if(!scalar_product(out->data + i*out->type->basic_size,data->data + i*data->type->basic_size,data->GA,scalar)){
                multivector_array_dealloc(out);
                return NULL;
            }
        }
    }
    return out;
}

//Element wise multiplication by a scalar array
static PyMvObject* multivector_scalar_array_product(PyMvObject *left, PyMvObject *right){
    PyMvObject *out = NULL;
    gascalarfunc scalar_product = NULL;
    PyMvObject *data = NULL;
    PyMvObject *scalar = NULL;

    if(!strcmp("scalar",left->type->type_name)){
        scalar = left;
        data = right;
    }else if(!strcmp("scalar",right->type->type_name)){
        scalar = right;
        data = left;
    }else return NULL;

    Py_ssize_t *strides, *shapes, ns,size;
    int scalar_inc = 1, data_inc = 1;
    if(scalar->strides[0] == 1){ // The scalar mvarray is only one mv
        strides = data->strides;
        shapes = data->shapes;
        ns = data->ns;
        size = data->strides[0];
        scalar_inc = 0; // Do not increment the scalar array indices
    }else if(data->strides[0] == 1){ // The non scalar mvarray is only one mv
        strides = scalar->strides;
        shapes = scalar->shapes;
        ns = scalar->ns;
        size = scalar->strides[0];
        data_inc = 0; // Do not increment the non-scalar array indices
    }else{ // Same shape mvarrays
        strides = data->strides;
        shapes = data->shapes;
        ns = data->ns;
        size = data->strides[0];
    }

    out = new_mvarray_inherit_type(data->GA, ns, strides, shapes, data->type);
    
    scalar_product = data->type->math_funcs->scalar_product;
    if(!scalar_product){
        multivector_array_dealloc(out);
        return NULL;
    }
    for(Py_ssize_t i = 0; i < size; i++){
        ScalarMultivector *scalar_mv = INDEX_DATA(scalar, i*scalar_inc);
        if(!scalar_product(INDEX_DATA(out, i),INDEX_DATA(data, i*data_inc),data->GA,*scalar_mv)){
            multivector_array_dealloc(out);
            return NULL;
        }
    }
    
    return out;
}

static int compare_shapes(PyObject *data0, PyObject *data1){
    PyMvObject *mv0 = (PyMvObject*)data0;
    PyMvObject *mv1 = (PyMvObject*)data1;
    if(mv0->ns != mv1->ns) return 0;
    if(mv0->strides[0] != mv1->strides[0]) return 0;
    for(Py_ssize_t i = 0; i < mv0->ns; i++){
        if(mv0->shapes[i] != mv1->shapes[i]) return 0;
    }
    return 1;
}

static int compare_multivector_types(PyMvObject *arg0, PyMvObject *arg1){
    return !strcmp(arg0->type->type_name,arg1->type->type_name);
}



static PyObject *multivector_tprod(PyObject *cls, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"", "","","ptype", NULL};
    PyObject *arg0 = NULL;
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    ProductType ptype = ProductType_geometric;

    char *type = NULL, *basis = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|s", kwlist, &arg0, &arg1, &arg2, &type))
		return NULL;
    
    ptype = string_to_product_type(type);

    if(!PyObject_IsInstance(arg0,cls) && !PyObject_IsInstance(arg1,cls) && !PyObject_IsInstance(arg2,cls)){
        PyErr_SetString(PyExc_ValueError,"Arguments must be multivectors!");
        return NULL;
    }
    if(!compare_multivector_types((PyMvObject*)arg0, (PyMvObject*)arg1) || !compare_multivector_types((PyMvObject*)arg0, (PyMvObject*)arg2)){
        PyErr_SetString(PyExc_TypeError,"Multivector arrays must be of the same type!");
        return NULL;
    }

    if(!compare_shapes(arg0, arg1) || !compare_shapes(arg0, arg2)){
        PyErr_SetString(PyExc_TypeError,"Multivector arrays must have the same shape!");
        return NULL;
    }
    
    PyMvObject* mv0 = (PyMvObject*)arg0;
    PyMvObject* mv1 = (PyMvObject*)arg1;
    PyMvObject* mv2 = (PyMvObject*)arg2;

    gaternaryprodfunc tprod = mv0->type->math_funcs->ternary_product;

    if(!tprod){
        PyErr_SetString(PyExc_TypeError,"Ternary products are not available for this type!");
        return NULL;
    }

    PyMvObject *out = new_mvarray_from_mvarray(mv0);
    for(Py_ssize_t i = 0; i < mv0->strides[0]; i++){
        if(!tprod(INDEX_DATA(out, i),INDEX_DATA(mv0, i),INDEX_DATA(mv1, i),INDEX_DATA(mv2, i),mv0->GA,ptype)){
            multivector_array_dealloc(out);
            PyErr_SetString(PyExc_ValueError,"Error computing the ternary product!");
            return NULL;
        }
    }
    
    return (PyObject*)out;
}

static PyObject *multivector_product(PyObject *left, PyObject *right, ProductType ptype){
    PyMultivectorObject *data0 = NULL;
    PyMvObject *out = NULL;
    ga_float value = 0;
    int scalar_left = -1;
    int isleft_single = -1;
    
    if(get_scalar(right,&value)) // check if right is a scalar
        data0 = (PyMultivectorObject*)left,scalar_left = 0;
    else if(get_scalar(left,&value)) // check if left is a scalar
        data0 = (PyMultivectorObject*)right,scalar_left = 1;

    if(scalar_left != -1){ // One of the arguments is a scalar
        out = multivector_scalar_product(data0,value,ptype,scalar_left);
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Something wrong computing the product with a scalar!");
            return NULL;
        }
        return (PyObject*)out;
    }

    if(!PyObject_TypeCheck(left,Py_TYPE(right))){
        PyErr_SetString(PyExc_TypeError,"Operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(!compare_shapes(left, right)){
        if(((PyMvObject*)left)->strides[0] == 1)
            isleft_single = 1;
        else if(((PyMvObject*)right)->strides[0] == 1)
            isleft_single = 0;
        else{
            PyErr_SetString(PyExc_TypeError,"Multivector arrays must be of the same shape, or one of them is a single multivector");
            return NULL;
        }
    }
    PyMvObject *left_mv = (PyMultivectorObject*)left;
    PyMvObject *right_mv = (PyMultivectorObject*)right;
    if(!strcmp("scalar",left_mv->type->type_name) || !strcmp("scalar",right_mv->type->type_name)){
        out = multivector_scalar_array_product(left_mv,right_mv);
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Error taking the product with a scalar array!");
            return NULL;
        }
        return (PyObject*)out;
    }

    if(left_mv->GA == right_mv->GA || left_mv->type->generated || right_mv->type->generated)
        out = multivector_casttype_product(left_mv, right_mv, ptype,isleft_single);
    else 
        out = multivector_mixed_product(left_mv, right_mv, ptype,isleft_single);
    
    if(!out){
        PyErr_SetString(PyExc_TypeError,"Probably Incompatible Algebras!");
        return NULL;
    }
    return (PyObject*)out;
}

static PyMvObject* multivector_mixed_addsubtract(PyMvObject *left, PyMvObject *right, int sign, int isleft_single){
    int isleft_bigger = is_bigger_metric(left->GA,right->GA);
    PyAlgebraObject *GA = NULL;
    PyMvObject *out = NULL;
    Py_ssize_t left_inc = 1,right_inc = 1;
    Py_ssize_t size = left->strides[0];
    gamixedaddfunc mixed_add = NULL;

    if(isleft_bigger == -1){
        return NULL;
    }else {
        gaiterinitfunc iter_init_left = left->type->data_funcs->iter_init;
        gaiterinitfunc iter_init_right = right->type->data_funcs->iter_init;
        Py_ssize_t *strides, *shapes, ns;
        
        if(isleft_single == 1){
            strides = right->strides;
            shapes = right->shapes;
            ns = right->ns;
            size = right->strides[0];
            left_inc = 0;
        }else if(isleft_single == 0){
            strides = left->strides;
            shapes = left->shapes;
            ns = left->ns;
            size = left->strides[0];
            right_inc = 0;
        }else{
            strides = left->strides;
            shapes = left->shapes;
            ns = left->ns;
            size = left->strides[0];
        }

        if(isleft_bigger){
            out = new_multivector_array(left->GA,"sparse",ns,strides,shapes);
            GA = left->GA;
            mixed_add = left->mixed->add;
        }
        else{
            out = new_multivector_array(right->GA,"sparse",ns,strides,shapes);
            GA = right->GA;
            mixed_add = right->mixed->add;
        }

        if(isleft_single == 1) {
            size = right->strides[0];
            left_inc = 0;
        }
        else if(isleft_single == 0) {
            size = left->strides[0];
            right_inc = 0;
        }
        
        if(mixed_add){
            for(Py_ssize_t i = 0; i < size; i++){
                PyMultivectorIter iter_left = iter_init_left(left->data + i*left->type->basic_size*left_inc,left->type);
                PyMultivectorIter iter_right = iter_init_right(right->data + + i*right->type->basic_size*right_inc,right->type);
                
                if(!mixed_add(out->data + i*out->type->basic_size,&iter_left,&iter_right,GA,sign)){
                    return NULL;
                }
                PyMem_RawFree(iter_left.index);
                PyMem_RawFree(iter_right.index);
            }
            return out;
            
        }else 
            return NULL;
    }
}

static PyMvObject *multivector_casttype_addsubtract(PyMvObject *left, PyMvObject *right, int sign, int isleft_single){
    PyMvObject *left_cast = NULL;
    PyMvObject *right_cast = NULL;
    PyMvObject *out = NULL;
    PyAlgebraObject *GA = NULL;
    gaprodfunc add = NULL;
    int cast_left = -1;
    Py_ssize_t left_inc = 1,right_inc = 1;
    Py_ssize_t size = left->strides[0];

    Py_ssize_t *strides, *shapes, ns;
    if(isleft_single == 1){
        strides = right->strides;
        shapes = right->shapes;
        ns = right->ns;
        size = right->strides[0];
        left_inc = 0;
    }else if(isleft_single == 0){
        strides = left->strides;
        shapes = left->shapes;
        ns = left->ns;
        size = left->strides[0];
        right_inc = 0;
    }else{
        strides = left->strides;
        shapes = left->shapes;
        ns = left->ns;
        size = left->strides[0];
    }
    
    if(left->type->generated || right->type->generated){
        int isleft_bigger = is_bigger_metric(left->GA,right->GA);
        if(isleft_bigger == -1)
           return NULL;
        
        if(isleft_bigger){
            cast_left = 0;
            left_cast = left;
            right_cast = cast_mvarray_inherit_type(right,left->type);
            GA = left->GA;
            add = left->type->math_funcs->add;
            out = new_mvarray_inherit_type(left->GA, ns, strides, shapes, left->type);
        }else{
            cast_left = 1;
            left_cast = cast_mvarray_inherit_type(left,right->type);
            right_cast = right;
            GA = right->GA;
            add = right->type->math_funcs->add;
            out = new_mvarray_inherit_type(right->GA, ns, strides, shapes, right->type);
        }
    }else if(left->GA == right->GA){
        left_cast = left;
        right_cast = right;
        GA = left->GA;
        add = left->type->math_funcs->add;
        out = new_mvarray_inherit_type(left->GA, ns, strides, shapes, left->type);
    }else goto failure;

    if(isleft_single == 1) {
        size = right->strides[0];
        left_inc = 0;
    }
    else if(isleft_single == 0) {
        size = left->strides[0];
        right_inc = 0;
    }

    for(Py_ssize_t i = 0; i < size; i++){
        if(!add(out->data + i*out->type->basic_size,
                left_cast->data + i*left_cast->type->basic_size*left_inc,
                right_cast->data + i*right_cast->type->basic_size*right_inc,GA,sign))
            goto failure;
        
    }

    goto success;

failure:
    multivector_array_dealloc(out);
    out = NULL;

success:
    if(cast_left != -1){
        if(cast_left)
            multivector_array_dealloc(left_cast);
        else
            multivector_array_dealloc(right_cast);
    }

    return out;
}

// Addition by a scalar
static PyMvObject* multivector_scalar_add(PyMvObject *data, ga_float scalar, int sign){
    PyMvObject *out = NULL;
    gascalaraddfunc scalar_add = NULL;

    out = new_mvarray_inherit_type(data->GA, data->ns, data->strides, data->shapes, data->type);
    scalar_add = data->type->math_funcs->scalar_add;
    if(!scalar_add){
        multivector_array_dealloc(out);
        return NULL;
    }
    for(Py_ssize_t i = 0; i < data->strides[0]; i++){
        if(!scalar_add(out->data + i*out->type->basic_size,data->data + i*data->type->basic_size,data->GA,scalar,sign)){
            multivector_array_dealloc(out);
            return NULL;
        }
    }
    
    return out;
}

// Element wise addition by a scalar array
static PyMvObject* multivector_scalar_array_add(PyMvObject *left, PyMvObject *right, int sign){
    PyMvObject *data = NULL;
    PyMvObject *scalar = NULL;
    PyMvObject *out = NULL;
    gascalaraddfunc scalar_add = NULL;
    int sign_scalar = 1;
    int sign_data = 1;

    if(!strcmp("scalar",left->type->type_name)){
        data = (PyMultivectorObject*)right;
        scalar = (PyMultivectorObject*)left;
        sign_data = sign; // Multiply the multivectors by the sign
        
    }else if(!strcmp("scalar",right->type->type_name)){
        data = (PyMultivectorObject*)left;
        scalar = (PyMultivectorObject*)right;
        sign_scalar = sign; // Multiply the scalars by the sign
    }else return NULL;
    
    out = new_mvarray_inherit_type(data->GA, data->ns, data->strides, data->shapes, data->type);
    scalar_add = data->type->math_funcs->scalar_add;
    if(!scalar_add){
        multivector_array_dealloc(out);
        return NULL;
    }
    for(Py_ssize_t i = 0; i < data->strides[0]; i++){
        ScalarMultivector *scalar_mv = INDEX_DATA(scalar, i);
        if(!scalar_add(INDEX_DATA(out, i),INDEX_DATA(data, i),data->GA,sign_scalar*(*scalar_mv),sign_data)){
            multivector_array_dealloc(out);
            return NULL;
        }
    }
    
    return out;
}

static PyObject *multivector_add_subtract(PyObject *left, PyObject *right, int sign){
    PyMultivectorObject *data0 = NULL;
    PyMvObject *out = NULL;
    ga_float value = 0;
    int scalar_left = -1;
    int isleft_single = -1;
    
    if(get_scalar(right,&value)){ // check if right is a scalar
        data0 = (PyMultivectorObject*)left,scalar_left = 1;
        value *= sign; // multiply by sign
        sign = 1;
    }
    else if(get_scalar(left,&value)) // check if left is a scalar
        data0 = (PyMultivectorObject*)right,scalar_left = 0;

    if(scalar_left != -1){ // One of the arguments is a scalar
        out = multivector_scalar_add(data0,value,sign);
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Something wrong computing the sum/subtraction with a scalar!");
            return NULL;
        }
        return (PyObject*)out;
    }

    if(!PyObject_TypeCheck(left,Py_TYPE(right))){
        PyErr_SetString(PyExc_TypeError,"Operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(!compare_shapes(left, right)){
        if(((PyMvObject*)left)->strides[0] == 1)
            isleft_single = 1;
        else if(((PyMvObject*)right)->strides[0] == 1)
            isleft_single = 0;
        else{
            PyErr_SetString(PyExc_TypeError,"Multivector arrays must be of the same shape, or one of them is a single multivector");
            return NULL;
        }
    }
    PyMvObject *left_mv = (PyMultivectorObject*)left;
    PyMvObject *right_mv = (PyMultivectorObject*)right;

    if(!strcmp("scalar",left_mv->type->type_name) || !strcmp("scalar",right_mv->type->type_name)){
        out = multivector_scalar_array_add(left_mv, right_mv,sign);
        if(!out){
           PyErr_SetString(PyExc_TypeError,"Probably Incompatible Algebras!");
            return NULL;
        }
        return (PyObject*)out;
    }


    if(left_mv->GA == right_mv->GA || left_mv->type->generated || right_mv->type->generated)
        out = multivector_casttype_addsubtract(left_mv, right_mv, sign,isleft_single);
    else 
        out = multivector_mixed_addsubtract(left_mv, right_mv, sign,isleft_single);
    
    if(!out){
        PyErr_SetString(PyExc_TypeError,"Probably Incompatible Algebras!");
        return NULL;
    }
    return (PyObject*)out;
}


PyObject *multivector_outer_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_outer);
}

PyObject *multivector_inner_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_inner);
}

PyObject *multivector_geometric_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_geometric);
}

PyObject *multivector_regressive_product(PyObject *left, PyObject *right){
    return multivector_product(left,right,ProductType_regressive);
}

PyObject *multivector_cast(PyMultivectorObject *self, PyObject *args) {
    char *type_name = NULL;
    PyMultivectorObject *out;
    if(!PyArg_ParseTuple(args, "s", &type_name))
        return NULL;

    out = cast_mvarray_to_type(self,type_name);
    if(!out){
        PyErr_SetString(PyExc_ValueError, "Type name probably incorrect!!");
        return NULL;
    }
    return (PyObject*)out;
}

PyMvObject *multivector_scalar_grade_projection(PyMvObject *self){
    Py_ssize_t size = self->strides[0];
    PyMultivectorIter *iters = init_multivector_array_iters(self);
    PyMvObject *out = new_multivector_array(self->GA, "scalar", self->ns, self->strides, self->shapes);
    if(!out) return NULL;
    if(self->type->generated || self->type->ntype == MultivectorType_dense){
        for(Py_ssize_t i = 0; i < size; i++){
            ScalarMultivector *scalar = INDEX_DATA(out, i);
            iters[i].next(&iters[i]);
            if(iters[i].bitmap != 0) return NULL;
            *scalar = iters[i].value;
        }
    }else{
        for(Py_ssize_t i = 0; i < size; i++){
            ScalarMultivector *scalar = INDEX_DATA(out, i);
            *scalar = 0;
            while(iters[i].next(&iters[i]))
                if(iters[i].bitmap == 0)
                    *scalar += iters[i].value;
        }
    }

    free_multivector_iter(iters, size);
    return out;
}

PyObject *multivector_grade_project(PyMultivectorObject *self, PyObject *args, PyObject *Py_UNUSED(ignored)){
    int *grades = NULL;
    Py_ssize_t size = -1;
    PyMultivectorObject *out = NULL;

    if((size = parse_arguments_as_grades(self->GA, args, &grades)) < 0){ 
        PyErr_SetString(PyExc_TypeError, "Invalid grades or arguments are empty!!");
        return NULL;
    }

    if(size == 1 && grades[0] == 0 && self->strides[0] == 1){
        // If we want grade projection to scalars
        PyMultivectorIter *iter = init_multivector_iter(self,1);
        while(iter->next(iter)){
            if(iter->bitmap == 0){
                ga_float value = iter->value;
                PyMem_RawFree(grades);
                free_multivector_iter(iter,1);
                return (PyObject*)PyFloat_FromDouble((double)value);
            }
        }
        PyMem_RawFree(grades);
        free_multivector_iter(iter,1);
        return (PyObject*)PyFloat_FromDouble((double)0.0);

    }else if(size == 1 && grades[0] == 0){
        out = multivector_scalar_grade_projection(self); // convert to scalar type multivector
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Probably invalid scalar type!");
            return NULL;
        }
        PyMem_RawFree(grades);
        return (PyObject*)out;
    }

    gaunarygradefunc grade_project = self->type->math_funcs->grade_project;
    if(grade_project){
        out = new_mvarray_from_mvarray(self);
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Error creating array!");
            return NULL; // raise error
        }

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            if(!grade_project(INDEX_DATA(out,i),INDEX_DATA(self,i),self->GA,grades,size)){
                PyErr_SetString(PyExc_TypeError, "Error projecting multivector array to the specified grades!");
                multivector_array_dealloc(out);
                return NULL; // raise error
            }
        }
        
    }else{
        return NULL; // raise error
    }

    PyMem_RawFree(grades);
    return (PyObject*)out;
}

PyObject *multivector_add(PyObject *left, PyObject *right){
    return (PyObject*)multivector_add_subtract(left,right,1);
}

PyObject *multivector_subtract(PyObject *left, PyObject *right){
    return (PyObject*)multivector_add_subtract(left,right,-1);
}

PyObject *multivector_invert(PyMultivectorObject *self){
    gaunaryfunc reverse = self->type->math_funcs->reverse;
    PyMultivectorObject *out;
    if(reverse){
        out = new_mvarray_from_mvarray(self);
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Error creating array!");
            return NULL; // raise error
        }

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            if(!reverse(INDEX_DATA(out,i),INDEX_DATA(self,i),self->GA)){
                PyErr_SetString(PyExc_TypeError, "Error reversing multivector array!");
                multivector_array_dealloc(out);
                return NULL; // raise error
            }
        }
    }else
        return NULL; // raise error
    
    return (PyObject*)out;
}

PyObject* multivector_dual(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    gaunaryfunc dual = self->type->math_funcs->dual;
    PyMultivectorObject *out = NULL;
    if(dual){
        out = new_mvarray_from_mvarray(self);
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Error creating array!");
            return NULL; // raise error
        }

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            if(!dual(INDEX_DATA(out,i),INDEX_DATA(self,i),self->GA)){
                PyErr_SetString(PyExc_TypeError, "Error dualizing multivector array!");
                multivector_array_dealloc(out);
                return NULL; // raise error
            }
        }
    }else{
        PyErr_SetString(PyExc_TypeError, "OPeration not available for the specified type!");
        return NULL; // raise error
    }
    return (PyObject*)out;
}

PyObject* multivector_undual(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    gaunaryfunc undual = self->type->math_funcs->undual;
    PyMultivectorObject *out = NULL;
    if(undual){
        out = new_mvarray_from_mvarray(self);
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Error creating array!");
            return NULL; // raise error
        }

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            if(!undual(INDEX_DATA(out,i),INDEX_DATA(self,i),self->GA)){
                PyErr_SetString(PyExc_TypeError, "Error undualizing multivector array!");
                multivector_array_dealloc(out);
                return NULL; // raise error
            }
        }
    }else{
        PyErr_SetString(PyExc_TypeError, "Operation not available for the specified type!");
        return NULL; // raise error
    }
    return (PyObject*)out;
}

static PyObject *multivector_sign(PyMultivectorObject *self, ga_float value){
    gascalarfunc scalar_product = self->type->math_funcs->scalar_product;
    PyMultivectorObject *out = NULL;
    if(scalar_product){
        out = new_mvarray_from_mvarray(self);
        if(!out){
            PyErr_SetString(PyExc_TypeError, "Error creating array!");
            return NULL; // raise error
        }

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            if(!scalar_product(INDEX_DATA(out,i),INDEX_DATA(self,i),self->GA,value)){
                PyErr_SetString(PyExc_TypeError, "Error scalar_productizing multivector array!");
                multivector_array_dealloc(out);
                return NULL; // raise error
            }
        }
    }else{
        PyErr_SetString(PyExc_TypeError, "OPeration not available for the specified type!");
        return NULL; // raise error
    }
    return (PyObject*)out;
}

PyObject *multivector_negative(PyMultivectorObject *self){
    return multivector_sign(self,-1);
}

PyObject *multivector_positive(PyMultivectorObject *self){
    return multivector_sign(self,1);
}

int get_biggest_algebra_index(PyObject *cls, PyObject *args){
    PyAlgebraObject *biggest_ga;
    Py_ssize_t index = 0;
    int same_algebra_and_type = 1;
    PyMultivectorObject *data0;
    Py_ssize_t size = PyTuple_Size(args);
    PyObject *arg0 = PyTuple_GetItem(args,0);
    int ntype;
    if(!PyObject_IsInstance(arg0,cls)) return -1;
    data0 = (PyMultivectorObject*)arg0;
    biggest_ga = data0->GA;
    ntype = data0->type->ntype;

    for(Py_ssize_t i = 1; i < size; i++){
        PyObject *argi = PyTuple_GetItem(args,i);
        // check if objects are multivectors
        if(!PyObject_IsInstance(argi,cls)){
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.mvarray");
            return -1;
        }
        data0 = (PyMultivectorObject*)argi;
        // check if object are compatible
        if(biggest_ga != data0->GA){
            int is0_bigger = is_bigger_metric(data0->GA,biggest_ga);
            if(is0_bigger == -1) return -1;
            else if(is0_bigger == 1) biggest_ga = data0->GA, index = i;
            same_algebra_and_type = 0;
        }else if(data0->type->ntype != ntype) same_algebra_and_type = 0;
    }

    if(same_algebra_and_type) return -2;
    return index;
}

static PyObject *multivector_algebra(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    Py_XINCREF((PyObject*)self->GA);
    return (PyObject*)self->GA;
}

static PyObject *multivector_type(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    PyObject *array_repr = Py_BuildValue("s", ".mvarray.");
	PyObject *subtype = Py_BuildValue("s", self->type->type_name);
    PyObject *ga_repr = PyObject_Repr((PyObject*)self->GA);
    
    PyObject *out = PyUnicode_Concat(ga_repr, array_repr);
    out = PyUnicode_Concat(out, subtype);
    
    Py_XDECREF(array_repr);
    Py_XDECREF(subtype);
    Py_XDECREF(ga_repr);
    return out;
}

PyObject *multivector_sum(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    // In the mean time should also add axis as an optional argument (when einsum is implemented)

    gaatomicfunc atomic_add = self->type->math_funcs->atomic_add;
    if(!atomic_add){
        PyErr_SetString(PyExc_NotImplementedError,"The atomic sum operation for these types is not implemented");
        return NULL; 
    }
    PyMvObject *out = new_multivector_inherit_type(self->GA, self->type);
    if(!out){
        PyErr_SetString(PyExc_MemoryError,"Error creating new multivector");
        return NULL; 
    }
    if(!atomic_add(out->data,self->data,self->GA,self->strides[0])){
        multivector_array_dealloc(out);
        PyErr_SetString(PyExc_ValueError,"Error atomic adding multivectors");
        return NULL; 
    }
    return (PyObject*)out;
}

PyObject *multivector_atomic_product(PyMvObject *self, ProductType ptype){
    // In the mean time should also add axis as an optional argument (when einsum is implemented)

    gaatomicprodfunc atomic_product = self->type->math_funcs->atomic_product;
    if(!atomic_product){
        PyErr_SetString(PyExc_NotImplementedError,"The atomic product operation for these types is not implemented");
        return NULL; 
    }
    PyMvObject *out = new_multivector_inherit_type(self->GA, self->type);
    if(!out){
        PyErr_SetString(PyExc_MemoryError,"Error creating new multivector");
        return NULL; 
    }
    if(!atomic_product(out->data,self->data,self->GA,self->strides[0],ptype)){
        multivector_array_dealloc(out);
        PyErr_SetString(PyExc_ValueError,"Error atomic producting multivectors");
        return NULL; 
    }
    return (PyObject*)out;
}

PyObject *multivector_atomic_geometric_product(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    return multivector_atomic_product(self,ProductType_geometric);
}

PyObject *multivector_atomic_outer_product(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    return multivector_atomic_product(self,ProductType_outer);
}

PyObject *multivector_atomic_inner_product(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    return multivector_atomic_product(self,ProductType_inner);
}

PyObject *multivector_atomic_regressive_product(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    return multivector_atomic_product(self,ProductType_regressive);
}

PyObject *multivector_concat(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    PyMultivectorSubType *type = NULL;
    PyMvObject *mv, *out;
     
    if(size != 1){
        PyErr_SetString(PyExc_ValueError,"Number of arguments has to be one!");
        return 0;
    }
    
    // Consider only the first argument
    PyObject *list = PyTuple_GetItem(args,0);
    if(!PyList_Check(list)){
        PyErr_SetString(PyExc_ValueError,"First argument must be a list!");
        return 0;
    }
    size = PyList_Size(list);
    for(Py_ssize_t i = 0; i < size; i++){
        PyObject *argi = PyList_GetItem(list,i);
        if(!PyObject_IsInstance(argi,cls)){
            PyErr_SetString(PyExc_ValueError,"Arguments must be multivectors!");
            return 0;
        }
        mv = (PyMvObject*)argi;
        if(!type) type = mv->type;
        if(type != mv->type){
            // Cast to the type of the biggest algbra
            PyErr_SetString(PyExc_NotImplementedError,"Mixed type concatenation is still not implemented!");
            return 0;
        }
        if(mv->strides[0] != 1){
            // Need to deal with shapes but it is essentially the same
            PyErr_SetString(PyExc_NotImplementedError,"Concatenation of arrays is not implemented!");
            return 0;
        }
    }
    Py_ssize_t strides[2] = {size,1};
    Py_ssize_t shape = size;
    
    out = new_mvarray_inherit_type(mv->GA,1,strides,&shape,mv->type);

    for(Py_ssize_t i = 0; i < size; i++){
        mv = (PyMvObject*)PyList_GetItem(list,i);
        gaiterinitfunc iter_init = mv->type->data_funcs->iter_init;
        PyMultivectorIter iter = iter_init(mv->data, mv->type);
        gacastfunc cast = out->type->data_funcs->cast;
        if(!cast(&iter,INDEX_DATA(out,i),out->GA)){
            PyErr_SetString(PyExc_MemoryError,"Error copying data!");
            multivector_array_dealloc(out);
            return 0;
        }
        PyMem_RawFree(iter.index);
    }
    return (PyObject*)out;
}


// Element wise division by a scalar array
static PyMvObject* multivector_scalar_array_divide(PyMvObject *data, PyMvObject *scalar){
    PyMvObject *out = NULL;
    gascalarfunc scalar_product = NULL;

    out = new_mvarray_inherit_type(data->GA, data->ns, data->strides, data->shapes, data->type);
    
    scalar_product = data->type->math_funcs->scalar_product;
    if(!scalar_product){
        multivector_array_dealloc(out);
        return NULL;
    }
    for(Py_ssize_t i = 0; i < data->strides[0]; i++){
        ScalarMultivector *scalar_mv = INDEX_DATA(scalar, i);
        if(!scalar_product(INDEX_DATA(out, i),INDEX_DATA(data, i),data->GA,1.0/(*scalar_mv))){
            multivector_array_dealloc(out);
            return NULL;
        }
    }
    
    return out;
}


// Division of a scalar array
static PyMvObject* multivector_scalar_divide(PyMvObject *scalar_array, ga_float scalar){
    PyMvObject *out = NULL;

    out = new_mvarray_inherit_type(scalar_array->GA, scalar_array->ns, scalar_array->strides, scalar_array->shapes, scalar_array->type);
    ScalarMultivector *scalar_data = (ScalarMultivector*)scalar_array->data;
    ScalarMultivector *scalar_out = (ScalarMultivector*)out->data;
    
    for(Py_ssize_t i = 0; i < scalar_array->strides[0]; i++)
        scalar_out[i] = scalar/scalar_data[i];
    
    return out;
}

PyObject *multivector_divide(PyObject *left, PyObject *right){
    PyMvObject *data0 = NULL;
    PyMvObject *out = NULL;
    ga_float value = 0;
    int scalar_left = -1;
    
    if(get_scalar(right,&value)) // check if right is a scalar
        data0 = (PyMvObject*)left,scalar_left = 0; // Left is not a scalar
    else if(get_scalar(left,&value)) {
        data0 = (PyMvObject*)right,scalar_left = 1;// Left is a scalar
        // If right is not a scalar
        if(strcmp("scalar",data0->type->type_name)){ // Then right must be a scalar array
            PyErr_SetString(PyExc_NotImplementedError,"Division by a multivector is still not implemented!!");
            return NULL;
        }
    }

    if(scalar_left == 0){ // Right is a scalar, Left is a multivector
        out = multivector_scalar_product(data0,1.0/value,ProductType_geometric,scalar_left);
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Something wrong computing the division with a scalar!");
            return NULL;
        }
    }else if(scalar_left == 1){ // Left is a scalar, Right is a scalar array
        out = multivector_scalar_divide(data0,value);
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Something wrong computing the division with a scalar!");
            return NULL;
        }
    }else {
        PyMvObject *scalar = (PyMvObject*)right;
        if(!strcmp("scalar",scalar->type->type_name)){  // The right multivector is a scalar array, Left is a multivector
            data0 = (PyMvObject*)left;
            out = multivector_scalar_array_divide(data0,scalar);
        }else{
            PyErr_SetString(PyExc_NotImplementedError,"Division by a multivector is still not implemented!!");
            return NULL;
        }
        if(!out){
            PyErr_SetString(PyExc_TypeError,"Something wrong computing the division with a scalar!");
            return NULL;
        }
    }
    return (PyObject*)out;
}


static int check_arguments(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    if(size > 1 || size == 0) {
        PyErr_SetString(PyExc_ValueError,"Number of arguments can only be one!");
        return 0;
    }
    PyObject *arg0 = PyTuple_GetItem(args,0);
    if(!PyObject_IsInstance(arg0,cls) && !PyFloat_Check(arg0) && !PyLong_Check(arg0)){
        PyErr_SetString(PyExc_ValueError,"Argument must be either a multivector, a float or an integer!");
        return 0;
    }
    

    return 1;
}

// Apply an element wise operation to a scalar array
static PyObject* multivector_scalar_array_operation(PyObject *self, scalarop op){
    PyMvObject *scalar_array = (PyMvObject*)self;
    if(PyLong_Check(self))
        return (PyObject*)PyFloat_FromDouble(op(PyLong_AsDouble(self)));
    else if(PyFloat_Check(self))
        return (PyObject*)PyFloat_FromDouble(op(PyFloat_AsDouble(self)));
    

    if(strcmp(scalar_array->type->type_name,"scalar")){
        PyErr_SetString(PyExc_ValueError,"Argument must be a scalar multivector");
        return NULL;
    }

    PyMvObject *out = NULL;

    out = new_mvarray_inherit_type(scalar_array->GA, scalar_array->ns, scalar_array->strides, scalar_array->shapes, scalar_array->type);
    ScalarMultivector *scalar_data = (ScalarMultivector*)scalar_array->data;
    ScalarMultivector *scalar_out = (ScalarMultivector*)out->data;
    
    for(Py_ssize_t i = 0; i < scalar_array->strides[0]; i++)
        scalar_out[i] = op(scalar_data[i]);
    
    return (PyObject*)out;
}

static PyObject *multivector_sqrt(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,sqrt);
}

static PyObject *multivector_cos(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,cos);
}

static PyObject *multivector_sin(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,sin);
}

static ScalarMultivector sign(ScalarMultivector scalar){
  return (scalar > 0) - (scalar < 0);
}

static PyObject *multivector_signum(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,sign);
}

static PyObject *multivector_absolute(PyObject *self){
    return (PyObject*)multivector_scalar_array_operation(self,fabs);
}

static PyObject *multivector_exp(PyObject *cls, PyObject *args){
    
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,exp);
}

static PyObject *multivector_cosh(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,cosh);
}

static PyObject *multivector_sinh(PyObject *cls, PyObject *args){
    if(!check_arguments(cls,args))
        return NULL;
    PyObject *scalar_array = (PyObject*)PyTuple_GetItem(args,0);
    return (PyObject*)multivector_scalar_array_operation(scalar_array,sinh);
}
/*
// This code is for comparing multivectors (concretely scalar multivectors)

 #define LESS_THAN(left,right) left < right

static ScalarMultivector less_than(ScalarMultivector left, ScalarMultivector right){
    return left < right;
}

static ScalarMultivector less_than_equal(ScalarMultivector left, ScalarMultivector right){
    return left <= right;
}

static ScalarMultivector equal(ScalarMultivector left, ScalarMultivector right){
    return left == right;
}

static ScalarMultivector not_equal(ScalarMultivector left, ScalarMultivector right){
    return left != right;
}

static ScalarMultivector greater_than(ScalarMultivector left, ScalarMultivector right){
    return left > right;
}

static ScalarMultivector greater_than_equal(ScalarMultivector left, ScalarMultivector right){
    return left >= right;
}


PyObject *multivector_richcompare(PyMvObject *self, PyObject *other, int op){
    ga_float value = 0;
    int other_isscalar = 1;
    scalarcomp comp = NULL;
    PyMvObject *out = NULL;

    // Not a scalar array
    if(strcmp(self->type->type_name,"scalar")){
        PyErr_SetString(PyExc_ValueError,"Argument must be a scalar multivector array");
        return NULL;
    }

    // Not a float and not an int
    if(!get_scalar(other,&value)){
        other_isscalar = 0;
        if(!PyObject_IsInstance(other,(PyObject*)self)){
            PyErr_SetString(PyExc_ValueError,"Other must be a multivector or a scalar");
            return 0;
        }
    }
    switch(op){
        case Py_LT: // Less than
            comp = less_than;
            break;
        case Py_LE: // Less than or equal
            comp = less_than_equal;
            break;
        case Py_EQ: // Equal
            comp = equal;
            break;
        case Py_NE:
            comp = not_equal;
            break;
        case Py_GT:
            comp = greater_than;
            break;
        case Py_GE:
            comp = greater_than_equal;
            break;
        default:
            return NULL;
    }
    if(other_isscalar){
        out = new_mvarray_inherit_type(self->GA, self->ns, self->strides, self->shapes, self->type);
        ScalarMultivector *scalar_data = (ScalarMultivector*)self->data;
        ScalarMultivector *scalar_out = (ScalarMultivector*)out->data;

        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            scalar_out[i] = comp(scalar_data[i],value);
        }
        return (PyObject*)out;
    } else{
        PyMvObject *other_mv = (PyMvObject*)other;
        if(!compare_shapes((PyObject*)self,(PyObject*)other_mv)){    
            PyErr_SetString(PyExc_TypeError,"Multivector arrays must be of the same shape");
            return NULL;
        }
    
        if(strcmp(other_mv->type->type_name,"scalar")){
            PyErr_SetString(PyExc_ValueError,"Argument must be a scalar multivector array");
            return NULL;
        }

        out = new_mvarray_inherit_type(self->GA, self->ns, self->strides, self->shapes, self->type);
        ScalarMultivector *scalar_left = (ScalarMultivector*)self->data;
        ScalarMultivector *scalar_right = (ScalarMultivector*)other_mv->data;
        ScalarMultivector *scalar_out = (ScalarMultivector*)out->data;
        for(Py_ssize_t i = 0; i < self->strides[0]; i++){
            scalar_out[i] = comp(scalar_left[i],scalar_right[i]);
        }
        return (PyObject*)out;
    }
}

// This code is to implement the mapping protocol

// Given an array of indices computes an iterator that iterates for that constant indices
static PyMultipleArrayIter init_single_array_iter_sequence(PyMvObject *self, Py_ssize_t *index, Py_ssize_t idx_size){
    PyMultipleArrayIter iter;
    Py_ssize_t size = self->ns - idx_size;
    iter.arrays = (PyMvBasicArray*)PyMem_RawMalloc(sizeof(PyMvBasicArray));
    iter.arrays->data = self->data;
    iter.arrays->data0 = self->data;
    // Get the item
    for(Py_ssize_t i = 0; i < idx_size; i++)
        iter.arrays->data0 += index[i]*self->strides[i+1]*self->type->basic_size;
    
    if(size >= 1){
        iter.arrays->strides = (Py_ssize_t*)PyMem_RawMalloc((self->ns + 1 - idx_size)*sizeof(Py_ssize_t));
        for(Py_ssize_t i = idx_size; i < self->ns + 1; i++)
            iter.arrays->strides[i-idx_size] = self->strides[i];
        
        iter.shapes = (Py_ssize_t*)PyMem_RawMalloc((self->ns - idx_size)*sizeof(Py_ssize_t));
        iter.index = (Py_ssize_t*)PyMem_RawMalloc((self->ns-idx_size)*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < self->ns-idx_size; i++){
            iter.index[i] = 0;
            iter.shapes[i] = self->shapes[i+idx_size];
        }
    }else{
        // Do not iterate
        iter.arrays->strides = NULL;
        iter.shapes = NULL;
        iter.index = NULL;
    }
    iter.arrays->ns = self->ns - idx_size;
    iter.arrays->basic_size = self->type->basic_size;
    iter.nm = 1; // The number of arrays
    
    iter.dim = -1;
    iter.ns = self->ns - idx_size;
    
    return iter;
}

# Get item/items from index/indices x[5] or x[1,2,3]
static PyObject* multivector_get_item(PyMvObject *self, PyObject *key){
    Py_ssize_t *index = NULL;
    Py_ssize_t idx_size = 0;
    PyMvObject *out = NULL;


    if(PyLong_Check(key)){
        index = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
        *index = (int)PyLong_AsLong(key);
    }else if(PyList_Check(key)){
        idx_size = PyTuple_GET_SIZE(key);
        index = (Py_ssize_t*)PyMem_RawMalloc(idx_size*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < idx_size; i++){
            PyObject *idx = PyTuple_GET_ITEM(key, i);
            index[i] = (int)PyLong_AsLong(idx);
        }
    }

    PyMultipleArrayIter iter = init_single_array_iter_sequence(self,index,idx_size);
}
*/
/*
PyObject* multivector_atomic_add(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    gaatomicfunc add = NULL;
    gamixedatomicfunc mixed_add = NULL;
    PyMultivectorObject *data0, *data1, *out;
    PyMultivectorObject *data_array;
    Py_ssize_t index;

    if(size <= 1){
        PyErr_SetString(PyExc_ValueError,"number of arguments must be at least two");
        return NULL;
    }
    if((index = get_biggest_algebra_index(cls,args)) == -1)
        return NULL;

    if(index == -2){
        if(size == 2){
            gaaddfunc binary_add;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            binary_add = data0->type->math_funcs->add;
            if(binary_add){
                return (PyObject*)binary_add(data0,data1,1);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The binary sum operation for these types is not implemented");
                return NULL;
            }
        }
    }

    data_array = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++)
        data_array[i] = *((PyMultivectorObject*)PyTuple_GetItem(args,i));

    if(index >= 0){ // dispatch mixed type operations
        mixed_add = data_array[index].mixed->atomic_add;
        if(mixed_add){
            out = mixed_add(data_array,size,&data_array[index]);
        }else{
            PyMem_RawFree(data_array);
            PyErr_SetString(PyExc_NotImplementedError,"The atomic mixed sum operation for these types is not implemented");
            return NULL; // raise not implemented error
        }
    }
    else{
        add = data_array->type->math_funcs->atomic_add;
        if(add){
            out = add(data_array,size);
        }else{
            PyMem_RawFree(data_array);
            PyErr_SetString(PyExc_NotImplementedError,"The atomic sum operation for these types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    PyMem_RawFree(data_array);
    return (PyObject*)out;
}




static PyObject* multivector_atomic_product(PyObject *cls, PyObject *args, ProductType ptype){
    Py_ssize_t size = PyTuple_Size(args);
    gaatomicprodfunc product = NULL;
    gamixedatomicprodfunc mixed_product = NULL;
    PyMultivectorObject *data0 = NULL, *data1 = NULL, *data2 = NULL, *out = NULL;
    PyMultivectorObject *data_array = NULL;
    Py_ssize_t index;

    if(size <= 1){
        PyErr_SetString(PyExc_ValueError,"number of arguments must be at least two");
        return NULL;
    }
    // check if objects are multivectors
    if((index = get_biggest_algebra_index(cls,args)) == -1)
        return NULL;

    if(index == -2){
        if(size == 2){
            gaprodfunc binary_product;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            binary_product = data0->type->math_funcs->product;
            if(binary_product){
                return (PyObject*)binary_product(data0,data1,ptype);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The binary product operation for these types is not implemented");
                return NULL;
            }
        }else if(size == 3){
            gaternaryprodfunc ternary_product;
            data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
            data1 = (PyMultivectorObject*)PyTuple_GetItem(args,1);
            data2 = (PyMultivectorObject*)PyTuple_GetItem(args,2);

            ternary_product = data0->type->math_funcs->ternary_product;
            if(ternary_product){
                return (PyObject*)ternary_product(data0,data1,data2,ptype);
            }else{
                PyErr_SetString(PyExc_NotImplementedError,"The ternary product operation for these types is not implemented");
                return NULL;
            }
        }
    }

    data_array = (PyMultivectorObject*)PyMem_RawMalloc(size*sizeof(PyMultivectorObject));
    for(Py_ssize_t i = 0; i < size; i++)
        data_array[i] = *((PyMultivectorObject*)PyTuple_GetItem(args,i));

    if(index >= 0){ // dispatch mixed type operations
        mixed_product = data_array[index].mixed->atomic_product;
        if(mixed_product){
            out = mixed_product(data_array,size,&data_array[index],ptype);
        }else{
            PyMem_RawFree(data_array);
            PyErr_SetString(PyExc_NotImplementedError,"The atomic mixed sum operation for these types is not implemented");
            return NULL; // raise not implemented error
        }
    }else{
        product = data_array->type->math_funcs->atomic_product;
        if(product){
            out = product(data_array,size,ptype);
        }else{
            PyMem_RawFree(data_array);
            PyErr_SetString(PyExc_NotImplementedError,"The atomic product operation for these types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    PyMem_RawFree(data_array);
    return (PyObject*)out;
}


PyObject* multivector_atomic_geometric_product(PyObject *cls, PyObject *args){
    return multivector_atomic_product(cls,args,ProductType_geometric);
}

PyObject* multivector_atomic_outer_product(PyObject *cls, PyObject *args){
    return multivector_atomic_product(cls,args,ProductType_outer);
}

PyObject* multivector_exponential(PyObject *cls, PyObject *args){
    Py_ssize_t size = PyTuple_Size(args);
    gaunaryfunc exponential;
    PyMultivectorObject *data0 = NULL;

    if(size != 1){
        PyErr_SetString(PyExc_ValueError,"number of arguments must be one");
        return NULL;
    }
    // check if objects are multivectors
    if(get_biggest_algebra_index(cls,args) == -1){
        PyErr_SetString(PyExc_ValueError,"the input must be a multivector");
        return NULL;
    }

    data0 = (PyMultivectorObject*)PyTuple_GetItem(args,0);
    exponential = data0->type->math_funcs->exp;
    if(exponential){
        return (PyObject*)exponential(data0);
    }else{
        PyErr_SetString(PyExc_NotImplementedError,"The exponential operation for this type is not implemented");
        return NULL;
    }
    return NULL;
}
*/
int get_multivector_type_table(PyAlgebraObject *ga, char *name, PyMultivectorSubType **subtype){
    
    if(ga->types == NULL) return 0;

    for(Py_ssize_t i = 0; i < ga->tsize; i++){
        if(!strncmp(name,ga->types[i].type_name,strlen(name))){
            *subtype = &ga->types[i];
            return 1;
        }
    }
    return 0;
}

int check_multivector_mixed_type_table(PyMultivectorObject *mv,char *name){
    char *mname;
    Py_ssize_t i = 0;
    if(!strcmp(*mv->mixed->type_names,"any"))
        return 1;
    // The type_names array is null terminated
    while((mname = mv->mixed->type_names[i])){
        if(!strncmp(mname,name,strlen(name)))
            return 1;
        i++;
    }
    return 0;
}

Py_ssize_t parse_list_as_values(PyObject *values, ga_float **values_float) {
	if (!PyList_Check(values)) {
		PyErr_SetString(PyExc_TypeError, "values must be a list");
		return -1;
	}
	Py_ssize_t size = PyList_Size(values);
	if (size <= 0)
		return -1;
	*values_float = (ga_float *)PyMem_RawMalloc(size * sizeof(ga_float));
	for (Py_ssize_t i = 0; i < size; i++) {
		PyObject *value_i = PyList_GetItem(values, i);
		if (PyFloat_Check(value_i))
			(*values_float)[i] = (ga_float)PyFloat_AsDouble(value_i);
		else if (PyLong_Check(value_i))
			(*values_float)[i] = (ga_float)PyLong_AsLong(value_i);
		else {
			//PyErr_SetString(PyExc_TypeError,
			//								"Elements of the list of values must be float or");
			PyMem_RawFree(*values_float);
			return -1;
		}
	}
	return size;
}

int parse_list_as_basis_grades(PyAlgebraObject ga, int *grades, int **bitmaps, Py_ssize_t gsize) {
	// Given an array of grades computes an array of bitmaps
	
	Py_ssize_t size = ga.asize;
	Py_ssize_t psize = 0;
	Py_ssize_t *grade_bool = get_grade_bool(grades, gsize, MAX_GRADE((&ga)) + 1);
	
	for (Py_ssize_t i = 0; i < size; i++) {
		if (grade_bool[GRADE(i)])
			psize++;
	}

	*bitmaps = (int *)PyMem_RawMalloc(psize * sizeof(int));

	Py_ssize_t j = 0;
		for (Py_ssize_t i = 0; i < size; i++) {
			if (grade_bool[GRADE(i)] && j < psize) {
				(*bitmaps)[j] = i;
				j++;
			} else if (j > psize) {
				break;
			}
		}
	PyMem_RawFree(grade_bool);
	return psize;
}

static int tuple_to_long_array(PyObject *key, Py_ssize_t **index, Py_ssize_t **pos, Py_ssize_t *size){
    Py_ssize_t tsize = PyTuple_Size(key);
    *index = (Py_ssize_t*)PyMem_RawMalloc(tsize*sizeof(Py_ssize_t));
    *pos = (Py_ssize_t*)PyMem_RawMalloc(tsize*sizeof(Py_ssize_t));
    // *size = tsize;

    Py_ssize_t j = 0;
    for(Py_ssize_t i = 0; i < tsize; i++){
        PyObject *item = PyTuple_GetItem(key,i);
        if(PyLong_Check(item)){
            Py_ssize_t indexi = (Py_ssize_t)PyLong_AsLong(item);
            (*index)[j] = indexi;
            (*pos)[j] = i;
            j++;// only increment if it is an integer
        } else if(PySlice_Check(item)){
            Py_ssize_t start, stop, step;
            // Py_ssize_t slice_length, length;
            if(PySlice_Unpack(item, &start, &stop, &step) < 0){
                return 0;
            }
            // printf("slice: %zd,%zd,%zd\n",start,stop - PY_SSIZE_T_MAX,step);
            if(start == 0 && stop == PY_SSIZE_T_MAX && step == 1){
                // Check if slice is ':'
                // Do nothing in this case
            }else{
                // if it is a different type of slice
            }
            // slice_length = PySlice_AdjustIndices(length, &start, &stop, step);
            // return 0;
        }
    }
    *size = j;
    return 1;
}

PyObject *PyMultivector_getitem(PyMvObject *self, PyObject *key){
    PyMvObject *new = NULL;
    PyMultipleArrayIter iter;
    Py_ssize_t *index = NULL, *pos = NULL, size = 0;
    if(self->ns == 0){
        PyErr_SetString(PyExc_IndexError, "A single multivector is not indexable.");
        return NULL;
    }
    if(PyLong_Check(key)) {
        Py_ssize_t index = (Py_ssize_t)PyLong_AsLong(key);
        Py_ssize_t pos = 0;
        iter = init_arrays_iter_index(self,&pos,&index,1,&new);
        
    } else if(PyTuple_Check(key)){
        tuple_to_long_array(key,&index,&pos,&size);
        iter = init_arrays_iter_index(self,pos,index,size,&new);
    }

    if(iter.arrays == NULL){
        PyErr_SetString(PyExc_IndexError,"Index out of range.");
        PyMem_RawFree(index);
        PyMem_RawFree(pos);
        return NULL;
    }

    gacastfunc cast = self->type->data_funcs->cast;
    gaiterinitfunc iter_init = self->type->data_funcs->iter_init;

    do{
        // Initialize iterator for the multivector
        PyMultivectorIter mviter = iter_init(iter.arrays[0].data,self->type);
        // copy the multivector data using the cast function
        if(!cast(&mviter,iter.arrays[1].data,new->GA)){
            free_arrays_iter_index(iter);
            PyMem_RawFree(mviter.index);
            return NULL;
        } PyMem_RawFree(mviter.index);
        
    } while(multiple_arrays_iter_next(&iter));

    free_arrays_iter_index(iter);
    PyMem_RawFree(index);
    PyMem_RawFree(pos);

    return (PyObject*)new;
}

PyObject *multivector_shape(PyMvObject *self, PyObject *Py_UNUSED(ignored)){
    PyObject *shape = PyTuple_New(self->ns);
    for(Py_ssize_t i = 0; i < self->ns; i++)
        PyTuple_SetItem(shape, i, (PyObject*)PyLong_FromLong(self->shapes[i]));
    
    return shape;
}

Py_ssize_t PyMultivector_length(PyMvObject *self){
    if(!self->ns)
        return 0;
    
    return self->shapes[0];
}

// Define mapping protocol methods
static PyMappingMethods PyMultivectorMappingMethods = {
    (lenfunc)PyMultivector_length,              /* mp_length */
    (binaryfunc)PyMultivector_getitem, /* mp_subscript */
    // (objobjargproc)PyMultivector_setitem,        /* mp_ass_subscript */
    (objobjargproc)NULL,        /* mp_ass_subscript */
};

static PyNumberMethods PyMultivectorNumberMethods = {
		.nb_multiply = (binaryfunc)multivector_geometric_product,
		.nb_xor = (binaryfunc)multivector_outer_product,
		.nb_and = (binaryfunc)multivector_regressive_product,
		.nb_or = (binaryfunc)multivector_inner_product,
		.nb_add = (binaryfunc)multivector_add,
		.nb_subtract = (binaryfunc)multivector_subtract,
		.nb_invert = (unaryfunc)multivector_invert,
		.nb_negative = (unaryfunc)multivector_negative,
		.nb_positive = (unaryfunc)multivector_positive,
        .nb_true_divide = (binaryfunc)multivector_divide,
        .nb_absolute = (unaryfunc)multivector_absolute,
};


// PyDoc_STRVAR(add_doc, "adds a bunch of multivectors.");
PyDoc_STRVAR(dual_doc, "dualizes the multivector.");
PyDoc_STRVAR(undual_doc, "undualizes the multivector.");
// PyDoc_STRVAR(product_doc, "multiplies a bunch of multivectors.");
// PyDoc_STRVAR(exponential_doc, "takes the exponential of multivectors.");
PyDoc_STRVAR(list_doc,"Returns a list with each coefficient of the multivector. The basis blades are multivectors");
PyDoc_STRVAR(listasbitmap_doc,"Returns a list with each coefficient of the multivector. The basis blades are represented as bitmaps.");
PyDoc_STRVAR(cast_doc, "Casts the multivector to the specified type");
PyDoc_STRVAR(grade_doc, "Returns the grade of a multivector");
PyDoc_STRVAR(algebra_doc, "Returns the algebra of a multivector");
PyDoc_STRVAR(type_doc, "Returns a string indicating the type of the multivector");
PyDoc_STRVAR(shape_doc, "Returns a tuple indicating the shape of the multivector array.");
PyDoc_STRVAR(sum_doc, "Sums all of the multivectors in the array together");
PyDoc_STRVAR(prod_doc, "Geometric multiplies all of the multivectors in the array together, \
                        \nthe first element on the array is the leftmost element to be geometric multipled");
PyDoc_STRVAR(outer_prod_doc, "Outer multiplies all of the multivectors in the array together, \
                              \nthe first element on the array is the leftmost element to be geometric multipled");
PyDoc_STRVAR(inner_prod_doc, "Inner multiplies all of the multivectors in the array together, \
                              \nthe first element on the array is the leftmost element to be geometric multipled.\
                              \nCarefull with precedence, leftmost products are computed first!!");
PyDoc_STRVAR(regressive_prod_doc, "Regressive multiplies all of the multivectors in the array together, \
                                   \nthe first element on the array is the leftmost element to be geometric multipled.\
                                   \nCarefull with precedence, leftmost products are computed first!!");



PyDoc_STRVAR(cosh_doc, "Hyperbolic Cosine of multivectors.");
PyDoc_STRVAR(sinh_doc, "Hyberbolic Sine of multivectors.");
PyDoc_STRVAR(exp_doc, "Exponential of multivectors.");

PyDoc_STRVAR(cos_doc, "Element wise cosine of scalar multivectors.");
PyDoc_STRVAR(sin_doc, "Element wise sine of scalar multivectors.");
PyDoc_STRVAR(sqrt_doc, "Element wise square root of scalar multivectors.");
PyDoc_STRVAR(sign_doc, "Element wise sign of scalar multivectors.");
PyDoc_STRVAR(concat_doc, "Concatenates a list of multivectors.");
PyDoc_STRVAR(tprod_doc, "Computes the product of three multivectors together.");

PyMethodDef multivector_methods[] = {
        {"GA",(PyCFunction)multivector_algebra,METH_NOARGS,algebra_doc},
        {"type",(PyCFunction)multivector_type,METH_NOARGS,type_doc},
        {"shape",(PyCFunction)multivector_shape,METH_NOARGS,shape_doc},
		{"dual", (PyCFunction)multivector_dual, METH_NOARGS, dual_doc},
		{"undual", (PyCFunction)multivector_undual, METH_NOARGS, undual_doc},
        {"sum",(PyCFunction)multivector_sum, METH_NOARGS, sum_doc},
        {"prod",(PyCFunction)multivector_atomic_geometric_product, METH_NOARGS, prod_doc},
        {"outer_prod",(PyCFunction)multivector_atomic_outer_product, METH_NOARGS, outer_prod_doc},
        {"inner_prod",(PyCFunction)multivector_atomic_inner_product, METH_NOARGS, inner_prod_doc},
        {"regressive_prod",(PyCFunction)multivector_atomic_regressive_product, METH_NOARGS, regressive_prod_doc},
		//{"add", (PyCFunction)multivector_atomic_add, METH_VARARGS | METH_CLASS, add_doc},
		//{"geometric_product", (PyCFunction)multivector_atomic_geometric_product,METH_VARARGS | METH_CLASS, product_doc},
		//{"outer_product", (PyCFunction)multivector_atomic_outer_product,METH_VARARGS | METH_CLASS, product_doc},
		//{"exp", (PyCFunction)multivector_exponential, METH_VARARGS | METH_CLASS,exponential_doc},
		{"tolist", (PyCFunction)multivector_list, METH_VARARGS, list_doc},
        {"tolist_as_bitmap",(PyCFunction)multivector_list_as_bitmap, METH_VARARGS, listasbitmap_doc},
		{"cast", (PyCFunction)multivector_cast, METH_VARARGS, cast_doc},
		{"grade", (PyCFunction)multivector_grade, METH_NOARGS, grade_doc},
        {"cos",   (PyCFunction)multivector_cos, METH_VARARGS | METH_CLASS, cos_doc},
        {"sin",   (PyCFunction)multivector_sin, METH_VARARGS | METH_CLASS, sin_doc},
        {"cosh",   (PyCFunction)multivector_cosh, METH_VARARGS | METH_CLASS, cosh_doc},
        {"sinh",   (PyCFunction)multivector_sinh, METH_VARARGS | METH_CLASS, sinh_doc},
        {"exp",   (PyCFunction)multivector_exp, METH_VARARGS | METH_CLASS, exp_doc},
        {"sqrt",   (PyCFunction)multivector_sqrt, METH_VARARGS | METH_CLASS, sqrt_doc},
        {"sign",   (PyCFunction)multivector_signum, METH_VARARGS | METH_CLASS, sign_doc},
        {"concat", (PyCFunction)multivector_concat, METH_VARARGS | METH_CLASS, concat_doc},
        {"tprod", (PyCFunction)multivector_tprod, METH_VARARGS | METH_CLASS | METH_KEYWORDS, tprod_doc},
		{NULL},
};

PyTypeObject PyMultivectorType = {
		PyVarObject_HEAD_INIT(NULL, 0).tp_name = "gasparse.mvarray",
		.tp_doc = PyDoc_STR(
				"Builds a multivector in different types (sparse,dense,blades)"),
		.tp_basicsize = sizeof(PyMultivectorObject),
		.tp_itemsize = 0,
        // .tp_richcompare = (richcmpfunc)multivector_richcompare,
		.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
		.tp_dealloc = (destructor)multivector_array_dealloc,
		.tp_repr = (reprfunc)multivector_repr,
		.tp_str = (reprfunc)multivector_repr,
		.tp_call = (ternaryfunc)multivector_grade_project,
		.tp_new = NULL,
		.tp_as_number = &PyMultivectorNumberMethods,
		.tp_methods = multivector_methods,
        .tp_as_mapping = &PyMultivectorMappingMethods,
};


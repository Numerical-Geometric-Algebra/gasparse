#include "multivector_object.h"
#include "pyerrors.h"
#include "pyport.h"
#include "types.h"

// Multivector Array Initializers
static PyMvObject *init_multivector_array(PyAlgebraObject *GA, Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes){
	if(!GA) return NULL; 
    PyMvObject *array_obj = (PyMvObject*)PyMem_RawMalloc(sizeof(PyMvObject));
    if(!array_obj) return NULL;

    array_obj->shapes = (Py_ssize_t*)PyMem_RawMalloc(ndims*sizeof(Py_ssize_t));
    array_obj->strides = (Py_ssize_t*)PyMem_RawMalloc((ndims+1)*sizeof(Py_ssize_t));
    
    // Copy strides and copy shapes
    for(Py_ssize_t i =0; i < ndims + 1; i++){
        if(i < ndims)
            array_obj->shapes[i] = shapes[i];
        array_obj->strides[i] = strides[i];
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

static int alloc_multivector_data(PyMultivectorObject *obj){
    // Alloc memory for the data
    if(obj->strides != NULL){
        obj->data = (void*)PyMem_RawMalloc(obj->strides[0]*obj->type->basic_size);
        if(!obj->data) return 0;
    }
    return 1;
}

PyMultivectorObject *new_multivector_array(PyAlgebraObject *GA, char *type,  Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes){
    PyMultivectorObject *self = init_multivector_array(GA, ndims, strides, shapes);
    if(!self) return NULL;
    if(!get_multivector_type_table(GA, type, &self->type)) return NULL;
    if(!alloc_multivector_data(self)) return NULL;
    return self;
}

PyMultivectorObject *new_mvarray_inherit_type(PyAlgebraObject *GA,  Py_ssize_t ndims, Py_ssize_t *strides, Py_ssize_t *shapes, PyMultivectorSubType *type){
    PyMultivectorObject *self = init_multivector_array(GA, ndims, strides, shapes);
    if(!self) return NULL;
    self->type = type;
    if(!alloc_multivector_data(self)) return NULL;
    return self;
}

// Single multivector initializers and destructors
PyMultivectorObject *init_multivector(PyAlgebraObject *GA){
    // Initializes a single multivector, also allocs memory for the data
    Py_ssize_t* strides = (Py_ssize_t*)PyMem_RawMalloc(sizeof(Py_ssize_t));
    *(strides) = 1;
    PyMultivectorObject *self = init_multivector_array(GA,0,strides,NULL);
    if(!self) return NULL;
    return self;
}

PyMultivectorObject *new_multivector(PyAlgebraObject *GA, char *type){
    PyMultivectorObject *self = init_multivector(GA);
    if(!self) return NULL;
    if(!get_multivector_type_table(GA, type, &self->type)) return NULL;
    if(!alloc_multivector_data(self)) return NULL;
    return self;
}

PyMultivectorObject *new_multivector_inherit_type(PyAlgebraObject *GA, PyMultivectorSubType *type){
    PyMultivectorObject *self = init_multivector(GA);
    if(!self || !type) return NULL;
    self->type = type;
    // Can only allocate memory after type is set
    if(!alloc_multivector_data(self)) return NULL;
    return self;
}

static void multivector_array_dealloc(PyMvObject *self){
	void *data = self->data;
    gafreefunc free_type = self->type->data_funcs->free;
    if(free_type)
        for(Py_ssize_t i = 0; i < self->strides[0]; i++)
            free_type(self->data + i*self->type->basic_size);

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

int cast_multivectors(PyMvObject *from, PyMvObject *to){
    PyMultivectorIter *iters = init_multivector_array_iters(from);
    gacastfunc cast = to->type->data_funcs->cast;
    for(Py_ssize_t i = 0; i < from->strides[0]; i++){
        if(!cast(iters+i,to->data,to->GA)){
            free_multivector_array_iter(iters,from->strides[0]);
            return 0;
        }
    }
    return 1;
}

PyMvObject *cast_multivector_inherit_type(PyMvObject *from, PyMultivectorSubType *type){
    PyMvObject *to = new_mvarray_inherit_type(from->GA,from->ns,from->strides,from->shapes,type);
    if(!cast_multivectors(from,to)){
        multivector_array_dealloc(to);
        return NULL;
    }
    return to;
}

PyMvObject *cast_multivector_to_type(PyMvObject *from, char *type){
    PyMvObject *to = new_multivector_array(from->GA,type,from->ns,from->strides,from->shapes);
    if(!cast_multivectors(from,to)){
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
        if(size <= 0 && size != nbasis) return -1;
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

static int multivector_array_iter_next(PyMvArrayIter *iter, Py_ssize_t dim){
    if(dim < 0) return 1; // Ignore iteration
    if(dim >= iter->ns) return 0;
    
    iter->data += iter->strides[dim+1]*iter->basic_size;
    iter->index[dim]++;
    if(iter->index[dim] >= iter->shapes[dim]){
        iter->data -= iter->shapes[dim]*iter->strides[dim+1]*iter->basic_size;
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

static PyMultipleArrayIter init_artificial_mtp_arrays_iter(PyMvObject *self){
    PyMultipleArrayIter iter;
    iter.array_iter = (PyMvArrayIter*)PyMem_RawMalloc(sizeof(PyMvArrayIter));
    iter.array_iter->index = (Py_ssize_t*)PyMem_RawMalloc(self->ns*sizeof(Py_ssize_t));
    iter.array_iter->data = self->data;
    iter.array_iter->strides = self->strides;
    iter.array_iter->shapes = self->shapes;
    iter.array_iter->next = multivector_array_iter_next;
    iter.array_iter->ns = self->ns;
    iter.array_iter->basic_size = self->type->basic_size;
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
        PyErr_SetString(PyExc_ValueError,"Error iterating nested lists");
        return NULL;
    }

    Py_ssize_t *shapes_ = NULL;
    if(ndims > 1){
        shapes_ =  PyMem_RawMalloc((ndims-1)*sizeof(Py_ssize_t));
        for(Py_ssize_t i = 0; i < ndims-1; i++) // Discard the innermost shape
            shapes_[i] = shapes[i];
    }
    PyMvObject *mv_array = new_multivector_array(self,type_name,ndims-1,strides,shapes_);
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
    
    for(Py_ssize_t i = 0; i < strides[0]; i++){
        if(!init(mv_array->data + i*mv_array->type->basic_size,self,bitmaps_int,values_float_array[i],bsize)){
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


int is_bigger_metric(PyAlgebraObject *ga0, PyAlgebraObject *ga1){
    Py_ssize_t size = METRIC_SIZE(ga0) < METRIC_SIZE(ga1) ?  METRIC_SIZE(ga0) :  METRIC_SIZE(ga1);
    for(Py_ssize_t i = 0; i < size; i++)
        if(ga0->metric[i] != ga1->metric[i])
            return -1;
    return METRIC_SIZE(ga0) > METRIC_SIZE(ga1);
}

char *type_iter_repr(PyMultivectorIter *iter, PrintTypeMV ptype, Py_ssize_t dsize){
    char *out_str;
    if(ptype == PrintTypeMV_reduced){
        if(dsize){
            char **str_blade = (char**)PyMem_RawMalloc(dsize*sizeof(char*));
            Py_ssize_t len = 0;
            char sep[] = " + ";

            Py_ssize_t i = 0;
            while(iter->next(iter)){
                // only skip small values if its not sparse or blades type
                if(iter->type != MultivectorType_sparse && iter->type != MultivectorType_blades)
                    // maybe should add this as an option when creating the algebra
                    if(ABS(iter->value) < 1e-12) continue; // don't print small values

                char *value = PyOS_double_to_string((double)iter->value,'g',16,0,NULL);
                if(iter->bitmap){
                    char *bitmap = bitmap_to_string(iter->bitmap);
                    Py_ssize_t size = strlen(bitmap) + strlen(value)+3;
                    str_blade[i] = (char*)PyMem_RawMalloc(size*sizeof(char));
                    PyOS_snprintf(str_blade[i],size,"%s^%s",value,bitmap);
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
    }else if(ptype == PrintTypeMV_normal){
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
        };
    }
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
*/
PyObject *multivector_repr(PyMvObject *self){
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
            PyMultivectorIter iter = iter_init(arr_iter.array_iter->data,self->type);
            mv_str = type_iter_repr(&iter,ptype,iter.niters);
            
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
        PyMultivectorIter iter = iter_init(self->data,self->type);
        mv_str = type_iter_repr(&iter,ptype,iter.niters);
        strcpy(out_str,mv_str);
        PyMem_RawFree(mv_str);
        PyMem_RawFree(iter.index);
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

PyObject* multivector_list(PyMvObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"grades","bitmap",NULL};
    PyObject *grades = NULL;
    int *grades_int = NULL;
    PyObject *list = NULL;
    PyObject *bitmap = NULL;
    PyMultivectorObject *dense;
    PyMultivectorIter iter;
	Py_ssize_t *grade_bool = NULL;
	Py_ssize_t size;
    int free_mv = 1;

    int as_bitmap = 0;
    
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "|Op", kwlist, &grades,&as_bitmap))
        return NULL;

    

    // Cast to dense if the multivector is not dense and is not a generated type
    if(!self->type->generated && strcmp(self->type->type_name,"dense")){

        dense = new_multivector(self->GA,"dense");

        if(!dense){
            PyErr_SetString(PyExc_TypeError,"Error populating types table");
            return NULL;
        }
        gacastfunc cast = dense->type->data_funcs->cast;
        if(!cast){ // Check if cast  function is available
            // Free the data0 multivector
            multivector_array_dealloc(dense);
            PyErr_SetString(PyExc_TypeError,"cast function not available for this type");
            return NULL;
        }
        
        
        if(!cast(self->data,dense->data,dense->GA)){
            // Free the dense multivector
            multivector_array_dealloc(dense);
            PyErr_SetString(PyExc_TypeError,"Error casting the multivector");
            return NULL;
        }

    }else {
        free_mv = 0;
        dense = self;
    }
    
    if(grades){
        size = parse_list_as_grades(self->GA, grades,&grades_int);
        if(size <= 0){
            PyErr_SetString(PyExc_TypeError, "Error parsing grades");
            return NULL;
        }

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
        grade_bool  = (Py_ssize_t*)PyMem_RawMalloc((MAX_GRADE(self->GA) +1)*sizeof(Py_ssize_t));
		for(Py_ssize_t i = 0; i < MAX_GRADE(self->GA) +1; i++) 
			grade_bool[i] = 1;
		size = self->GA->asize;
    }

	list = PyList_New(size);
    bitmap = PyList_New(size);
    iter = dense->type->data_funcs->iter_init(dense->data,dense->type);
	
	Py_ssize_t j = 0;
	ga_float basis_value = 1;
    while(iter.next(&iter)){
        if(grade_bool[GRADE(iter.bitmap)] && j < size){
			PyObject *value = PyFloat_FromDouble(iter.value);
			PyList_SetItem(list,j,value);
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
                mv->data = (void*)PyMem_RawMalloc(mv->type->basic_size);
				init(mv->data,self->GA,&iter.bitmap,&basis_value,1);
				PyList_SetItem(bitmap,j,(PyObject*)mv);
			}
			j++;
		}else if (j > size) {
        	break;
      	}
    }
	
    PyMem_RawFree(iter.index);

    PyObject *tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple,0,list);
    PyTuple_SetItem(tuple,1,bitmap);
	PyMem_RawFree(grade_bool);
    if(free_mv)
        Py_XDECREF((PyObject*)dense);

    return tuple;
}

PyObject* multivector_grade(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    int grade = -1;
    PyMultivectorIter iter = self->type->data_funcs->iter_init(self->data,self->type);

    while(iter.next(&iter)){

        if(grade == -1){ //First iteration
            if(iter.value != 0.0)
                grade = iter.grade;
        }else if(grade != iter.grade){
            if(iter.value != 0.0){
                PyMem_RawFree(iter.index);
                Py_RETURN_NONE;
            }
        }
    }if(grade == -1) grade = 0;

    PyMem_RawFree(iter.index);
    return PyLong_FromLong(grade);
}

static PyObject *multivector_product(PyObject *left, PyObject *right, ProductType ptype){
    PyMultivectorObject *data0 = NULL, *data1 = NULL, *def = NULL, *out = NULL;
    ga_float value = 0;
    int is_left = -1;
    gaprodfunc product;
    gamixedprodfunc mixed_product;
    gascalarfunc scalar_product;

    if(get_scalar(right,&value)) // check if right is a scalar
        data0 = (PyMultivectorObject*)left,is_left=1;
    else if(get_scalar(left,&value)) // check if left is a scalar
        data0 = (PyMultivectorObject*)right,is_left=0;

    // One of the arguments is scalar apply multiplication by scalar
    if(data0){
        // return 0 if inner product with scalar
        if(ptype == ProductType_inner){
            out = new_multivector_inherit_type(data0->GA,data0->type);
            data0->type->data_funcs->init(out->data,data0->GA,NULL,NULL,0); // initialize empty multivector
            return (PyObject*)out;
        }else if(ptype == ProductType_regressive){
            // convert value to multivector and then apply the product
            ga_float *pvalue = (ga_float*)PyMem_RawMalloc(sizeof(ga_float));
            int *pbitmap = (int*)PyMem_RawMalloc(sizeof(int));
            *pvalue = value; *pbitmap = 0;
            data1 = new_multivector_inherit_type(data0->GA,data0->type);
            data1->data = (void*)PyMem_RawMalloc(data1->type->basic_size);
            data0->type->data_funcs->init(data1->data,data0->GA,pbitmap,pvalue,1);
            product = data0->type->math_funcs->product;
            if(product){
                if(is_left){
                    if(!product(out->data,data0->data,data1->data,data0->GA,ptype))
                        return NULL;
                    else
                        return (PyObject*)out;
                } else
                    if(!product(out->data,data0->data,data1->data,data0->GA,ptype))
                        return NULL;
                    else
                        return (PyObject*)out;
            }else{
                return NULL;
            }
            Py_XDECREF((PyObject*)data1);
            PyMem_RawFree(pvalue);
            PyMem_RawFree(pbitmap);
            return (PyObject*)out;
        }
        // multiply by scalar
        scalar_product = data0->type->math_funcs->scalar_product;
        if(scalar_product){
            // Allocate single multivector and data
            PyMvObject *out = new_multivector_inherit_type(data0->GA, data0->type);
            if(!scalar_product(out->data,data0->data,data0->GA,value))
                return NULL;
        }else{
            PyErr_SetString(PyExc_NotImplementedError,"The scalar product for this types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(PyObject_TypeCheck(left,Py_TYPE(right))){
        data0 = (PyMultivectorObject*)left;
        data1 = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(data0->GA != data1->GA){
        int is0_bigger;// METRIC_SIZE(data0->GA) > METRIC_SIZE(data1->GA)
        if((is0_bigger = is_bigger_metric(data0->GA,data1->GA)) == -1){
            PyErr_SetString(PyExc_TypeError,"operands must have overlaping metric");
            return NULL;
        }
        if(is0_bigger) mixed_product = data0->mixed->product,def = data0; // data0's GA is bigger
        else           mixed_product = data1->mixed->product,def = data1; // data1's GA is bigger
        if(mixed_product){
            return (PyObject*)mixed_product(data0,data1,def,ptype);
        }else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(data0->type->ntype == data1->type->ntype){
        product = data0->type->math_funcs->product;
        if(product){
            PyMvObject *out = new_multivector_inherit_type(data0->GA, data0->type);
            if(!product(out->data,data0->data,data1->data,out->GA,ptype))
                return NULL;
            else
                return (PyObject*)out;
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
            return NULL; // raise not implemented error
        }
    } else{
        mixed_product = data0->mixed->product;
        if(mixed_product){
            return (PyObject*)mixed_product(data0,data1,data0,ptype);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }

    }

    return NULL;
}

static PyObject *multivector_add_subtract(PyObject *left, PyObject *right, int sign){
    PyMultivectorObject *data0 = NULL, *data1 = NULL, *def = NULL;
    ga_float value = 0;
    gaaddfunc add;
    gamixedaddfunc mixed_add;
    gascalaraddfunc scalar_add;

    if(get_scalar(right,&value)){ // check if right is a scalar
        data0 = (PyMultivectorObject*)left;
        value *= sign; // multiply by sign
        sign = 1;
    }else if(get_scalar(left,&value)){ // check if left is a scalar
        data0 = (PyMultivectorObject*)right;
    }

    if(data0){
        // add a scalar
        scalar_add = data0->type->math_funcs->scalar_add;
        if(scalar_add){
            PyMvObject *out = new_multivector_inherit_type(data0->GA, data0->type);
            if(!scalar_add(out->data,data0->data,out->GA,value,sign))
                return NULL;
        }else{
            PyErr_SetString(PyExc_NotImplementedError,"The scalar product for this types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(PyObject_TypeCheck(left,Py_TYPE(right))){
        data0 = (PyMultivectorObject*)left;
        data1 = (PyMultivectorObject*)right;
    }else{
        PyErr_SetString(PyExc_TypeError,"operands must be of the same type or int or ga_float");
        return NULL;
    }
    if(data0->GA != data1->GA){
        int is0_bigger;// METRIC_SIZE(data0->GA) > METRIC_SIZE(data1->GA)
        if((is0_bigger = is_bigger_metric(data0->GA,data1->GA)) == -1){
            PyErr_SetString(PyExc_TypeError,"operands must have overlapping metric");
            return NULL;
        }
        if(is0_bigger) mixed_add = data0->mixed->add, def = data0; // data0's GA is bigger
        else           mixed_add = data1->mixed->add, def = data1; // data1's GA is bigger
        if(mixed_add){
            return (PyObject*)mixed_add(data0,data1,def,sign);
        }else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    if(data0->type->ntype == data1->type->ntype){
        add = data0->type->math_funcs->add;
        if(add){
            PyMvObject *out = new_multivector_inherit_type(data0->GA, data0->type);
            if(!add(out->data,data0->data,data1->data,out->GA,sign))
                return NULL;
            else
                return (PyObject*)out;
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for these types is not implemented");
            return NULL; // raise not implemented error
        }
    } else{
        mixed_add = data0->mixed->add;
        if(mixed_add){
            return (PyObject*)mixed_add(data0,data1,data0,sign);
        } else {
            PyErr_SetString(PyExc_NotImplementedError,"The product for mixed types is not implemented");
            return NULL; // raise not implemented error
        }
    }

    return NULL;
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
    PyMultivectorObject *data0;
    if(!PyArg_ParseTuple(args, "s", &type_name))
        return NULL;

    data0 = new_multivector(self->GA,type_name);
    if(!data0){
        PyErr_SetString(PyExc_TypeError,"Error populating types table");
        return NULL;
    }
    gacastfunc cast = data0->type->data_funcs->cast;
    if(!cast){
        // Free the data0 multivector
        multivector_array_dealloc(data0);
        PyErr_SetString(PyExc_TypeError,"cast function not available for this type");
        return NULL;
    }
    
    if(!cast(self,data0)){
        // Free the data0 multivector
        multivector_array_dealloc(data0);
        PyErr_SetString(PyExc_TypeError,"Error casting the multivector");
        return NULL;
    }
    // Free the data0 multivector
    return (PyObject*)data0;
}



PyObject *multivector_grade_project(PyMultivectorObject *self, PyObject *args, PyObject *kwds){
    static char *kwlist[] = {"grades",NULL};
    int *grades = NULL;
    PyObject *grades_obj = NULL;
    Py_ssize_t size = -1;
    PyMultivectorObject *out = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &grades_obj))
        return NULL; // raise error

    size = parse_list_as_grades(self->GA,grades_obj,&grades);
    if(size <= 0) return NULL;

    if(size == 1 && grades[0] == 0){
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
    }

    gaunarygradefunc grade_project = self->type->math_funcs->grade_project;
    if(grade_project){
        out = grade_project(self,grades,size);
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
        out = reverse(self);
    }else{
        return NULL; // raise error
    }
    return (PyObject*)out;
}

PyObject* multivector_dual(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    gaunaryfunc dual = self->type->math_funcs->dual;
    if(dual){
        return (PyObject*)dual(self);
    }else {
        return NULL; // raise error
    }
    return NULL;
}

PyObject* multivector_undual(PyMultivectorObject *self, PyObject *Py_UNUSED(ignored)){
    gaunaryfunc undual = self->type->math_funcs->undual;
    if(undual){
        return (PyObject*)undual(self);
    }else {
        return NULL; // raise error
    }
    return NULL;
}

static PyObject *multivector_sign(PyMultivectorObject *self, ga_float value){
    gascalarfunc scalar_product = self->type->math_funcs->scalar_product;
    PyMultivectorObject *out;
    if(scalar_product){
        out = scalar_product(self,value);
    }else{
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
            PyErr_SetString(PyExc_TypeError,"objects must be an instance of gasparse.multivector");
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
			//								"Elements of the list of values must be flot or");
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

};

PyDoc_STRVAR(add_doc, "adds a bunch of multivectors.");
PyDoc_STRVAR(dual_doc, "dualizes the multivector.");
PyDoc_STRVAR(undual_doc, "undualizes the multivector.");
PyDoc_STRVAR(product_doc, "multiplies a bunch of multivectors.");
PyDoc_STRVAR(exponential_doc, "takes the exponential of multivectors.");
PyDoc_STRVAR(list_doc,
						 "Returns a list with each coefficient of the multivector.");
PyDoc_STRVAR(cast_doc, "Casts the multivector to the specified type");
PyDoc_STRVAR(grade_doc, "Returns the grade of a multivector");

PyMethodDef multivector_methods[] = {
		{"dual", (PyCFunction)multivector_dual, METH_NOARGS, dual_doc},
		{"undual", (PyCFunction)multivector_undual, METH_NOARGS, undual_doc},
		{"add", (PyCFunction)multivector_atomic_add, METH_VARARGS | METH_CLASS,
		 add_doc},
		{"geometric_product", (PyCFunction)multivector_atomic_geometric_product,
		 METH_VARARGS | METH_CLASS, product_doc},
		{"outer_product", (PyCFunction)multivector_atomic_outer_product,
		 METH_VARARGS | METH_CLASS, product_doc},
		{"exp", (PyCFunction)multivector_exponential, METH_VARARGS | METH_CLASS,
		 exponential_doc},
		//{"list", (PyCFunction)multivector_list, METH_NOARGS, list_doc},
		{"list", (PyCFunction)multivector_list, METH_VARARGS | METH_KEYWORDS,
		 list_doc},
		{"cast", (PyCFunction)multivector_cast, METH_VARARGS, cast_doc},
		{"grade", (PyCFunction)multivector_grade, METH_NOARGS, grade_doc},
		{NULL},
};

PyTypeObject PyMultivectorType = {
		PyVarObject_HEAD_INIT(NULL, 0).tp_name = "gasparse.multivector",
		.tp_doc = PyDoc_STR(
				"Builds a multivector in different types (sparse,dense,blades)"),
		.tp_basicsize = sizeof(PyMultivectorObject),
		.tp_itemsize = 0,
		.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
		.tp_dealloc = (destructor)multivector_array_dealloc,
		.tp_repr = (reprfunc)multivector_repr,
		.tp_str = (reprfunc)multivector_repr,
		.tp_call = (ternaryfunc)multivector_grade_project,
		.tp_new = NULL,
		.tp_as_number = &PyMultivectorNumberMethods,
		.tp_methods = multivector_methods};


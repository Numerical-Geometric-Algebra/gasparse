-----------------
Multivector Types
-----------------

One of the main goals of this project is being able to define new data structures for multivectors (the type of multivector) for which specialized operations can be dispatched for better performance. 
Different approaches to represent multivectors in a data structure can result in better or worse performance and can also enable different types of optimizations. 
There are multiple ways to represent multivectors, for example we can allways represent multivector using duality :math:`x=a+bI`, 
or a two blade as the wedge product of two vectors :math:`x = a\wedge b`. For :math:`x=a+bI` we would store :math:`a` and :math:`b` in direct space 
which can reduce the size of the cayley table making computations more efficient. For :math:`x = a\wedge b` we would store the vectors :math:`a` and :math:`b` in the data structure, 
in an `n` dimensional algebra the bivector basis yields :math:`{n\choose 2}` coefficients while two vectors need :math:`2n` basis elements. Note that :math:`{n\choose 2}\geqslant 2n` for :math:`n\geqslant 5`.

Depending on the algebra or on the data that is going to be handled we can chose from three different multivector types. The available multivector types are:

.. code-block:: c

    // Use dtype='sparse' to chose this type
    typedef struct SparseMultivector{
        int *bitmap;
        ga_float *value;
        Py_ssize_t size;
    }SparseMultivector;

    // Use dtype='blades' to chose this type
    typedef struct BladesMultivector{
        SparseMultivector *data;
        Py_ssize_t *grade;
        Py_ssize_t size;
    }BladesMultivector;

    // Use dtype='dense' to chose this type
    typedef struct DenseMultivector{
        ga_float *value;
        Py_ssize_t size; 
    }DenseMultivector;


``DenseMultivector`` are indexed by bitmaps. The ``BladesMultivector`` are indexed by grades where for each grade we have a ``SparseMultivector`` structure. 
The ``SparseMultivector`` represent multivectors by ``value`` and by its corresponding ``bitmap``, bitmaps indicate the basis multivector that corresponds to that value. 
The bit `i` indicates if there is a vector in the basis multivector. For example considering ``b'101101'`` the corresponding basis blade is :math:`e_{1346} = e_1\wedge e_3\wedge e_4\wedge e_6`. 
Using this representation for bitmaps makes it very easy to compute bitmaps resultant from the geometric product as ``b1 ^ b2`` where ``^`` is the bitwise xor operator. Yet signs are much more difficult 
to compute since it involves determining how many times we need to swap basis vectors to get to the canonical basis and determining the signs resulted from squaring vectors.

For the sparse multivectors we used an approach that allocates memory for the entire basis of the multivector and uses the indices as bitmaps. Thus when we get a bitmap as the result of an operation we 
simply take that as an index e.g. ``dense.bitmap[dense0->bitmap[i]] = dense0.bitmap[i]``. Then we ignore all the unvisited indices by iterating over the entire multivector and checking if the bitmap 
is equal to ``-1``.

.. note:: 

    Iterating over the entire geometric algebra is expensive. For an `n` dimensional geometric algebra we have to loop over :math:`2^n` values and bitmaps.


To avoid iterating over all basis elements of the geoemetric algebra we implemented a linked list and array of pointers strategy to (i.) 
add a new element to the list when a new basis is encountered (ii.) to index using bitmaps the pointers to each element of the linked list. 
That is we initialize an empty array ``addr`` of the size of the algebra. Then for each ``bitmap`` it will be assigned a pointer ``graph`` to an element of the linked list ``addr[bitmap]=graph``,
for each new basis element with associated bitmap. The main advantage of this approach is that we do not need to iterate over all elements of the geometric algebra.

.. Summarizing the "`algorithm`"

.. #. Allocating memory for a new element ``graph = graph_new(); graph->value = value; graph->bitmap = bitmap; addr[bitmap] = graph``
.. #. Repeated basis elements are summed using the address array ``addr[bitmap]->value += value``



.. To understand how we compute products for each of the different types of multivectors we show some functions for the different types


.. API
.. ---

.. We have defined some helper functions and macros. 

.. Index data from a multivector array using the macro `INDEX_DATA`

.. .. code-block:: c
    
..     #define INDEX_DATA(s,i) ((s)->data + (i)*(s)->type->basic_size)
    
.. Example usage (product of a multivector with a scalar multivector)

.. .. code-block:: c
    
..     out = new_mvarray_from_mvarray(data);
..     scalar_product = data->type->math_funcs->scalar_product;
..     for(Py_ssize_t i = 0; i < size; i++){
..         ScalarMultivector *scalar_mv = INDEX_DATA(scalar, i);
..         if(!scalar_product(INDEX_DATA(out, i),INDEX_DATA(data, i),data->GA,*scalar_mv)){
..             multivector_array_dealloc(out);
..             return NULL;
..         }
..     }

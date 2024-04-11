-----------------
Multivector Types
-----------------

One of the main goals of this project is being able to define new data structures for multivectors (the type of multivector) for which specialized operations can be dispatched for better performance. 
Different approaches to represent multivectors in a data structure can result in better or worse performance and can also enable different types of optimizations. 
There are multiple ways to represent multivectors, for example we can allways represent multivector using duality :math:`x=a+bI`, 
or a two blade as the wedge product of two vectors :math:`x = a\wedge b`. For :math:`x=a+bI` we would store :math:`a` and :math:`b` in direct space 
which can reduce the size of the cayley table making computations more efficient. For :math:`x = a\wedge b` we would store the vectors :math:`a` and :math:`b` in the data structure, 
in an `n` dimensional algebra the bivector basis yields :math:`{n\choose 2}` coefficients while two vectors need :math:`2n` basis elements. Note that :math:`{n\choose 2}\geqslant 2n` for :math:`n\geqslant 5`.


The class where multivector arrays reside is ``gasparse.mvarray`` but instances of this class can be of a subtype. The available subtypes are

- sparse (``SparseMultivector``) 
- dense (``DenseMultivector``)
- graded sparse (``BladesMultivector``)
- graded dense **NOT IMPLEMENTED YET**
- scalar (``ga_float``)


Bitmaps
=======

Bitmaps are integer values that are used to represent basis blades. We consider that the bitmaps represent basis blades in canonical order, that is increasing order from left to right. 
Each bit in a bitmap represents a basis vector, if the bit is zero then that basis vector is not present in the basis blade, if the bit is zero then it is present. 
A simple example is shown: the bitmap ``0b101101`` the corresponding basis blade is :math:`e_{1346} = e_1\wedge e_3\wedge e_4\wedge e_6`. 

Using this representation for bitmaps makes it very easy to compute bitmaps resultant from the geometric product as ``bitmap_left ^ bitmap_right`` where ``^`` is the bitwise xor operator. 
Yet signs are much more difficult to compute since it involves determining how many times we need to swap basis vectors to get to the canonical order and determining the signs resulted from squaring vectors.


The ``C`` types
===============

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

    // NOT IMPLEMENTED YET
    typedef struct GradedDenseMultivector{ 
        ga_float **values;
        Py_ssize_t *grades;
        Py_ssize_t size; // The number of different grades
    }GradedDenseMultivector;


``DenseMultivector`` are indexed by bitmaps. The ``BladesMultivector`` are indexed by grades where for each grade we have a ``SparseMultivector`` structure. 
The ``SparseMultivector`` represent multivectors by ``value`` and by its corresponding ``bitmap``, 
For the sparse multivectors we used an approach that allocates memory for the entire basis of the multivector and uses the indices as bitmaps. Thus when we get a bitmap as the result of an operation we 
simply take that as an index e.g. ``dense.bitmap[dense0->bitmap[i]] = dense0.bitmap[i]``. Then we ignore all the unvisited indices by iterating over the entire multivector and checking if the bitmap 
is equal to ``-1``.

.. note:: 

    Iterating over the entire geometric algebra is expensive. For an `n` dimensional geometric algebra we have to loop over :math:`2^n` values and bitmaps.

Computation Modes
=================

One of the main motivations for writting this package is that after describing the general "algorithm" to compute the different operations we can allways find new ways to improve the performance
of said "algorithms". There are several ways to implement operations between multivectors, several data structures where defined in order to store the operations `"information"`, included data structures are

- Table of positions and signs (only for ``DenseMultivector`` and ``GradedDenseMultivector``   **NOT IMPLEMENTED YET**)  `The Table of positions and signs`_
- A two dimensional array that maps bitmaps to signs and bitmaps (for ``SparseMultivector``) `The Bitmap Map`_

As for the different computation strategies we have

- Linked lists strategy. The output of operations is stored in a linked list. Pointers to nodes are stored in an array indexed by bitmaps. (``compute_mode='large'``) `Linked Lists`_.
- Entire algebra basis alocation. The output is directly stored on an array indexed by bitmaps. The output is then compressed to a sparse data structure. (``compute_mode='generic'``).
- Code generated strategy. Using the table of positions and signs code generate the different operations (``compute_mode='generated'``).

.. _`The bitmap Map`:

The bitmap Map
""""""""""""""

Signs and bitmaps of products are determined by indexing the sign array and the bitmap array of the map. 

.. code-block:: c
    
    sign_out = map.sign[bitmap_left][bitmap_right];
    bitmap_out = map.bitmap[bitmap_left][bitmap_right];
    value_out = sign_out*value_left*value_right;

.. _`The Table of positions and signs`:

The Table of positions and signs
""""""""""""""""""""""""""""""""

An array of a data structure containing positions and signs is created then the product is computed as a loop of the following

.. code-block:: c

    value_out[table[i].position[2]] = table[i].sign*value_left[table[i].position[0]]
                                      *value_right[table[i].position[1]];

By defining different values for the positions and signs of the table we can define the different types of products in this way. 

.. _`Linked Lists`:

Linked Lists
""""""""""""

For the ``"large"`` computation mode we consider using linked lists for storing the resulting output of the product between basis elements. 
We also define an array of pointers to items of the linked list, indexed by bitmaps. That is, when the result of taking the product between two basis elements is some bitmap ``bitmap`` 
we store the item ``item`` of the linked list at position ``bitmap`` in the array ``addr`` thus 

.. code-block:: c

    if(is_item_new(addr,bitmap)){
        item->next = item_alloc_memory();
        item = item->next;
        item->bitmap = bitmap;
        item->value = value_out;
        addr[bitmap] = item;
    }else
        addr[bitmap] += value_out;
    
the ``is_item_new`` function is simply ``!addr[bitmap]``. If at position bitmap ``addr`` points to null then it means that this item never apperead before.


By the extra cost of allocating memory, assigning the item in the adress array and checking if the item is new, this strategy avoids iterating over all the basis elements of a geometric algebra. 
Which large geometric algebras take huge benifits from.


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
